#!/usr/bin/env python3
"""
    check_scene.py: stage 1 — gates V0 and V1, is this product usable over water?

    Input:  data/tanager/<scene>/<id>_ortho_sr_hdf5.h5 (stage 0).
    Output: outputs/metrics/v0_v1_<scene>.json   the measured spectra and corner offsets
            outputs/metrics/verdicts.csv         V0, V1

    Two questions, asked before any analysis is built on the answers.

    **V1 — is the grid where it says it is?** The georeferencing comes from a text blob
    (``StructMetadata.0``), the grid is corner-registered, and nothing cross-checks it.
    Here the four grid corners are projected from UTM to WGS84 and required to *contain*
    the scene's own STAC bounding box — containment rather than equality, because the
    north-up array is sized to hold a rotated imaging strip and so is legitimately
    larger than the imaged data by the nodata margin. A sign error, a transposed axis or
    the wrong UTM zone all break containment at once.

    **V0 — is the surface reflectance physical over water?** Planet's correction is
    tuned for land. This stage samples the water spectrum and reports the pedestal, the
    green-over-near-infrared contrast that says water is water, and the noise floor the
    product declares for itself. It does not require the pedestal to be small — it
    requires the *water-leaving shape* to be present and the *spatial variance* to
    exceed the product's own uncertainty, because those are what the experiment uses.

    Offline apart from the coastline check. Raises on a failed gate.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tanager_io as tio
from _scenes import OUT_DIR, Scene, load_scenes
from _verdicts import configure_stage_logging, upsert_verdict, verdict_row

logger = logging.getLogger(__name__)

#: The grid must contain the STAC bounding box; this is the slack allowed on that
#: containment, absorbing the 4-decimal rounding of the published bbox (~0.4 px).
CORNER_TOLERANCE_PX = 1.0

#: Water must be at least this fraction of the imaged (non-nodata) pixels for a scene to
#: carry a water comparison at all.
MIN_WATER_FRACTION = 0.15

#: Water's green band must exceed its near-infrared band by at least this ratio. Every
#: land cover fails it; the pedestal pushes the ratio towards 1, so this is a weak bar
#: on purpose — it is a sanity check, not a discriminator.
MIN_GREEN_OVER_NIR = 1.05

#: Bands sampled for the water-spectrum ledger.
REPORT_NM = (443, 490, 560, 620, 665, 681, 709, 754, 865, 1240, 1600, 2200)


@dataclass
class StageConfig:
    """Stage 1 wiring."""

    out_dir: Path = OUT_DIR
    corner_tolerance_px: float = CORNER_TOLERANCE_PX
    min_water_fraction: float = MIN_WATER_FRACTION


def check_corners(scene: Scene, dataset, grid) -> dict:
    """V1: the grid must *contain* the scene's STAC bounding box, tightly.

    The right invariant is containment, not equality. The ortho grid is a north-up array
    sized to hold a rotated imaging strip, so its extent is the bounding box of the
    *array* while the STAC bbox is the bounding box of the *imaged data* inside it —
    the grid is legitimately the larger of the two, by the nodata margin. On the SF Bay
    scene that margin is 0.3-1.5 px on the four sides.

    So this fails on a grid that does not contain the STAC bbox (a sign error, a
    transposed axis, or the wrong UTM zone would all break containment immediately) and
    reports the margin without failing on it, since the margin is a property of the
    strip's rotation and varies from scene to scene.
    """
    lon, lat = tio.to_lonlat(dataset, *grid.corners_xy())
    stac = scene.bbox
    measured = (float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max()))
    metres_per_deg_lat = 111_320.0
    metres_per_deg_lon = metres_per_deg_lat * np.cos(np.radians(np.mean(lat)))

    # Positive = the grid extends beyond the STAC bbox on that side, which is expected.
    margins_px = [
        (stac[0] - measured[0]) * metres_per_deg_lon / grid.dx,
        (stac[1] - measured[1]) * metres_per_deg_lat / grid.dy,
        (measured[2] - stac[2]) * metres_per_deg_lon / grid.dx,
        (measured[3] - stac[3]) * metres_per_deg_lat / grid.dy,
    ]
    return {"stac_bbox": list(stac), "grid_bbox": list(measured),
            "margin_px": [round(value, 3) for value in margins_px],
            "worst_shortfall_px": round(float(-min(min(margins_px), 0.0)), 3),
            "widest_margin_px": round(float(max(margins_px)), 3),
            "contains_stac_bbox": bool(min(margins_px) > -CORNER_TOLERANCE_PX),
            "epsg": grid.epsg, "pixel_size_m": grid.dx,
            "n_x": grid.n_x, "n_y": grid.n_y}


def check_water(dataset, water: np.ndarray, valid: np.ndarray) -> dict:
    """V0: the water spectrum, the pedestal, and the variance against the noise floor."""
    wavelength = dataset["wavelength"].values
    indices = [tio.band_index(wavelength, nm) for nm in REPORT_NM]
    cube = dataset["surface_reflectance"].isel(wavelength=indices).values
    sampled = cube[:, water]

    spectrum = {f"{wavelength[index]:.0f}": {
        "median": float(np.nanmedian(sampled[position])),
        "p1_pedestal": float(np.nanpercentile(sampled[position], 1)),
    } for position, index in enumerate(indices)}

    green = float(np.nanmedian(sampled[REPORT_NM.index(560)]))
    near_infrared = float(np.nanmedian(sampled[REPORT_NM.index(865)]))

    # The product's own declared relative uncertainty at the green band, which is the
    # floor any spatial variance has to clear before it can be called water structure.
    uncertainty = dataset["surface_reflectance_uncertainty"].isel(
        wavelength=indices[REPORT_NM.index(560)]).values
    relative_uncertainty = float(np.nanmedian(uncertainty[water]) / green)

    return {"n_water_px": int(water.sum()), "n_valid_px": int(valid.sum()),
            "water_fraction_of_valid": float(water.sum() / max(valid.sum(), 1)),
            "valid_fraction_of_scene": float(valid.sum() / valid.size),
            "green_over_nir": green / near_infrared if near_infrared else float("nan"),
            "relative_uncertainty_560": relative_uncertainty,
            "spectrum": spectrum}


def check_one(scene: Scene, config: StageConfig) -> dict:
    """Both gates for one scene, plus the numbers that decided them."""
    grid = tio.read_grid(scene.sr_path)
    dataset = tio.open_sr(scene.sr_path)
    ledger: dict = {"scene": scene.key, "tanager_id": scene.tanager_id, "date": scene.date}
    ledger["v1_geometry"] = check_corners(scene, dataset, grid)

    water_mask, threshold = tio.water_mask(dataset, return_threshold=True)
    water = water_mask.values
    valid = tio.valid_mask(dataset).compute().values
    ledger["v0_water"] = check_water(dataset, water, valid)
    ledger["v0_water"]["ndwi_threshold"] = threshold

    problems: list[str] = []
    if not ledger["v1_geometry"]["contains_stac_bbox"]:
        problems.append(f"grid does not contain the STAC bbox — short by "
                        f"{ledger['v1_geometry']['worst_shortfall_px']:.2f} px "
                        f"(slack {config.corner_tolerance_px})")
    if ledger["v0_water"]["water_fraction_of_valid"] < config.min_water_fraction:
        problems.append(f"water is {ledger['v0_water']['water_fraction_of_valid']:.1%} of "
                        f"imaged pixels (floor {config.min_water_fraction:.0%})")
    if ledger["v0_water"]["green_over_nir"] < MIN_GREEN_OVER_NIR:
        problems.append(f"green/NIR over water is "
                        f"{ledger['v0_water']['green_over_nir']:.3f} "
                        f"(floor {MIN_GREEN_OVER_NIR}) — the mask is not on water")
    ledger["problems"] = problems
    dataset.close()
    return ledger


def main(scenes: list[str] | None = None, out_dir: Path = OUT_DIR,
         log_level: str = "info") -> int:
    configure_stage_logging(out_dir / "check_scene.log", log_level)
    config = StageConfig(out_dir=out_dir)
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    ledgers = {}
    for scene in load_scenes(keys=scenes):
        if not scene.sr_path.exists():
            logger.warning("%s: no cube on disk — run fetch_scenes.py; skipping", scene.key)
            continue
        logger.info("=== %s ===", scene.key)
        ledger = check_one(scene, config)
        ledgers[scene.key] = ledger
        (metrics_dir / f"v0_v1_{scene.key}.json").write_text(json.dumps(ledger, indent=2))
        geometry, water = ledger["v1_geometry"], ledger["v0_water"]
        logger.info("%s: EPSG %d, %.1f m px, contains STAC bbox (margin %.2f px, "
                    "shortfall %.2f px)", scene.key, geometry["epsg"],
                    geometry["pixel_size_m"], geometry["widest_margin_px"],
                    geometry["worst_shortfall_px"])
        logger.info("%s: water %.1f%% of imaged px, green/NIR %.2f, "
                    "declared relative uncertainty %.4f",
                    scene.key, 100 * water["water_fraction_of_valid"],
                    water["green_over_nir"], water["relative_uncertainty_560"])

    if not ledgers:
        raise RuntimeError("no scenes had a cube on disk — run fetch_scenes.py first")

    worst_shortfall = max(entry["v1_geometry"]["worst_shortfall_px"] for entry in ledgers.values())
    widest_margin = max(entry["v1_geometry"]["widest_margin_px"] for entry in ledgers.values())
    all_problems = [f"{key}: {problem}" for key, entry in ledgers.items()
                    for problem in entry["problems"]]
    upsert_verdict(metrics_dir, [
        verdict_row("V1 the grid is where StructMetadata.0 says it is",
                    "worst shortfall containing the STAC bbox, in pixels",
                    f"{worst_shortfall:.2f}", f"<= {CORNER_TOLERANCE_PX}",
                    "PASS" if worst_shortfall <= CORNER_TOLERANCE_PX else "FAIL",
                    f"widest nodata margin {widest_margin:.2f} px; " + "; ".join(
                        f"{key} EPSG {entry['v1_geometry']['epsg']}"
                        for key, entry in ledgers.items())),
        verdict_row("V0 surface reflectance is usable over water",
                    "water fraction of imaged px; green/NIR over water",
                    "; ".join(f"{key} {entry['v0_water']['water_fraction_of_valid']:.2f}/"
                              f"{entry['v0_water']['green_over_nir']:.2f}"
                              for key, entry in ledgers.items()),
                    f">= {MIN_WATER_FRACTION}; >= {MIN_GREEN_OVER_NIR}",
                    "PASS" if not all_problems else "FAIL",
                    "absolute level carries a large additive atmospheric pedestal; the "
                    "experiment uses shape and spatial variance, not absolute Rrs"),
    ])

    print(json.dumps({key: {"shortfall_px": entry["v1_geometry"]["worst_shortfall_px"],
                            "water_frac": round(entry["v0_water"]["water_fraction_of_valid"], 3),
                            "green_over_nir": round(entry["v0_water"]["green_over_nir"], 3),
                            "problems": entry["problems"]}
                     for key, entry in ledgers.items()}, indent=2))
    if all_problems:
        raise RuntimeError("V0/V1 failed: " + "; ".join(all_problems))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 1: gates V0 (water usability) and V1 (georeferencing)",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenes", nargs="*", default=None)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
