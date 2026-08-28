#!/usr/bin/env python3
"""
    aggregate.py: exp 35 stage 2 — bin 30 m Tanager water pixels into coarse sensor pixels.

    Input:  data/tanager/<scene>/<id>_ortho_sr_hdf5.h5, PACE_OCI.*.L2.OC_{AOP,BGC}.*.nc,
            olci_l3.nc  (stage 0); gates V0/V1 (stage 1).
    Output: outputs/<scene>_<sensor>_agg.parquet   one row per coarse pixel
            outputs/metrics/g2_<scene>.json        counts and the V2 check
            outputs/metrics/verdicts.csv           V2

    One row per coarse (PACE or OLCI) pixel, carrying:

    * how many Tanager water pixels fell inside it, and what fraction of it they cover;
    * the mean and standard deviation of every Tanager band in the ocean-color range;
    * the mean, standard deviation, and extremes of each *derived 30 m scalar* — the
      band ratios. These cannot be recovered from per-band moments, because the mean of
      a ratio is not the ratio of the means, and it is exactly this difference the
      experiment is about;
    * the coarse sensor's own reflectances and chlorophyll for the same pixel.

    **V2, the linchpin gate.** Tanager's variance is only credible if it agrees with an
    independent sensor at a scale both can see. Tanager is therefore aggregated to the
    OLCI 300 m grid and the *between-pixel* variability of the two is compared. If they
    agree, the within-pixel variance at 30 m — which OLCI cannot see at all — is
    licensed. If they do not, nothing downstream is.

    That comparison is two numbers, not one: a correlation, which says the two see the
    same pattern, and a regression slope, which says whether they see the same *amount*
    of it. Correlation alone is scale-invariant and would pass a Tanager field with ten
    times OLCI's spread — see :func:`check_v2`.

    Assignment is nearest-center (a Voronoi partition) via a KD-tree in the scene's UTM
    meters, rather than building and rasterising footprint polygons: for a smooth swath
    the two are equivalent, and this one has no topology to get wrong.

    Offline. ~3 min per scene, dominated by reading ~100 bands off disk.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tanager_io as tio
import xarray as xr
from _scenes import OUT_DIR, Scene, load_scenes
from _verdicts import configure_stage_logging, upsert_verdict, verdict_row
from pyproj import Transformer
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)

#: Tanager bands whose per-pixel moments are kept. The ocean-color range; beyond ~900 nm
#: water is black and the bands carry only the atmosphere.
BAND_RANGE_NM = (400.0, 900.0)

#: Derived 30 m scalars, as (name, numerator nm, denominator nm). Ratios, because a ratio
#: divides out a multiplicative error and — once the additive pedestal is removed — is
#: what carries the water's optical signal.
RATIOS = (("blue_green", 490.0, 560.0),
          ("green_blue", 560.0, 490.0),
          ("red_green", 665.0, 560.0),
          ("green_red", 560.0, 665.0))

#: A coarse pixel is used only if this fraction of it is Tanager water.
MIN_WATER_FRACTION = 0.90

#: V2 passes if Tanager-aggregated and OLCI's own between-pixel variability correlate
#: at least this well across the scene.
V2_MIN_CORRELATION = 0.60

#: Both sides of V2 are now remote-sensing reflectance in sr^-1 — Tanager's via
#: :func:`tanager_io.to_rrs`, which divides Planet's dimensionless surface reflectance
#: rho by pi (rho = pi * Rrs) after removing the additive pedestal. Parity is therefore a
#: regression slope of 1. Until 2026-08-25 this stage subtracted the pedestal inline and
#: never divided, so it was regressing rho on Rrs and parity was a slope of pi; nothing
#: downstream noticed, because every statistic below is a CV or a correlation and both
#: are scale-free. That is exactly why the pipeline now goes through the library.
V2_PARITY_SLOPE = 1.0

#: The bands the primary 560/490 ratio is built from, and the only ones V2 gates on.
#: 443 nm is measured and reported because it is the weakest band and hiding that would
#: misrepresent the product — but nothing downstream is computed from it.
V2_GATED_BANDS = ("490", "560")

#: V2 also requires Tanager's between-pixel spread to sit within this factor of OLCI's
#: once the sr^-1 difference is divided out. **Added on revision**, not pre-registered:
#: the original gate was a correlation, which is scale-invariant and therefore cannot
#: license a claim about variance — and every claim downstream of V2 is about variance.
V2_MAX_SPREAD_FACTOR = 2.5

#: Read this many Tanager bands at a time. Matches the file's own 14-band chunking.
BAND_BLOCK = 14


@dataclass
class StageConfig:
    """Stage 2 wiring."""

    out_dir: Path = OUT_DIR
    band_range_nm: tuple[float, float] = BAND_RANGE_NM
    min_water_fraction: float = MIN_WATER_FRACTION
    sensors: tuple[str, ...] = ("olci", "pace")
    ratios: tuple[tuple[str, float, float], ...] = field(default_factory=lambda: RATIOS)


def coarse_centres(scene: Scene, sensor: str) -> tuple[np.ndarray, np.ndarray, xr.Dataset]:
    """Longitude and latitude of every coarse pixel center, flattened, plus its dataset.

    OLCI arrives as a regular lat/lon grid and PACE as a swath with two-dimensional
    navigation arrays; both are reduced to flat center vectors so the caller does not
    care which it has.
    """
    if sensor == "olci":
        dataset = xr.open_dataset(scene.olci_path)
        if "time" in dataset.dims:
            dataset = dataset.isel(time=0)
        lon2d, lat2d = np.meshgrid(dataset["lon"].values, dataset["lat"].values)
        return lon2d.ravel(), lat2d.ravel(), dataset

    navigation = xr.open_dataset(scene.pace_aop_path, group="navigation_data")
    geophysical = xr.open_dataset(scene.pace_aop_path, group="geophysical_data")
    biogeochemistry = xr.open_dataset(scene.pace_bgc_path, group="geophysical_data")
    dataset = xr.merge([geophysical, biogeochemistry], compat="override", join="override")
    dataset = dataset.assign_coords(latitude=navigation["latitude"],
                                    longitude=navigation["longitude"])
    return (navigation["longitude"].values.ravel(),
            navigation["latitude"].values.ravel(), dataset)


def assign_labels(dataset: xr.Dataset, water: np.ndarray, lon: np.ndarray, lat: np.ndarray,
                  ) -> tuple[np.ndarray, np.ndarray, int]:
    """Nearest coarse-pixel index for every Tanager water pixel.

    Returns ``(labels, keep, n_coarse)`` where ``labels`` indexes the flattened coarse
    grid, ``keep`` selects the Tanager water pixels that fell within one coarse pixel
    spacing of a center (the rest are outside the swath), and ``n_coarse`` is the size of
    the flattened coarse grid.
    """
    to_utm = Transformer.from_crs("EPSG:4326", dataset.attrs["crs"], always_xy=True)
    finite = np.isfinite(lon) & np.isfinite(lat)
    coarse_x, coarse_y = to_utm.transform(lon[finite], lat[finite])
    tree = cKDTree(np.column_stack([coarse_x, coarse_y]))

    grid_x, grid_y = np.meshgrid(dataset["x"].values, dataset["y"].values)
    distance, nearest = tree.query(np.column_stack([grid_x[water], grid_y[water]]), k=1)

    # A Tanager pixel further from every coarse center than the coarse spacing is not
    # inside any coarse pixel — it is off the edge of the swath.
    spacing = float(np.median(tree.query(np.column_stack([coarse_x, coarse_y]), k=2)[0][:, 1]))
    keep = distance <= spacing
    labels = np.flatnonzero(finite)[nearest[keep]]
    return labels, keep, lon.size


def moments(labels: np.ndarray, values: np.ndarray, n_bins: int) -> tuple[np.ndarray, ...]:
    """Count, mean and standard deviation of ``values`` grouped by ``labels``.

    ``np.bincount`` over the label array is what makes this affordable: one pass per
    band, constant memory in the number of bands, no grouping object built.
    """
    finite = np.isfinite(values)
    label = labels[finite]
    value = values[finite]
    count = np.bincount(label, minlength=n_bins).astype("float64")
    total = np.bincount(label, weights=value, minlength=n_bins)
    total_square = np.bincount(label, weights=value * value, minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
        variance = np.maximum(total_square / count - mean * mean, 0.0)
    return count, mean, np.sqrt(variance)


def aggregate_scene(scene: Scene, sensor: str, config: StageConfig) -> tuple[pd.DataFrame, dict]:
    """One row per coarse pixel for one scene and one coarse sensor."""
    dataset = tio.open_sr(scene.sr_path)
    water_mask = tio.water_mask(dataset)
    water = water_mask.values
    valid = tio.valid_mask(dataset).compute().values

    lon, lat, coarse = coarse_centres(scene, sensor)
    labels, keep, n_coarse = assign_labels(dataset, water, lon, lat)
    logger.info("%s/%s: %d water px -> %d assigned, %d coarse pixels on the grid",
                scene.key, sensor, int(water.sum()), int(keep.sum()), n_coarse)

    # How much of each coarse pixel the Tanager scene actually covers: water pixels over
    # every Tanager pixel (water or not) that the same coarse pixel claims.
    all_labels, all_keep, _ = assign_labels(dataset, valid, lon, lat)
    n_valid = np.bincount(all_labels, minlength=n_coarse).astype("float64")
    n_water = np.bincount(labels, minlength=n_coarse).astype("float64")
    with np.errstate(invalid="ignore", divide="ignore"):
        water_fraction = np.where(n_valid > 0, n_water / n_valid, np.nan)

    frame = pd.DataFrame({"coarse_index": np.arange(n_coarse), "lon": lon, "lat": lat,
                          "n_water": n_water, "n_valid": n_valid,
                          "frac_water": water_fraction})

    wavelength = dataset["wavelength"].values
    in_range = np.flatnonzero((wavelength >= config.band_range_nm[0])
                              & (wavelength <= config.band_range_nm[1]))
    logger.info("%s/%s: %d bands in %.0f-%.0f nm", scene.key, sensor, in_range.size,
                *config.band_range_nm)
    for start in range(0, in_range.size, BAND_BLOCK):
        block = in_range[start:start + BAND_BLOCK]
        # A block at a time, through the library. Both helpers act on whatever Dataset
        # they are handed, so subsetting the wavelength axis keeps the 426-band cube out
        # of memory without reimplementing either of them here — which is how the
        # pipeline came to be comparing surface reflectance against Rrs in the first
        # place. `.load()` reads the block once for both calls.
        subset = dataset[["surface_reflectance"]].isel(wavelength=block).load()
        rrs = tio.to_rrs(subset, offset=tio.dark_spectrum_offset(subset, water_mask)).values
        for position, band_index in enumerate(block):
            values = rrs[position][water][keep]
            _, mean, deviation = moments(labels, values, n_coarse)
            name = f"{wavelength[band_index]:.0f}"
            frame[f"rrs_mean_{name}"] = mean
            frame[f"rrs_std_{name}"] = deviation

    # The derived scalars, computed per 30 m pixel *before* aggregation. The pedestal is
    # estimated per band, with no coupling between bands, so taking it from a two-band
    # subset gives exactly what the whole-cube call gives for those two bands.
    for name, numerator_nm, denominator_nm in config.ratios:
        pair = dataset[["surface_reflectance"]].isel(wavelength=[
            tio.band_index(wavelength, numerator_nm),
            tio.band_index(wavelength, denominator_nm)]).load()
        ratio_field = tio.band_ratio(pair, numerator_nm, denominator_nm,
                                     source=tio.to_rrs(pair, water=water_mask))
        ratio = ratio_field.values[water][keep]
        count, mean, deviation = moments(labels, ratio, n_coarse)
        with np.errstate(invalid="ignore", divide="ignore"):
            frame[f"{name}_mean"] = mean
            frame[f"{name}_std"] = deviation
            frame[f"{name}_cv"] = np.where(mean != 0, deviation / np.abs(mean), np.nan)
        frame[f"{name}_n"] = count

    frame = attach_coarse_values(frame, coarse, sensor, n_coarse)
    usable = frame["frac_water"] >= config.min_water_fraction
    # The coarse grid's 2-D shape, recorded so stage 3 can reshape a flat index without
    # reopening the source granule. A PACE L2 file is ~200 MB and stage 3 wants nothing
    # else from it, which made the aggregates unusable on their own.
    grid_shape = ((coarse.sizes["lat"], coarse.sizes["lon"]) if sensor == "olci"
                  else tuple(int(size) for size in coarse["latitude"].shape))
    ledger = {"scene": scene.key, "sensor": sensor,
              "coarse_shape": [int(size) for size in grid_shape],
              "n_coarse_total": int(n_coarse),
              "n_coarse_with_water": int((frame["n_water"] > 0).sum()),
              "n_coarse_usable": int(usable.sum()),
              "min_water_fraction": config.min_water_fraction,
              "tanager_water_px": int(water.sum()),
              "tanager_px_assigned": int(keep.sum()),
              "median_px_per_coarse": float(frame.loc[usable, "n_water"].median())
              if usable.any() else 0.0}
    dataset.close()
    return frame.loc[frame["n_water"] > 0].reset_index(drop=True), ledger


def attach_coarse_values(frame: pd.DataFrame, coarse: xr.Dataset, sensor: str,
                         n_coarse: int) -> pd.DataFrame:
    """Add the coarse sensor's own retrievals on the same rows."""
    wanted = ["CHL"] if sensor == "olci" else ["chlor_a"]
    for name, variable in coarse.data_vars.items():
        if variable.size != n_coarse:
            continue
        if name.startswith("RRS") or name.startswith("Rrs") or name in wanted:
            frame[f"coarse_{name}"] = np.asarray(variable.values).ravel()
    return frame


def check_v2(frame: pd.DataFrame, sensor: str) -> dict:
    """V2: does Tanager's variability agree with OLCI's, at the scale OLCI can see?

    Tanager's per-coarse-pixel mean is its estimate of what OLCI measured. If the two
    disagree *spatially* across the scene then Tanager's variance is not the water's,
    and no within-pixel claim can be built on it.

    Two statistics, because one of them is not enough. The **correlation** says the two
    instruments see the same pattern. It is scale-invariant, so it would read 1.0 on a
    Tanager field with ten times OLCI's spread — and the whole experiment is a claim
    about spread. So the **regression slope** is reported alongside it. Both sides are
    Rrs in sr^-1 (see :data:`V2_PARITY_SLOPE`), so agreement means a slope of 1 and the
    slope is directly the factor by which Tanager overstates OLCI's 300 m variability.

    Every band is reported. The gate binds only on :data:`V2_GATED_BANDS`, the pair the
    primary 560/490 ratio is built from.
    """
    if sensor != "olci":
        return {"applicable": False}
    coarse_columns = [name for name in frame.columns if name.startswith("coarse_RRS")]
    if not coarse_columns:
        return {"applicable": False, "note": "no OLCI reflectance columns"}

    usable = frame[frame["frac_water"] >= MIN_WATER_FRACTION]
    results = {}
    for column in ("coarse_RRS443", "coarse_RRS490", "coarse_RRS560"):
        if column not in usable:
            continue
        nanometres = column.replace("coarse_RRS", "")
        tanager_column = next((name for name in usable.columns
                               if name.startswith("rrs_mean_")
                               and abs(float(name.split("_")[-1]) - float(nanometres)) <= 3), None)
        if tanager_column is None:
            continue
        pair = usable[[tanager_column, column]].dropna()
        if len(pair) < 20:
            continue
        tanager_values = pair[tanager_column].to_numpy()
        coarse_values = pair[column].to_numpy()
        slope = float(np.polyfit(coarse_values, tanager_values, 1)[0])
        results[nanometres] = {
            "n": int(len(pair)),
            "correlation": float(pair[tanager_column].corr(pair[column])),
            "slope_tanager_on_olci": slope,
            "spread_factor": slope / V2_PARITY_SLOPE,
            "sd_tanager": float(tanager_values.std()),
            "sd_olci": float(coarse_values.std()),
        }
    gated = {nm: entry for nm, entry in results.items() if nm in V2_GATED_BANDS}
    worst_correlation = min((entry["correlation"] for entry in gated.values()),
                            default=float("nan"))
    # Furthest from parity in log space, reported with its sign intact: a Tanager field
    # with *half* OLCI's spread is as much of a failure as one with twice, and averaging
    # or taking a bare maximum would hide which way it went.
    spreads = [entry["spread_factor"] for entry in gated.values() if entry["spread_factor"] > 0]
    worst_spread = max(spreads, key=lambda factor: abs(np.log(factor))) if spreads else float("nan")
    return {"applicable": True, "per_band": results,
            "gated_bands": list(V2_GATED_BANDS),
            "worst_gated_correlation": worst_correlation,
            "worst_gated_spread_factor": worst_spread,
            "n_usable_pixels": int(len(usable))}


def main(scenes: list[str] | None = None, sensors: list[str] | None = None,
         out_dir: Path = OUT_DIR, log_level: str = "info",
         verify_only_gate: bool = False) -> int:
    configure_stage_logging(out_dir / "aggregate.log", log_level)
    config = StageConfig(out_dir=out_dir,
                         sensors=tuple(sensors) if sensors else ("olci", "pace"))
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    if verify_only_gate:
        return verify_only(scenes, config.sensors, out_dir)

    ledgers: dict[str, dict] = {}
    for scene in load_scenes(keys=scenes):
        if not scene.sr_path.exists():
            logger.warning("%s: no cube on disk — skipping", scene.key)
            continue
        scene_ledger: dict = {}
        for sensor in config.sensors:
            path = scene.olci_path if sensor == "olci" else scene.pace_aop_path
            if not path.exists():
                logger.warning("%s/%s: %s missing — skipping", scene.key, sensor, path.name)
                continue
            logger.info("=== %s / %s ===", scene.key, sensor)
            frame, entry = aggregate_scene(scene, sensor, config)
            frame.to_parquet(out_dir / f"{scene.key}_{sensor}_agg.parquet", index=False)
            entry["v2"] = check_v2(frame, sensor)
            scene_ledger[sensor] = entry
            logger.info("%s/%s: %d usable coarse px (>=%.0f%% water), median %.0f Tanager "
                        "px each", scene.key, sensor, entry["n_coarse_usable"],
                        100 * config.min_water_fraction, entry["median_px_per_coarse"])
        ledgers[scene.key] = scene_ledger
        (metrics_dir / f"g2_{scene.key}.json").write_text(json.dumps(scene_ledger, indent=2))

    write_v2_verdict(ledgers, metrics_dir)

    print(json.dumps({key: {sensor: {"usable": entry["n_coarse_usable"],
                                     "median_px": entry["median_px_per_coarse"],
                                     "v2": entry.get("v2", {}).get("worst_gated_correlation")}
                            for sensor, entry in scene_entry.items()}
                      for key, scene_entry in ledgers.items()}, indent=2, default=str))
    return 0


def write_v2_verdict(ledgers: dict, metrics_dir: Path) -> None:
    """Gate and record V2 on both statistics, across every scene.

    The gate binds on the *worst* scene, not the best: one well-behaved scene cannot
    license the variance claim on the others.
    """
    assessed = {key: entry["olci"]["v2"] for key, entry in ledgers.items()
                if "olci" in entry and entry["olci"]["v2"].get("applicable")}
    correlations = [entry["worst_gated_correlation"] for entry in assessed.values()
                    if np.isfinite(entry.get("worst_gated_correlation", np.nan))]
    spreads = [entry["worst_gated_spread_factor"] for entry in assessed.values()
               if np.isfinite(entry.get("worst_gated_spread_factor", np.nan))]
    worst_correlation = min(correlations) if correlations else float("nan")
    worst_spread = (max(spreads, key=lambda factor: abs(np.log(factor)))
                    if spreads else float("nan"))
    correlation_ok = bool(correlations) and worst_correlation >= V2_MIN_CORRELATION
    spread_ok = bool(spreads) and abs(np.log(worst_spread)) <= np.log(V2_MAX_SPREAD_FACTOR)
    passed = correlation_ok and spread_ok

    per_band = "; ".join(
        f"{key} " + ", ".join(f"{nm} r={band['correlation']:.2f} "
                              f"x{band['spread_factor']:.2f}"
                              for nm, band in sorted(entry["per_band"].items()))
        for key, entry in assessed.items())
    upsert_verdict(metrics_dir, [verdict_row(
        "V2 Tanager reproduces OLCI's variability at the scale OLCI can see",
        "per band: correlation, and regression slope (Tanager's spread as a factor of "
        "OLCI's, both in sr-1) — worst gated band, worst scene",
        f"r {worst_correlation:.3f}; spread x{worst_spread:.2f}"
        if correlations else "not computed",
        f"r >= {V2_MIN_CORRELATION}; spread within {V2_MAX_SPREAD_FACTOR}x of parity",
        "PASS" if passed else ("FAIL" if correlations else "NOT ASSESSABLE"),
        (f"{per_band}. Gated on {'/'.join(V2_GATED_BANDS)} nm, the bands the 560/490 "
         "ratio is built from; 443 nm is reported but not gated. Planet ships surface "
         "reflectance, which tanager_io.to_rrs divides by pi, so both sides are Rrs in "
         "sr-1 and parity is a slope of 1. The spread statistic was added on revision: "
         "correlation is scale-invariant and cannot license a claim about variance."))])
    if correlations and not passed:
        raise RuntimeError(
            f"V2 failed: worst gated correlation {worst_correlation:.3f} "
            f"(needs >= {V2_MIN_CORRELATION}), worst spread factor {worst_spread:.2f} "
            f"(needs within {V2_MAX_SPREAD_FACTOR}x of parity)")


def verify_only(scenes: list[str] | None, sensors: tuple[str, ...], out_dir: Path) -> int:
    """Re-gate V2 from the aggregates already on disk, without re-binning anything.

    Stage 2's binning is the expensive part (~3 min per scene and sensor) and its output
    is deterministic, so a change to the *gate* should not have to pay for it again.
    """
    metrics_dir = out_dir / "metrics"
    ledgers: dict = {}
    for scene in load_scenes(keys=scenes):
        path = metrics_dir / f"g2_{scene.key}.json"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing — run exp35-aggregate without --verify-only")
        scene_ledger = json.loads(path.read_text())
        for sensor in sensors:
            aggregate_path = out_dir / f"{scene.key}_{sensor}_agg.parquet"
            if sensor not in scene_ledger or not aggregate_path.exists():
                continue
            scene_ledger[sensor]["v2"] = check_v2(pd.read_parquet(aggregate_path), sensor)
        ledgers[scene.key] = scene_ledger
        path.write_text(json.dumps(scene_ledger, indent=2))
        logger.info("%s: V2 re-gated from %s", scene.key, aggregate_path.name)
    if not ledgers:
        raise RuntimeError("no scene ledgers found — run exp35-aggregate first")
    write_v2_verdict(ledgers, metrics_dir)
    print(json.dumps({key: entry.get("olci", {}).get("v2", {}).get("per_band", {})
                      for key, entry in ledgers.items()}, indent=2, default=str))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exp 35 stage 2: bin Tanager water pixels into PACE/OLCI footprints",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenes", nargs="*", default=None)
    parser.add_argument("--sensors", nargs="*", default=None, choices=["olci", "pace"])
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    parser.add_argument("--verify-only-gate", action="store_true",
                        help="re-gate V2 from the aggregates already on disk, without "
                             "re-binning (the binning is deterministic and slow)")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
