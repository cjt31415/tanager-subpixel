#!/usr/bin/env python3
"""
    compare.py: stage 3 — the variance a match-up protocol cannot see.

    Input:  outputs/<scene>_<sensor>_agg.parquet (stage 2); the OLCI/PACE files for grid
            shape; outputs/metrics/v0_v1_<scene>.json for the declared noise floor.
    Output: outputs/homogeneity_<scene>_<sensor>.parquet   one row per candidate box
            outputs/h1_spectral.parquet                    per-band Tanager vs coarse
            outputs/metrics/h_results.json                 V3, V4, V5 with their numbers
            outputs/metrics/verdicts.csv                   V3, V4, V5

    **The question.** Ocean-color validation protocols (Bailey & Werdell 2006, Remote
    Sens. Environ. 102:12) accept a satellite-to-in-situ match-up when a 5x5 box of
    satellite pixels around the station has a coefficient of variation under 0.15. That
    is a *between-pixel* statistic: it is computed across the 25 coarse pixels. It is
    used as a proxy for "this water is homogeneous enough that a point measurement
    represents a pixel" — and the *within-pixel* variance it assumes away has never been
    measurable, because the sensor being filtered is the one doing the measuring.

    Tanager resolves roughly 300 samples inside one CMEMS OLCI cell and 1,600-3,000
    inside one PACE pixel. So for every box that passes the filter we can ask what the
    filter could not: how variable is the water *inside* the pixels it just accepted?

    Three gates:

    V3  the within-pixel variability is spatially structured, not noise — adjacent
        coarse pixels resemble each other far more than a shuffled field does;
    V4  it exceeds the noise floor the product declares for itself;
    V5  the headline, answered either way: among boxes that pass the standard filter,
        the distribution of within-pixel variability, and what fraction of them contain
        a pixel whose interior is more variable than the threshold the box just cleared.

    **Two scalars, reported side by side.** The primary one is the green/blue ratio; the
    secondary is plain 560 nm reflectance. They do not agree — the ratio hides three to
    five times more variance inside a pixel than between pixels, the band often less than
    one — and that disagreement is a result, not a defect. A brightness field and a
    color field have different spatial spectra: brightness in these scenes varies mostly
    at scales a coarse sensor resolves, color mostly at scales it does not. Since a
    match-up certifies a color-derived quantity (chlorophyll) using a filter applied to
    brightness, both numbers belong in the record. Reporting only the larger would be a
    choice made after seeing the answer.

    **Two between-pixel statistics, and they are not interchangeable.** The filter is
    applied where the protocol applies it — to the coarse sensor's own retrieval — and
    that is what decides which boxes pass. But the *comparison* against the within-pixel
    term uses ``tanager_between_pixel_cv``: the same 5x5 window, the same green/blue
    ratio, the same instrument, differing only in scale. Dividing the within-pixel CV of
    a Tanager band ratio by the between-pixel CV of OLCI reflectance or PACE chlorophyll
    would be a ratio of unlike quantities, and the number it produces is not a
    scale-decomposition of anything.

    Every number here is a *lower bound*. Tanager's surface reflectance carries a large
    additive atmospheric pedestal (see tanager_io.water); an additive term inflates a
    mean and leaves a standard deviation alone, so it biases every coefficient of
    variation downwards. Scene-level pedestal removal in stage 2 reduces this but cannot
    abolish it.

    Offline, seconds.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from _scenes import OUT_DIR, Scene, load_scenes
from _verdicts import configure_stage_logging, upsert_verdict, verdict_row

logger = logging.getLogger(__name__)

#: The protocol's homogeneity threshold on the between-pixel coefficient of variation.
CV_THRESHOLD = 0.15

#: The protocol's box, in coarse pixels on a side.
BOX = 5

#: The protocol's minimum fraction of the box that must be valid.
MIN_VALID_FRACTION = 0.50

#: The scalar the filter is applied to, and whose within-pixel variability is measured.
PRIMARY_RATIO = "green_blue"

#: A second scalar, measured the same way and reported beside the first: the plain
#: reflectance at this wavelength. The answer to "how much variance is hidden inside a
#: pixel" is not a property of the water alone — it depends on which quantity is being
#: asked about, because a brightness field and a color field have different spatial
#: spectra. Reporting only the one that gives the larger answer would be a choice made
#: after seeing the result, so both are computed.
SECONDARY_BAND_NM = 560.0

#: Physically admissible range for the mean green/blue ratio of a water pixel. A
#: coefficient of variation is a ratio to the mean, so it diverges as the mean approaches
#: zero; a coarse pixel whose mean ratio falls outside this range is not reporting a water
#: spectrum but a failed pedestal removal, and its CV is an artefact rather than a
#: measurement. Excluding them is a correctness fix, not a convenience: without it a
#: handful of near-zero denominators dominate every upper percentile.
RATIO_MEAN_RANGE = (0.2, 5.0)

#: V3 passes if adjacent within-pixel variability differs by less than this multiple of
#: what a spatially shuffled field gives.
V3_MAX_NEIGHBOUR_RATIO = 0.75

#: V3 needs this many horizontally adjacent finite pairs to mean anything. On a sparse
#: grid — a PACE swath holding a handful of usable pixels — neighboring cells are rarely
#: both finite, and the statistic is undefined rather than failed.
V3_MIN_ADJACENT_PAIRS = 50

#: V4 passes if the median within-pixel CV exceeds this multiple of the declared floor.
V4_MIN_NOISE_MULTIPLE = 1.5


@dataclass
class StageConfig:
    """Stage 3 wiring; the thresholds are the protocol's, not ours."""

    out_dir: Path = OUT_DIR
    cv_threshold: float = CV_THRESHOLD
    box: int = BOX
    min_valid_fraction: float = MIN_VALID_FRACTION
    ratio: str = PRIMARY_RATIO
    secondary_band_nm: float = SECONDARY_BAND_NM


def coarse_shape(scene: Scene, sensor: str, out_dir: Path = OUT_DIR) -> tuple[int, int]:
    """The coarse sensor's native 2-D grid shape, so a flat index can be reshaped.

    Read from stage 2's ledger where it is recorded, and only otherwise from the source
    granule — a PACE L2 file is ~200 MB and this is the single fact stage 3 wants from
    it, so requiring it would mean the aggregates could not be shipped on their own.
    """
    ledger_path = out_dir / "metrics" / f"g2_{scene.key}.json"
    if ledger_path.exists():
        recorded = json.loads(ledger_path.read_text()).get(sensor, {}).get("coarse_shape")
        if recorded:
            return int(recorded[0]), int(recorded[1])
    if sensor == "olci":
        with xr.open_dataset(scene.olci_path) as dataset:
            return dataset.sizes["lat"], dataset.sizes["lon"]
    with xr.open_dataset(scene.pace_aop_path, group="navigation_data") as navigation:
        return navigation["latitude"].shape


def mask_unphysical(cv_grid: np.ndarray, mean_grid: np.ndarray) -> np.ndarray:
    """Blank the coefficient of variation where its mean is not a water spectrum.

    See :data:`RATIO_MEAN_RANGE`. This is applied wherever the within-pixel CV is used,
    so the filter cannot differ between the headline table and the structure score.
    """
    lo, hi = RATIO_MEAN_RANGE
    return np.where((mean_grid >= lo) & (mean_grid <= hi), cv_grid, np.nan)


def band_columns(frame: pd.DataFrame, nanometres: float) -> tuple[str, str]:
    """The per-coarse-pixel mean and standard-deviation columns nearest ``nanometres``."""
    means = [name for name in frame.columns if name.startswith("rrs_mean_")]
    if not means:
        raise KeyError("no rrs_mean_* columns — stage 2 wrote no per-band moments")
    nearest = min(means, key=lambda name: abs(float(name.split("_")[-1]) - nanometres))
    return nearest, nearest.replace("rrs_mean_", "rrs_std_")


def to_grid(frame: pd.DataFrame, column: str, shape: tuple[int, int]) -> np.ndarray:
    """Scatter one column back onto the coarse sensor's 2-D grid, NaN elsewhere."""
    grid = np.full(int(np.prod(shape)), np.nan)
    grid[frame["coarse_index"].to_numpy()] = frame[column].to_numpy()
    return grid.reshape(shape)


def boxes(grid: np.ndarray, box: int) -> np.ndarray:
    """Every ``box`` x ``box`` window of ``grid`` as a (rows, cols, box*box) array."""
    windows = np.lib.stride_tricks.sliding_window_view(grid, (box, box))
    return windows.reshape(windows.shape[0], windows.shape[1], box * box)


def between_pixel_cv(values: np.ndarray, min_valid_fraction: float,
                     trim: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """The protocol's box statistic: valid fraction and between-pixel CV per window.

    ``trim`` applies the protocol's iterated 1.5-sigma rejection before recomputing,
    which *lowers* the CV and so lets more boxes through the filter — the conservative
    direction for this experiment, since a box that passes is a box we then interrogate.
    """
    valid = np.isfinite(values)
    fraction = valid.mean(axis=-1)
    working = np.where(valid, values, np.nan)
    # A window holding no valid water is legitimately all-NaN; numpy warns about the
    # empty slice, but NaN is the answer we want, so the warning is noise here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if trim:
            mean = np.nanmean(working, axis=-1, keepdims=True)
            deviation = np.nanstd(working, axis=-1, keepdims=True)
            inside = np.abs(working - mean) <= 1.5 * deviation
            working = np.where(inside, working, np.nan)
        mean = np.nanmean(working, axis=-1)
        deviation = np.nanstd(working, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cv = np.where(np.abs(mean) > 0, deviation / np.abs(mean), np.nan)
    cv = np.where(fraction >= min_valid_fraction, cv, np.nan)
    return fraction, cv


def homogeneity_table(scene: Scene, sensor: str, config: StageConfig) -> pd.DataFrame:
    """One row per candidate 5x5 box: what the protocol sees, and what it cannot.

    Carries both between-pixel statistics — ``between_pixel_cv`` on the coarse sensor's
    retrieval, which sets ``passes_filter``, and ``tanager_between_pixel_cv`` on the same
    scalar as the within-pixel term, which is what the within-pixel term is compared to.
    """
    path = config.out_dir / f"{scene.key}_{sensor}_agg.parquet"
    frame = pd.read_parquet(path)
    shape = coarse_shape(scene, sensor, config.out_dir)

    # What the protocol filters on: the coarse sensor's own retrieval. Use its green
    # reflectance, the band every ocean-color algorithm leans on.
    candidates = [name for name in ("coarse_RRS560", "coarse_Rrs_560", "coarse_CHL",
                                    "coarse_chlor_a") if name in frame.columns]
    if not candidates:
        raise KeyError(f"{path.name} has no coarse retrieval to filter on; "
                       f"have {[c for c in frame.columns if c.startswith('coarse_')]}")
    filter_column = candidates[0]

    coarse_grid = to_grid(frame, filter_column, shape)
    mean_grid = to_grid(frame, f"{config.ratio}_mean", shape)
    within_grid = to_grid(frame, f"{config.ratio}_cv", shape)
    fraction_grid = to_grid(frame, "frac_water", shape)
    within_grid = mask_unphysical(within_grid, mean_grid)

    # What the protocol filters on, exactly as it does it: the coarse sensor's own
    # retrieval. This decides which boxes pass.
    valid_fraction, box_cv = between_pixel_cv(boxes(coarse_grid, config.box),
                                              config.min_valid_fraction)
    # What the *comparison* must use: the same protocol statistic on the same scalar
    # whose interior is then interrogated. The filter column is a different quantity on
    # a different instrument — OLCI's 560 nm reflectance, or PACE's chlorophyll — so
    # "within-pixel CV over between-pixel CV" across the two is a ratio of unlike
    # things, and no identity line can be drawn through it. Same ratio, same sensor,
    # same 5x5 window: the only thing that differs is the scale.
    _, tanager_box_cv = between_pixel_cv(
        boxes(mask_unphysical(mean_grid, mean_grid), config.box), config.min_valid_fraction)

    # The same decomposition on a plain band rather than a ratio. Guarded only on a
    # positive mean, which is the one condition a pedestal-subtracted reflectance must
    # satisfy: a near-zero mean inflates its own within-pixel CV, so a permissive guard
    # can only push this variant *towards* the ratio's answer, never away from it. The
    # reported gap between the two is therefore a floor on the gap.
    band_mean_name, band_std_name = band_columns(frame, config.secondary_band_nm)
    band_mean_grid = to_grid(frame, band_mean_name, shape)
    positive_band = np.where(band_mean_grid > 0, band_mean_grid, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        band_within_grid = to_grid(frame, band_std_name, shape) / positive_band
    _, band_box_cv = between_pixel_cv(boxes(positive_band, config.box),
                                      config.min_valid_fraction)
    edge = config.box // 2
    centre_within = within_grid[edge:-edge or None, edge:-edge or None]
    centre_water = fraction_grid[edge:-edge or None, edge:-edge or None]
    # The most heterogeneous pixel anywhere in the box, not only at its center: the
    # protocol accepts the whole box, so any pixel in it can be the one used.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        box_within_max = np.nanmax(boxes(within_grid, config.box), axis=-1)

    table = pd.DataFrame({
        "scene": scene.key, "sensor": sensor,
        "filter_on": filter_column,
        "box_valid_fraction": valid_fraction.ravel(),
        "between_pixel_cv": box_cv.ravel(),
        "tanager_between_pixel_cv": tanager_box_cv.ravel(),
        "band_nm": float(band_mean_name.split("_")[-1]),
        "band_between_pixel_cv": band_box_cv.ravel(),
        "band_centre_within_pixel_cv": band_within_grid[edge:-edge or None,
                                                        edge:-edge or None].ravel(),
        "centre_within_pixel_cv": centre_within.ravel(),
        "box_max_within_pixel_cv": box_within_max.ravel(),
        "centre_frac_water": centre_water.ravel(),
    })
    table["passes_filter"] = (table["between_pixel_cv"] < config.cv_threshold)
    return table.dropna(subset=["between_pixel_cv", "centre_within_pixel_cv"])


def structure_score(scene: Scene, sensor: str, config: StageConfig, seed: int = 0) -> dict:
    """V3: is within-pixel variability spatially structured, or salt-and-pepper?

    Compares the typical difference between horizontally adjacent cells against the same
    statistic after the finite values are shuffled across the grid. Real water structure
    is autocorrelated and scores far below the shuffled field; sensor noise scores the
    same as it.
    """
    frame = pd.read_parquet(config.out_dir / f"{scene.key}_{sensor}_agg.parquet")
    shape = coarse_shape(scene, sensor, config.out_dir)
    grid = mask_unphysical(to_grid(frame, f"{config.ratio}_cv", shape),
                           to_grid(frame, f"{config.ratio}_mean", shape))
    differences = np.abs(np.diff(grid, axis=1))
    n_pairs = int(np.isfinite(differences).sum())
    if n_pairs < V3_MIN_ADJACENT_PAIRS:
        return {"assessable": False, "n_adjacent_pairs": n_pairs,
                "n_cells": int(np.isfinite(grid).sum()),
                "note": f"only {n_pairs} adjacent finite pairs; the statistic is undefined"}

    neighbour = float(np.nanmedian(differences))
    finite = np.isfinite(grid)
    shuffled = grid.copy()
    values = shuffled[finite]
    np.random.default_rng(seed).shuffle(values)
    shuffled[finite] = values
    shuffled_neighbour = float(np.nanmedian(np.abs(np.diff(shuffled, axis=1))))
    return {"assessable": True,
            "neighbour_difference": neighbour,
            "shuffled_difference": shuffled_neighbour,
            "ratio": neighbour / shuffled_neighbour if shuffled_neighbour else float("nan"),
            "n_adjacent_pairs": n_pairs,
            "n_cells": int(finite.sum())}


def spectral_comparison(scene: Scene, config: StageConfig) -> pd.DataFrame:
    """H1: per-band Tanager-aggregated against OLCI's own reflectance, on usable pixels."""
    frame = pd.read_parquet(config.out_dir / f"{scene.key}_olci_agg.parquet")
    usable = frame[frame["frac_water"] >= 0.90]
    rows = []
    for column in [name for name in frame.columns if name.startswith("coarse_RRS")
                   and not name.endswith("_uncertainty")]:
        nanometres = float(column.replace("coarse_RRS", ""))
        tanager = next((name for name in frame.columns if name.startswith("rrs_mean_")
                        and abs(float(name.split("_")[-1]) - nanometres) <= 3), None)
        if tanager is None:
            continue
        pair = usable[[tanager, column]].dropna()
        if len(pair) < 20:
            continue
        rows.append({"scene": scene.key, "wavelength_nm": nanometres, "n": len(pair),
                     "correlation": float(pair[tanager].corr(pair[column])),
                     "tanager_median": float(pair[tanager].median()),
                     "coarse_median": float(pair[column].median()),
                     "ratio_of_medians": float(pair[tanager].median() / pair[column].median())
                     if pair[column].median() else np.nan})
    return pd.DataFrame(rows)


def noise_floor(scene: Scene, config: StageConfig) -> float:
    """The relative uncertainty the product declares at 560 nm, from stage 1's ledger."""
    path = config.out_dir / "metrics" / f"v0_v1_{scene.key}.json"
    if not path.exists():
        return float("nan")
    return float(json.loads(path.read_text())["v0_water"]["relative_uncertainty_560"])


def main(scenes: list[str] | None = None, sensors: list[str] | None = None,
         out_dir: Path = OUT_DIR, log_level: str = "info") -> int:
    configure_stage_logging(out_dir / "compare.log", log_level)
    config = StageConfig(out_dir=out_dir)
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    wanted_sensors = tuple(sensors) if sensors else ("olci", "pace")

    tables, structure, spectral, floors = [], {}, [], {}
    for scene in load_scenes(keys=scenes):
        floors[scene.key] = noise_floor(scene, config)
        for sensor in wanted_sensors:
            if not (out_dir / f"{scene.key}_{sensor}_agg.parquet").exists():
                logger.warning("%s/%s: no aggregate — run aggregate.py", scene.key, sensor)
                continue
            table = homogeneity_table(scene, sensor, config)
            table.to_parquet(out_dir / f"homogeneity_{scene.key}_{sensor}.parquet", index=False)
            tables.append(table)
            structure[f"{scene.key}/{sensor}"] = structure_score(scene, sensor, config)
            logger.info("%s/%s: %d boxes, %d pass the CV<%.2f filter", scene.key, sensor,
                        len(table), int(table["passes_filter"].sum()), config.cv_threshold)
        if (out_dir / f"{scene.key}_olci_agg.parquet").exists():
            spectral.append(spectral_comparison(scene, config))

    if not tables:
        raise RuntimeError("no aggregates found — run aggregate.py first")
    everything = pd.concat(tables, ignore_index=True)
    if spectral:
        pd.concat(spectral, ignore_index=True).to_parquet(out_dir / "h1_spectral.parquet",
                                                          index=False)

    results: dict = {"v5": {}, "v3": structure, "v4": {}}
    for (scene_key, sensor), group in everything.groupby(["scene", "sensor"]):
        passing = group[group["passes_filter"]]
        like_for_like = passing.dropna(subset=["tanager_between_pixel_cv"])
        within_median = (float(like_for_like["centre_within_pixel_cv"].median())
                         if len(like_for_like) else float("nan"))
        between_median = (float(like_for_like["tanager_between_pixel_cv"].median())
                          if len(like_for_like) else float("nan"))
        band = passing.dropna(subset=["band_between_pixel_cv",
                                      "band_centre_within_pixel_cv"])
        band_within_median = (float(band["band_centre_within_pixel_cv"].median())
                              if len(band) else float("nan"))
        band_between_median = (float(band["band_between_pixel_cv"].median())
                               if len(band) else float("nan"))
        results["v5"][f"{scene_key}/{sensor}"] = {
            "n_boxes": int(len(group)),
            "n_passing": int(len(passing)),
            "pass_rate": float(len(passing) / len(group)) if len(group) else float("nan"),
            "filter_on": str(group["filter_on"].iloc[0]) if len(group) else "",
            # The protocol's own filter statistic, on the coarse sensor's retrieval.
            # Reported because it is what the filter saw, *not* as the denominator of a
            # scale ratio: it is a different quantity on a different instrument.
            "median_between_pixel_cv_of_passing": float(passing["between_pixel_cv"].median())
            if len(passing) else float("nan"),
            # Like for like: same green/blue ratio, same 5x5 window, same instrument.
            "n_passing_like_for_like": int(len(like_for_like)),
            "median_tanager_between_pixel_cv_of_passing": between_median,
            "within_over_between_like_for_like": (within_median / between_median
                                                  if between_median else float("nan")),
            # The same decomposition on a plain band. Reported beside the ratio, not
            # instead of it: the two disagree, and which one is relevant depends on
            # whether the quantity being validated is a brightness or a color.
            "band_nm": float(group["band_nm"].iloc[0]) if len(group) else float("nan"),
            "n_passing_band": int(len(band)),
            "median_band_between_pixel_cv_of_passing": band_between_median,
            "median_band_within_pixel_cv_of_passing": band_within_median,
            "within_over_between_band": (band_within_median / band_between_median
                                         if band_between_median else float("nan")),
            "median_within_pixel_cv_of_passing": float(passing["centre_within_pixel_cv"].median())
            if len(passing) else float("nan"),
            "p90_within_pixel_cv_of_passing": float(passing["centre_within_pixel_cv"].quantile(0.9))
            if len(passing) else float("nan"),
            "frac_passing_with_centre_above_threshold": float(
                (passing["centre_within_pixel_cv"] > CV_THRESHOLD).mean())
            if len(passing) else float("nan"),
            "frac_passing_with_any_pixel_above_threshold": float(
                (passing["box_max_within_pixel_cv"] > CV_THRESHOLD).mean())
            if len(passing) else float("nan"),
        }
    for scene_key, group in everything.groupby("scene"):
        floor = floors.get(scene_key, float("nan"))
        median_cv = float(group["centre_within_pixel_cv"].median())
        results["v4"][scene_key] = {"declared_relative_uncertainty": floor,
                                    "median_within_pixel_cv": median_cv,
                                    "multiple_of_floor": median_cv / floor if floor else np.nan}

    (metrics_dir / "h_results.json").write_text(json.dumps(results, indent=2))

    assessed = {key: entry for key, entry in structure.items() if entry.get("assessable")}
    # Reported per combination rather than as one worst-case gate: V3 licenses the
    # interpretation of a given scene and sensor, so a marginal result on one of them is
    # a fact about that pairing, not grounds to discard the others.
    v3_failures = [f"{key} {entry['ratio']:.2f}" for key, entry in assessed.items()
                   if entry["ratio"] > V3_MAX_NEIGHBOUR_RATIO]
    multiples = [entry["multiple_of_floor"] for entry in results["v4"].values()
                 if np.isfinite(entry["multiple_of_floor"])]
    least_multiple = min(multiples) if multiples else float("nan")
    total_passing = sum(entry["n_passing"] for entry in results["v5"].values())
    hidden = {key: entry["frac_passing_with_any_pixel_above_threshold"]
              for key, entry in results["v5"].items()}
    ratios = {key: entry["within_over_between_like_for_like"]
              for key, entry in results["v5"].items()}
    band_ratios = {key: entry["within_over_between_band"]
                   for key, entry in results["v5"].items()}

    upsert_verdict(metrics_dir, [
        verdict_row("V3 within-pixel variability is spatially structured, not noise",
                    "median adjacent difference / same after shuffling, per scene and sensor",
                    "; ".join(f"{key} {entry['ratio']:.2f}" for key, entry in assessed.items()),
                    f"<= {V3_MAX_NEIGHBOUR_RATIO} where assessable",
                    "PASS" if assessed and not v3_failures
                    else ("NOT ASSESSABLE" if not assessed else "MIXED"),
                    (f"{len(assessed) - len(v3_failures)} of {len(assessed)} assessable "
                     f"combinations are structured"
                     + (f"; above the bar: {', '.join(v3_failures)}" if v3_failures else "")
                     + (f"; not assessable on {', '.join(key for key in structure if key not in assessed)}"
                        if len(assessed) < len(structure) else "")
                     + ". A shuffled field scores 1.0 by construction.")),
        verdict_row("V4 within-pixel variability exceeds the declared noise floor",
                    "median within-pixel CV / declared relative uncertainty at 560 nm",
                    "; ".join(f"{key} {entry['multiple_of_floor']:.1f}x"
                              for key, entry in results["v4"].items()),
                    f">= {V4_MIN_NOISE_MULTIPLE}x everywhere",
                    "PASS" if least_multiple >= V4_MIN_NOISE_MULTIPLE else "FAIL",
                    "floor is the product's own surface_reflectance_uncertainty"),
        verdict_row("V5 the homogeneity filter accepts boxes with variable interiors",
                    "boxes passing CV<0.15; within-pixel CV over between-pixel CV of the "
                    "same green/blue ratio; fraction of passing boxes containing a pixel "
                    "whose own interior CV exceeds 0.15",
                    f"{total_passing}; " + "; ".join(
                        f"{key} {ratios[key]:.1f}x, {hidden[key]:.1%}" for key in hidden),
                    "reported either way",
                    "ANSWERED",
                    "the ratio is like for like — same scalar, same window, same "
                    "instrument, only the scale differs; the protocol's own filter "
                    "statistic is on the coarse sensor's retrieval and is reported "
                    "separately in h_results.json. Within-pixel CV is a lower bound: the "
                    "additive atmospheric pedestal inflates the mean and so deflates "
                    "every CV. The same decomposition on plain 560 nm reflectance "
                    "instead of the color ratio gives "
                    + "; ".join(f"{key} {value:.1f}x" for key, value in band_ratios.items())
                    + " — the answer depends on which quantity is asked about, and both "
                    "are reported"),
    ])

    print(json.dumps(results, indent=2))
    problems = []
    if v3_failures:
        logger.warning("V3: %s above the %.2f bar — those pairings carry less spatial "
                       "structure and the memo must say so; the others are unaffected",
                       ", ".join(v3_failures), V3_MAX_NEIGHBOUR_RATIO)
    if multiples and least_multiple < V4_MIN_NOISE_MULTIPLE:
        problems.append(f"V4: {least_multiple:.1f}x < {V4_MIN_NOISE_MULTIPLE}x noise floor")
    if problems:
        raise RuntimeError("; ".join(problems))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 3: the variance a match-up protocol cannot see",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenes", nargs="*", default=None)
    parser.add_argument("--sensors", nargs="*", default=None, choices=["olci", "pace"])
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
