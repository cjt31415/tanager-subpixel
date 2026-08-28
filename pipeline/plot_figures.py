#!/usr/bin/env python3
"""
    plot_figures.py: stage 5 — the figures.

    Input:  outputs/homogeneity_*.parquet, tasking_boxes.parquet (stages 3-4);
            outputs/metrics/v0_v1_<scene>.json (stage 1);
            data/tanager/<scene>/*.h5 for the resolution ladder.
    Output: outputs/figures/f1_homogeneity.png  the headline
            outputs/figures/f2_ladder.png       the same water at three resolutions
            outputs/figures/f3_pedestal.png     what the product does over water
            outputs/figures/f4_tasking.png      where a scene would buy the most

    Colors are the validated categorical palette (blue/orange/aqua, checked for
    color-vision separation and lightness band) with scene as hue and sensor as marker
    shape, so identity never rests on color alone.

    Offline. ~1 min, dominated by reading a few Tanager bands for the ladder.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
import tanager_io as tio
from _scenes import OUT_DIR, SCENE_INVENTORY, Scene, load_scenes
from _verdicts import configure_stage_logging
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)

#: Validated categorical palette, one hue per scene (dataviz slots 1-3).
SCENE_COLOUR = {"sf_bay": "#2a78d6", "lake_ontario": "#eb6834", "loreto": "#1baf7a"}
SCENE_LABEL = {"sf_bay": "San Francisco Bay", "lake_ontario": "Lake Ontario",
               "loreto": "Loreto, Gulf of California"}
#: Sensor as a second, non-color encoding.
SENSOR_MARKER = {"olci": "o", "pace": "^"}
SENSOR_LABEL = {"olci": "OLCI ~550 m", "pace": "PACE ~1.2 km"}

INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8880"
SURFACE = "#fcfcfb"
SEQUENTIAL = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab"]

#: The protocol's homogeneity threshold, drawn on the headline figure.
CV_THRESHOLD = 0.15


@dataclass
class StageConfig:
    """Stage 5 wiring."""

    out_dir: Path = OUT_DIR
    dpi: int = 200


def style() -> None:
    """Recessive axes, thin marks, ink-colored text — applied once."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK_MUTED, "axes.linewidth": 0.8,
        "axes.labelcolor": INK_SECONDARY, "axes.titlecolor": INK,
        "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_SECONDARY, "ytick.labelcolor": INK_SECONDARY,
        "text.color": INK, "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
        "legend.frameon": False,
    })


def load_homogeneity(out_dir: Path) -> pd.DataFrame:
    """Every homogeneity table stage 3 wrote."""
    frames = [pd.read_parquet(path) for path in sorted(out_dir.glob("homogeneity_*.parquet"))]
    if not frames:
        raise FileNotFoundError("no homogeneity tables — run compare.py first")
    return pd.concat(frames, ignore_index=True)


def figure_homogeneity(table: pd.DataFrame, config: StageConfig) -> Path:
    """F1, the headline: what the filter measures against what it cannot.

    A scatter, because the claim is about the *relationship* between two quantities —
    and specifically about which side of the identity line the cloud sits on.

    Both axes are the same green/blue ratio measured by the same instrument over the
    same 5x5 window, so the identity line is a real statement about scale. The protocol's
    own filter statistic — computed on the coarse sensor's retrieval, a different
    quantity on a different instrument — decides which boxes are drawn, and nothing more:
    plotting it on the x-axis would make the identity line meaningless.
    """
    style()
    figure, axes = plt.subplots(figsize=(7.2, 5.4))
    passing = table[table["passes_filter"]].dropna(subset=["tanager_between_pixel_cv"])

    limit = (0.004, 3.0)
    axes.plot(limit, limit, color=INK_MUTED, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    axes.annotate("equal variability\ninside and between pixels", xy=(0.9, 0.9),
                  xytext=(1.15, 0.55), fontsize=7.5, color=INK_MUTED, ha="left",
                  va="center", linespacing=1.4)
    axes.axhline(CV_THRESHOLD, color=INK_MUTED, linewidth=0.8, alpha=0.7, zorder=1)

    for (scene, sensor), group in passing.groupby(["scene", "sensor"]):
        if group.empty:
            continue
        axes.scatter(group["tanager_between_pixel_cv"], group["centre_within_pixel_cv"],
                     s=16, marker=SENSOR_MARKER[sensor], facecolor=SCENE_COLOUR[scene],
                     edgecolor=SURFACE, linewidth=0.4, alpha=0.55, zorder=3,
                     label=f"{SCENE_LABEL[scene]} · {SENSOR_LABEL[sensor]}")
        axes.scatter([group["tanager_between_pixel_cv"].median()],
                     [group["centre_within_pixel_cv"].median()],
                     s=150, marker=SENSOR_MARKER[sensor], facecolor=SCENE_COLOUR[scene],
                     edgecolor=INK, linewidth=1.4, zorder=5)

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(*limit)
    axes.set_ylim(*limit)
    axes.set_xlabel("between-pixel CV across the 5×5 box  —  the scale the filter screens on")
    axes.set_ylabel("within-pixel CV inside one pixel  —  the scale it cannot see")
    axes.set_title("The homogeneity a match-up protocol assumes, and the water it accepts",
                   fontsize=11.5, pad=32, loc="left")
    axes.text(0, 1.015, "Only boxes passing the standard CV < 0.15 filter, applied where "
              "the protocol applies it — to the coarse sensor's own retrieval.\nBoth axes "
              "are Tanager's green/blue ratio, so only the scale differs. Large outlined "
              "markers are per-group medians.",
              transform=axes.transAxes, fontsize=8, color=INK_SECONDARY, va="bottom",
              linespacing=1.6)
    axes.text(limit[1] * 0.42, CV_THRESHOLD * 1.14, "filter threshold 0.15", fontsize=7.5,
              color=INK_SECONDARY, va="bottom", ha="right")
    axes.grid(True, which="major", alpha=0.7, zorder=0)
    axes.legend(loc="upper left", fontsize=8, labelcolor=INK_SECONDARY,
                handletextpad=0.4, borderaxespad=0.8)

    above = (passing["centre_within_pixel_cv"] > CV_THRESHOLD).mean()
    ratio = (passing["centre_within_pixel_cv"].median()
             / passing["tanager_between_pixel_cv"].median())
    axes.text(0.985, 0.035, f"{above:.0%} of accepted boxes sit above the threshold they "
              f"just cleared\nmedian within-pixel CV is {ratio:.1f}× the between-pixel CV "
              f"of the same field\nn = {len(passing):,} boxes  ·  Tanager's CV is a lower bound",
              transform=axes.transAxes, fontsize=8.5, color=INK, ha="right", va="bottom",
              linespacing=1.6,
              bbox={"facecolor": SURFACE, "edgecolor": "#e6e5e0", "boxstyle": "round,pad=0.5"})

    path = config.out_dir / "figures" / "f1_homogeneity.png"
    figure.tight_layout()
    figure.savefig(path, dpi=config.dpi, bbox_inches="tight")
    plt.close(figure)
    return path


#: Block-averaging factors from Tanager's 30 m onto each coarse sensor's pixel.
LADDER_FACTORS = ((1, "Tanager  30 m"), (18, "OLCI  ~550 m"), (40, "PACE  ~1.2 km"))

#: Side of the square window the ladder and the variance split are computed over, in
#: Tanager pixels. 480 x 30 m = 14.4 km, comfortably wider than a 5x5 PACE box.
LADDER_PATCH_PIXELS = 480


def variance_between(array: np.ndarray, factor: int) -> float:
    """The part of the total variance that survives averaging into ``factor`` blocks.

    The law of total variance splits the 30 m variance exactly into a *between-cell*
    part, which a coarse sensor can see, and a *within-cell* part, which it cannot. The
    two sum to the total by construction, which is why this is the honest ladder
    statistic — a ratio of panel variances is not, and can exceed one.

    It is also the one headline number in this experiment that needs no threshold and no
    second instrument: one field, one window, two scales.
    """
    if factor == 1:
        return float(np.nanvar(array))
    rows_, cols_ = array.shape
    trimmed = array[:rows_ // factor * factor, :cols_ // factor * factor]
    shaped = trimmed.reshape(trimmed.shape[0] // factor, factor,
                             trimmed.shape[1] // factor, factor)
    counts = np.isfinite(shaped).sum(axis=(1, 3))
    with np.errstate(invalid="ignore"):
        means = np.nanmean(shaped, axis=(1, 3))
    valid = counts > 0
    grand = float(np.nansum(means[valid] * counts[valid]) / counts[valid].sum())
    return float(np.nansum(counts[valid] * (means[valid] - grand) ** 2)
                 / counts[valid].sum())


def block(array: np.ndarray, factor: int) -> np.ndarray:
    """``array`` averaged into non-overlapping ``factor`` x ``factor`` blocks."""
    rows_, cols_ = array.shape
    trimmed = array[:rows_ // factor * factor, :cols_ // factor * factor]
    shaped = trimmed.reshape(trimmed.shape[0] // factor, factor,
                             trimmed.shape[1] // factor, factor)
    with np.errstate(invalid="ignore"):
        return np.nanmean(shaped, axis=(1, 3))


def ladder_patch(scene_key: str, scalar: str = "ratio") -> tuple[np.ndarray, Scene]:
    """One scene's central water window, and the scene.

    ``scalar="ratio"`` gives the green/blue ratio — the same construction as stage 3's
    30 m scalar, pedestal-subtracted 560/490 and guarded to the physically admissible
    range — so the variance split and the CV table describe the same quantity.

    ``scalar="band"`` gives plain 560 nm reflectance, so the split can be reported for a
    brightness field as well as a color one. It needs no range guard: the split is a
    ratio of variances, and variance is shift-invariant, so the additive atmospheric
    pedestal cancels out of it exactly.
    """
    scene = next(s for s in load_scenes(keys=[scene_key]))
    dataset = tio.open_sr(scene.sr_path)
    water_mask = tio.water_mask(dataset)
    water = water_mask.values
    wavelength = dataset["wavelength"].values
    # The same two bands stage 2 uses, through the same library calls, so the split and
    # the CV table cannot drift apart in their pedestal or their units.
    pair = dataset[["surface_reflectance"]].isel(wavelength=[
        tio.band_index(wavelength, 560.0), tio.band_index(wavelength, 490.0)]).load()
    rrs = tio.to_rrs(pair, water=water_mask)
    if scalar == "band":
        field = np.where(water, rrs.isel(wavelength=0).values, np.nan)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            field = np.where(water, tio.band_ratio(pair, 560.0, 490.0, source=rrs).values,
                             np.nan)
        # The same physical guard stage 3 applies: a ratio outside this range is a failed
        # pedestal removal, not a water spectrum, and it would dominate both the color
        # stretch and any dispersion statistic computed from the panel.
        field = np.where((field >= 0.2) & (field <= 5.0), field, np.nan)

    # A window over the front: the densest water quarter of the scene.
    rows = np.flatnonzero(water.any(axis=1))
    cols = np.flatnonzero(water.any(axis=0))
    size = LADDER_PATCH_PIXELS
    row0 = max(int(np.median(rows)) - size // 2, 0)
    col0 = max(int(np.median(cols)) - size // 2, 0)
    patch = field[row0:row0 + size, col0:col0 + size]
    dataset.close()
    return patch, scene


def variance_splits(config: StageConfig) -> dict:
    """The exact between/within-cell variance split for every scene, both pixel sizes.

    The CV table answers the protocol's question and inherits the protocol's threshold.
    This answers the same question with no threshold at all — and for all three scenes,
    not only the one the ladder figure happens to draw.

    Computed for both scalars, because the CV table's two variants disagree and this is
    the statistic that says whether that disagreement is about *scale* or only about
    magnitude.
    """
    splits: dict = {}
    for scene in load_scenes():
        if not scene.sr_path.exists():
            logger.warning("%s: no cube on disk — no variance split", scene.key)
            continue
        entry: dict = {"window_km": LADDER_PATCH_PIXELS * 30 / 1000}
        for scalar in ("ratio", "band"):
            patch, _ = ladder_patch(scene.key, scalar)
            total = float(np.nanvar(patch))
            finite = int(np.isfinite(patch).sum())
            if not np.isfinite(total) or total <= 0 or finite < LADDER_PATCH_PIXELS:
                logger.warning("%s/%s: %d finite px, total variance %.3g — split undefined",
                               scene.key, scalar, finite, total)
                continue
            entry[scalar] = {
                "n_finite_30m_px": finite,
                "total_variance_30m": total,
                **{label.split()[0].lower(): {
                    "factor": factor,
                    "between_cell_fraction": variance_between(patch, factor) / total,
                    "within_cell_fraction": 1 - variance_between(patch, factor) / total}
                   for factor, label in LADDER_FACTORS if factor > 1},
            }
            logger.info("%s/%s: %.0f%% of the 30 m variance is hidden inside an OLCI "
                        "cell, %.0f%% inside a PACE pixel", scene.key, scalar,
                        100 * entry[scalar]["olci"]["within_cell_fraction"],
                        100 * entry[scalar]["pace"]["within_cell_fraction"])
        if "ratio" in entry:
            splits[scene.key] = entry
    if not splits:
        raise FileNotFoundError("no Tanager cubes on disk — cannot compute the split")
    path = config.out_dir / "metrics" / "variance_split.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(splits, indent=2))
    return splits


def figure_ladder(config: StageConfig, scene_key: str = "lake_ontario") -> Path:
    """F2: the same stretch of water at 30 m, ~550 m and ~1.2 km.

    Block-averaging one Tanager scene is the honest way to show the ladder: every panel
    is the same measurement, so the only thing changing between them is the pixel.
    """
    style()
    patch, scene = ladder_patch(scene_key)
    size = LADDER_PATCH_PIXELS

    cmap = mpl.colors.LinearSegmentedColormap.from_list("seq", SEQUENTIAL)
    finite = patch[np.isfinite(patch)]
    low, high = np.percentile(finite, [8, 92])
    total_variance = float(np.nanvar(patch))

    figure, axes_row = plt.subplots(1, 3, figsize=(10.4, 3.9))
    for axes, (factor, label) in zip(axes_row, LADDER_FACTORS, strict=True):
        image = block(patch, factor) if factor > 1 else patch
        axes.imshow(image, cmap=cmap, vmin=low, vmax=high, interpolation="nearest")
        axes.set_title(label, fontsize=10, color=INK, pad=8)
        axes.set_xticks([])
        axes.set_yticks([])
        for spine in axes.spines.values():
            spine.set_visible(True)
            spine.set_color(INK_MUTED)
        between = variance_between(patch, factor) / total_variance if total_variance else np.nan
        axes.text(0.5, -0.13, f"{image.shape[0]}×{image.shape[1]} cells\n"
                  f"{between:.0%} of the variance is between cells\n"
                  f"{1 - between:.0%} is hidden inside them",
                  transform=axes.transAxes, ha="center", fontsize=8.5,
                  color=INK if factor > 1 else INK_SECONDARY, linespacing=1.6)

    figure.suptitle(f"The same {size * 30 / 1000:.0f} km of water, three pixel sizes  ·  "
                    f"{SCENE_LABEL[scene_key]}, {scene.date}",
                    fontsize=11.5, color=INK, x=0.02, ha="left", y=1.02)
    figure.text(0.02, 0.945, "Green/blue reflectance ratio. Every panel is the same "
                "Tanager measurement, block-averaged — only the pixel changes. The split "
                "is exact (law of total variance);\nthe hidden fraction is an upper bound "
                "on hidden *water* variance, since the 30 m panel's pushbroom striping "
                "counts towards it.",
                fontsize=8, color=INK_SECONDARY, ha="left")
    path = config.out_dir / "figures" / "f2_ladder.png"
    figure.tight_layout()
    figure.savefig(path, dpi=config.dpi, bbox_inches="tight")
    plt.close(figure)
    return path


def figure_pedestal(config: StageConfig) -> Path:
    """F3: the measured water spectrum, and why absolute reflectance is not usable."""
    style()
    figure, axes = plt.subplots(figsize=(7.2, 4.6))
    for scene in load_scenes():
        path = config.out_dir / "metrics" / f"v0_v1_{scene.key}.json"
        if not path.exists():
            continue
        spectrum = json.loads(path.read_text())["v0_water"]["spectrum"]
        nanometres = np.array([float(key) for key in spectrum])
        median = np.array([entry["median"] for entry in spectrum.values()])
        order = np.argsort(nanometres)
        axes.plot(nanometres[order], median[order], color=SCENE_COLOUR[scene.key],
                  linewidth=2.0, marker="o", markersize=4.5,
                  markeredgecolor=SURFACE, markeredgewidth=0.6,
                  label=SCENE_LABEL[scene.key])

    axes.axhspan(0, 0.01, color="#e6e5e0", zorder=0)
    axes.text(2280, 0.013, "where water actually sits\nbeyond ~1000 nm", fontsize=7.5,
              color=INK_SECONDARY, ha="right", va="bottom", linespacing=1.4)
    axes.set_xlabel("wavelength (nm)")
    axes.set_ylabel("Tanager surface reflectance over water (median)")
    axes.set_title("A land-tuned atmospheric correction, seen from the water",
                   fontsize=11.5, pad=14, loc="left")
    axes.text(0, 1.015, "Liquid water is effectively black past 1000 nm. What remains is "
              "an additive path-reflectance pedestal of ~0.05–0.10.",
              transform=axes.transAxes, fontsize=8, color=INK_SECONDARY, va="bottom")
    axes.set_ylim(0, None)
    axes.grid(True, alpha=0.7, zorder=0)
    axes.legend(loc="upper right", fontsize=8.5, labelcolor=INK_SECONDARY)
    path = config.out_dir / "figures" / "f3_pedestal.png"
    figure.tight_layout()
    figure.savefig(path, dpi=config.dpi, bbox_inches="tight")
    plt.close(figure)
    return path


def figure_tasking(config: StageConfig) -> Path:
    """F4: the top float boxes, the open scenes, and the distance between them."""
    style()
    boxes = pd.read_parquet(config.out_dir / "tasking_boxes.parquet")
    ledger = json.loads((config.out_dir / "metrics" / "g4_tasking.json").read_text())
    scenes = pd.DataFrame(json.loads(
        SCENE_INVENTORY.read_text()))
    bbox = np.vstack(scenes["bbox"].to_numpy())

    figure, axes = plt.subplots(figsize=(10.4, 5.2))
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        plt.close(figure)
        figure = plt.figure(figsize=(10.4, 6.0))
        axes = figure.add_axes([0.02, 0.14, 0.96, 0.72], projection=ccrs.Robinson())
        axes.add_feature(cfeature.LAND, facecolor="#ecebe6", edgecolor="none", zorder=0)
        axes.add_feature(cfeature.COASTLINE, edgecolor=INK_MUTED, linewidth=0.4, zorder=1)
        axes.set_global()
        transform = {"transform": ccrs.PlateCarree()}
    except ImportError:  # pragma: no cover - cartopy is in the env, this is a courtesy
        axes.set_xlim(-180, 180)
        axes.set_ylim(-90, 90)
        transform = {}

    axes.scatter((bbox[:, 0] + bbox[:, 2]) / 2, (bbox[:, 1] + bbox[:, 3]) / 2,
                 s=14, marker="s", facecolor="none", edgecolor=INK_MUTED, linewidth=0.7,
                 zorder=3, label=f"open Tanager scenes (n = {len(scenes)})", **transform)

    ocean = boxes[(boxes["metric"] == "floats") & (boxes["window"] == "all")]
    axes.scatter((ocean["lon_west"] + ocean["lon_east"]) / 2,
                 (ocean["lat_min"] + ocean["lat_max"]) / 2,
                 s=ocean["n_floats"] * 6, marker="o", facecolor=SCENE_COLOUR["sf_bay"],
                 edgecolor=SURFACE, linewidth=0.8, alpha=0.85, zorder=5,
                 label="densest BGC-Argo boxes, one Tanager footprint", **transform)

    erie = boxes[boxes["metric"] == "in_situ_stations"]
    if not erie.empty:
        row = erie.iloc[0]
        axes.scatter([(row["lon_west"] + row["lon_east"]) / 2],
                     [(row["lat_min"] + row["lat_max"]) / 2], s=110, marker="*",
                     facecolor=SCENE_COLOUR["lake_ontario"], edgecolor=INK, linewidth=0.8,
                     zorder=6, label="proposed freshwater scene (western Lake Erie)",
                     **transform)

    top = ledger["ocean_top_by_floats"]
    axes.annotate(f"Ligurian Sea — {top['n_floats']} floats, "
                  f"{top['chla_profiles']} chlorophyll profiles\n"
                  f"nearest open scene {top['nearest_open_scene_km']:.0f} km away",
                  xy=((top["lon_west"] + top["lon_east"]) / 2,
                      (top["lat_min"] + top["lat_max"]) / 2),
                  xytext=(26, -52), textcoords="offset points", fontsize=8, color=INK,
                  ha="left", va="top", linespacing=1.5,
                  arrowprops={"arrowstyle": "-", "color": INK_MUTED, "linewidth": 0.8},
                  **transform)

    figure.text(0.02, 0.955, "Where a Tanager scene would buy the most", fontsize=12.5,
                color=INK, ha="left", va="top")
    figure.text(0.02, 0.905, "Marker area is distinct BGC-Argo floats inside one Tanager "
                f"footprint (0.30° × 0.20°). No open scene lies within "
                f"{ledger['min_nearest_scene_km_over_all_ocean_boxes']:.0f} km of any of them.",
                fontsize=8.5, color=INK_SECONDARY, ha="left", va="top")
    axes.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=8.5,
                labelcolor=INK_SECONDARY, ncols=3, columnspacing=2.2)
    path = config.out_dir / "figures" / "f4_tasking.png"
    figure.savefig(path, dpi=config.dpi, bbox_inches="tight")
    plt.close(figure)
    return path


def main(out_dir: Path = OUT_DIR, dpi: int = 200, log_level: str = "info") -> int:
    configure_stage_logging(out_dir / "plot_figures.log", log_level)
    config = StageConfig(out_dir=out_dir, dpi=dpi)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    table = load_homogeneity(out_dir)
    try:
        splits = variance_splits(config)
        logger.info("variance split written for %d scenes", len(splits))
    except FileNotFoundError as error:
        logger.warning("variance split skipped: %s", error)
    written = [figure_homogeneity(table, config), figure_pedestal(config)]
    if (out_dir / "tasking_boxes.parquet").exists():
        written.append(figure_tasking(config))
    else:
        logger.warning("no tasking_boxes.parquet — run tasking.py; skipping F4")
    try:
        written.append(figure_ladder(config))
    except (FileNotFoundError, StopIteration) as error:
        logger.warning("ladder skipped: %s", error)

    for path in written:
        logger.info("wrote %s (%.0f kB)", path.name, path.stat().st_size / 1e3)
    print(json.dumps([str(path.relative_to(out_dir)) for path in written], indent=2))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 5: figures",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
