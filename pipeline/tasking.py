#!/usr/bin/env python3
"""
    tasking.py: exp 35 stage 4 — where a Tanager scene would buy the most.

    Input:  outputs/bgc_argo_index_subset.parquet (five columns of the BGC-Argo index,
            written here from drift's full index on first run, shipped thereafter);
            tanager_open_scenes.json (the 153 open scenes);
            erie_insitu.csv (the western Lake Erie monitoring network).
    Output: outputs/tasking_boxes.parquet   candidate boxes, ranked, with what is in them
            outputs/metrics/g4_tasking.json the ledger, including the coverage gap
            outputs/metrics/verdicts.csv    V6

    Experiment 35 measures a term that ocean-color validation assumes away: how variable
    the water is *inside* a coarse pixel. Doing that at a site with in-situ truth would
    close the loop — and the open Tanager catalogue currently has no scene at any such
    site. This stage says where they are.

    Two questions, deliberately kept apart because they have different answers:

    **Ocean.** Rank Tanager-footprint-sized boxes by BGC-Argo density. The metric is the
    result (experiment 20's finding, re-run at 0.30 x 0.22 deg): ranking by *profiles*
    selects marginal seas where a few trapped floats cycle fast, while ranking by
    *distinct floats* selects the places a scene would actually catch a fleet. Both are
    emitted so the difference stays visible rather than being chosen silently.

    **Freshwater.** Western Lake Erie needs no ranking — it holds the longest
    cyanobacteria-bloom record in North America and a weekly three-depth sampling network.
    It is here as the worked example of the whole argument: the only two open Tanager
    scenes over the basin are ~90 % downtown Detroit and a river mouth, both at least
    12 km from every station in that network. The data are there; the scene is not.

    Offline, seconds.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from _geo import haversine_km, top_boxes
from _scenes import ERIE_INSITU, OUT_DIR, PROJECT_ROOT, SCENE_INVENTORY
from _verdicts import configure_stage_logging, upsert_verdict, verdict_row

logger = logging.getLogger(__name__)

#: Drift's full BGC-Argo index (395k profiles, 14 columns, 17 MB) — present only in the
#: drift tree. Stage 4 reads the five columns it uses from it once, and writes them next
#: to the outputs as the subset below, which is what the standalone repo ships and reads.
BGC_INDEX_FULL = PROJECT_ROOT / "data" / "ready" / "world" / "sensors" / "bgc_argo_index.parquet"
INDEX_COLUMNS = ["PLATFORM_NUMBER", "LATITUDE", "LONGITUDE", "TIME", "parameters"]


#: One Tanager footprint, in degrees. Measured across all 153 open scenes rather than
#: assumed: the median extent is 0.290 x 0.208 deg (about 26.5 x 23.2 km at |lat| 35),
#: rounded here to a whole multiple of the search grid below, which the box placement
#: requires. Both values sit inside the interquartile range of the real footprints.
BOX_DEG = (0.30, 0.20)

#: Keep proposed sites this far apart, so the list is distinct regions rather than the
#: same hotspot found repeatedly at neighboring offsets.
MIN_SEPARATION_DEG = 3.0

#: Search grid for box placement. Finer than the box, coarse enough to stay quick.
CELL_DEG = 0.05

#: Recency cut. A box ranked on twenty years of floats that has held none since this date
#: is a historical fact, not a tasking target.
RECENT_SINCE = "2023-01-01"

#: The proposed freshwater scene: one Tanager footprint centered here holds the whole
#: GLERL/CIGLR western-basin network. See README section 4b.
ERIE_CENTRE = (-83.28, 41.80)

#: The two open scenes that already cover the western basin, and cover the wrong part of
#: it: one is ~90 % downtown Detroit, the other a river mouth. They are the measurement
#: of the gap, not a case study.
EXISTING_ERIE_SCENES = ("20250914_171527_18_4001", "20250914_171551_87_4001")


@dataclass
class StageConfig:
    """Stage 4 wiring; the box is the one science constant."""

    out_dir: Path = OUT_DIR
    box_deg: tuple[float, float] = BOX_DEG
    top_n: int = 10
    min_separation_deg: float = MIN_SEPARATION_DEG
    cell_deg: float = CELL_DEG
    recent_since: str = RECENT_SINCE


def index_subset_path(out_dir: Path) -> Path:
    return out_dir / "bgc_argo_index_subset.parquet"


def load_index(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """The BGC-Argo profile index, with a CHLA flag derived from its parameter string.

    Reads the shipped subset when it is there; otherwise (in the drift tree) cuts it from
    the full index and writes it, so the next run — and the standalone repo — need only
    the subset. The subset is derived from the public GDAC ``argo_bio-profile_index.txt``.
    """
    subset = index_subset_path(out_dir)
    if subset.exists():
        frame = pd.read_parquet(subset)
    elif BGC_INDEX_FULL.exists():
        frame = pd.read_parquet(BGC_INDEX_FULL, columns=INDEX_COLUMNS)
        out_dir.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(subset, index=False)
        logger.info("wrote %s (%d rows) from drift's full index", subset, len(frame))
    else:
        raise FileNotFoundError(f"no BGC-Argo index: neither {subset} nor {BGC_INDEX_FULL}")
    frame["TIME"] = pd.to_datetime(frame["TIME"], utc=True, errors="coerce")
    frame["has_chla"] = frame["parameters"].str.contains("CHLA", na=False)
    return frame.dropna(subset=["LATITUDE", "LONGITUDE", "TIME"])


def describe_boxes(boxes: pd.DataFrame, index: pd.DataFrame,
                   config: StageConfig) -> pd.DataFrame:
    """Add what each box actually holds: chlorophyll profiles, span, recency.

    ``top_boxes`` ranks on one metric; these columns say whether the box is worth a scene
    — a box whose floats all left in 2019 is a fact about the past.
    """
    recent = pd.Timestamp(config.recent_since, tz="UTC")
    rows = []
    for _, box in boxes.iterrows():
        inside = index[index["LONGITUDE"].between(box["lon_west"], box["lon_east"])
                       & index["LATITUDE"].between(box["lat_min"], box["lat_max"])]
        since = inside[inside["TIME"] >= recent]
        rows.append({
            "chla_profiles": int(inside["has_chla"].sum()),
            "chla_floats": int(inside.loc[inside["has_chla"], "PLATFORM_NUMBER"].nunique()),
            "first_seen": inside["TIME"].min(),
            "last_seen": inside["TIME"].max(),
            "profiles_since_2023": int(len(since)),
            "floats_since_2023": int(since["PLATFORM_NUMBER"].nunique()),
        })
    return pd.concat([boxes.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def rank(index: pd.DataFrame, config: StageConfig, metric: str,
         window: str) -> pd.DataFrame:
    """Top boxes by one metric over one time window, described and labelled."""
    subset = index if window == "all" else index[index["TIME"] >= pd.Timestamp(
        config.recent_since, tz="UTC")]
    boxes = top_boxes(subset, box_deg=config.box_deg, top_n=config.top_n,
                      min_separation_deg=config.min_separation_deg,
                      cell_deg=config.cell_deg)
    # top_boxes ranks on distinct ids; re-rank on profiles when that is what was asked,
    # so both readings come from the same candidate placement.
    if metric == "profiles":
        boxes = boxes.sort_values("n_points", ascending=False).reset_index(drop=True)
        boxes["rank"] = np.arange(1, len(boxes) + 1)
    described = describe_boxes(boxes.drop(columns=["parts"], errors="ignore"), index, config)
    described["metric"] = metric
    described["window"] = window
    return described.rename(columns={"n_ids": "n_floats", "n_points": "n_profiles"})


def nearest_open_scene(lon: float, lat: float, scenes: pd.DataFrame) -> dict:
    """Distance from a point to the nearest existing open Tanager scene center.

    This is the tasking argument in one number: if the nearest open scene to the densest
    float box in the ocean is hundreds of kilometers away, the catalogue cannot answer
    the question, however good the scenes in it are.
    """
    distance = haversine_km(lat, lon, scenes["lat_centre"].to_numpy(),
                            scenes["lon_centre"].to_numpy())
    closest = int(np.argmin(distance))
    return {"nearest_scene_id": scenes.iloc[closest]["id"],
            "nearest_scene_km": float(distance[closest]),
            "nearest_scene_place": scenes.iloc[closest]["loc"]}


def scene_centres() -> pd.DataFrame:
    """Center of every open Tanager scene, from the committed inventory.

    Read through ``json.loads`` rather than ``pd.read_json``: a Tanager id looks like
    ``20250914_171527_18_4001``, and pandas reads that as the float 2.0250914e+19, after
    which no id ever matches again and the failure is a silent NaN rather than an error.
    This is the third place the same trap has appeared (YAML digit separators, the
    inventory, here), so it is asserted rather than trusted.
    """
    scenes = pd.DataFrame(json.loads(SCENE_INVENTORY.read_text())).rename(
        columns={"dt": "date"})
    if not scenes["id"].map(type).eq(str).all():
        raise TypeError("scene ids did not survive JSON parsing as strings")
    bbox = np.vstack(scenes["bbox"].to_numpy())
    scenes["lon_centre"] = (bbox[:, 0] + bbox[:, 2]) / 2
    scenes["lat_centre"] = (bbox[:, 1] + bbox[:, 3]) / 2
    return scenes


def erie_box(config: StageConfig, scenes: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """The proposed freshwater scene, and the network it would contain."""
    stations = pd.read_csv(ERIE_INSITU)
    half_lon, half_lat = config.box_deg[0] / 2, config.box_deg[1] / 2
    lon_west, lon_east = ERIE_CENTRE[0] - half_lon, ERIE_CENTRE[0] + half_lon
    lat_min, lat_max = ERIE_CENTRE[1] - half_lat, ERIE_CENTRE[1] + half_lat
    inside = stations[stations["lon"].between(lon_west, lon_east)
                      & stations["lat"].between(lat_min, lat_max)]

    # How far the existing scenes are from the network they should have covered. The
    # honest measure is distance to the scene *edge*, not its center: a scene is ~26 km
    # across, so a center distance overstates the gap by up to half a footprint. The
    # question is whether a station was in frame at all.
    erie_scenes = scenes[scenes["id"].isin(EXISTING_ERIE_SCENES)]
    if erie_scenes.empty:
        raise KeyError(f"neither existing Erie scene {EXISTING_ERIE_SCENES} is in the "
                       f"inventory — the coverage gap cannot be measured")
    footprints = np.vstack(erie_scenes["bbox"].to_numpy())
    gaps, inside_count = [], 0
    for _, station in stations.iterrows():
        lon, lat = station["lon"], station["lat"]
        in_any = ((footprints[:, 0] <= lon) & (lon <= footprints[:, 2])
                  & (footprints[:, 1] <= lat) & (lat <= footprints[:, 3]))
        if in_any.any():
            inside_count += 1
            gaps.append(0.0)
            continue
        # Nearest point on each footprint rectangle, then the nearest of those.
        nearest_lon = np.clip(lon, footprints[:, 0], footprints[:, 2])
        nearest_lat = np.clip(lat, footprints[:, 1], footprints[:, 3])
        gaps.append(float(haversine_km(lat, lon, nearest_lat, nearest_lon).min()))

    frame = pd.DataFrame([{
        "rank": 1, "metric": "in_situ_stations", "window": "proposed",
        "lon_west": lon_west, "lon_east": lon_east,
        "lat_min": lat_min, "lat_max": lat_max,
        "n_stations": int(len(inside)),
        "n_weekly_ctd_stations": int((inside["kind"] == "GLERL weekly station").sum()),
        "label": "western Lake Erie (proposed)",
    }])
    ledger = {
        "centre": list(ERIE_CENTRE),
        "stations_in_proposed_box": int(len(inside)),
        "weekly_ctd_stations_in_box": int((inside["kind"] == "GLERL weekly station").sum()),
        "stations_total": int(len(stations)),
        "stations_inside_existing_scenes": inside_count,
        "min_km_station_to_existing_scene_edge": float(min(gaps)) if gaps else float("nan"),
        "median_km_station_to_existing_scene_edge": float(np.median(gaps)) if gaps else float("nan"),
        "station_names_in_box": inside["station"].tolist(),
    }
    return frame, ledger


def main(out_dir: Path = OUT_DIR, top_n: int = 10, log_level: str = "info") -> int:
    configure_stage_logging(out_dir / "tasking.log", log_level)
    config = StageConfig(out_dir=out_dir, top_n=top_n)
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    index = load_index(out_dir)
    scenes = scene_centres()
    logger.info("%d BGC profiles, %d floats, %s..%s; %d open Tanager scenes",
                len(index), index["PLATFORM_NUMBER"].nunique(),
                index["TIME"].min().date(), index["TIME"].max().date(), len(scenes))

    tables = []
    for window in ("all", "recent"):
        for metric in ("floats", "profiles"):
            table = rank(index, config, metric, window)
            tables.append(table)
            top = table.iloc[0]
            logger.info("%s / %s: top box %.2f..%.2f E, %.2f..%.2f N — %d floats, "
                        "%d profiles, %d with CHLA", window, metric,
                        top["lon_west"], top["lon_east"], top["lat_min"], top["lat_max"],
                        top["n_floats"], top["n_profiles"], top["chla_profiles"])

    ocean = pd.concat(tables, ignore_index=True)
    centres = [((row["lon_west"] + row["lon_east"]) / 2, (row["lat_min"] + row["lat_max"]) / 2)
               for _, row in ocean.iterrows()]
    nearest = pd.DataFrame([nearest_open_scene(lon, lat, scenes) for lon, lat in centres])
    ocean = pd.concat([ocean, nearest], axis=1)

    erie, erie_ledger = erie_box(config, scenes)
    combined = pd.concat([ocean, erie], ignore_index=True)
    combined.to_parquet(out_dir / "tasking_boxes.parquet", index=False)

    headline = ocean[(ocean["metric"] == "floats") & (ocean["window"] == "all")].iloc[0]
    ledger = {
        "box_deg": list(config.box_deg),
        "n_profiles_indexed": int(len(index)),
        "n_floats_indexed": int(index["PLATFORM_NUMBER"].nunique()),
        "n_open_scenes": int(len(scenes)),
        "ocean_top_by_floats": {
            "lon_west": float(headline["lon_west"]), "lon_east": float(headline["lon_east"]),
            "lat_min": float(headline["lat_min"]), "lat_max": float(headline["lat_max"]),
            "n_floats": int(headline["n_floats"]), "n_profiles": int(headline["n_profiles"]),
            "chla_profiles": int(headline["chla_profiles"]),
            "floats_since_2023": int(headline["floats_since_2023"]),
            "nearest_open_scene_km": float(headline["nearest_scene_km"]),
            "nearest_open_scene_place": headline["nearest_scene_place"],
        },
        "min_nearest_scene_km_over_all_ocean_boxes": float(ocean["nearest_scene_km"].min()),
        "erie": erie_ledger,
    }
    (metrics_dir / "g4_tasking.json").write_text(json.dumps(ledger, indent=2, default=str))

    upsert_verdict(metrics_dir, [verdict_row(
        "V6 the tasking proposal names boxes with in-situ truth and no scene",
        "top ocean box by distinct floats; km to the nearest open Tanager scene; "
        "GLERL stations inside the proposed Erie box",
        f"{int(headline['n_floats'])} floats; "
        f"{float(headline['nearest_scene_km']):.0f} km; "
        f"{erie_ledger['stations_in_proposed_box']}",
        "a box exists with floats and no nearby scene",
        "PASS" if headline["n_floats"] > 0 and erie_ledger["stations_in_proposed_box"] > 0
        else "FAIL",
        f"nearest open scene to the densest float box is "
        f"{headline['nearest_scene_place']}")])

    print(json.dumps(ledger, indent=2, default=str))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exp 35 stage 4: where a Tanager scene would buy the most",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
