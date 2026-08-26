#!/usr/bin/env python3
"""
    fetch_scenes.py: exp 35 stage 0 — download every scene and its coincident match-ups.

    Input:  scenes.yaml (the scene table); ~/.netrc (NASA Earthdata);
            ~/.copernicusmarine (Copernicus Marine).
    Output: data/tanager/<scene>/<id>_ortho_sr_hdf5.h5   Tanager surface reflectance
            data/tanager/<scene>/PACE_OCI.*.L2.OC_{AOP,BGC}.*.nc
            data/tanager/<scene>/olci_l3.nc              CHL + RRS, cut to the scene bbox
            outputs/metrics/g0_fetch.json                sizes, times and what was skipped

    NETWORK, ~1.5 GB per scene, ~10 min on a fast link. Resumable: a file whose size
    already matches the server's Content-Length is skipped, and a partial download
    continues with an HTTP Range request rather than starting over.

    The two OLCI paths are different on purpose (DESIGN §13.1): chlorophyll comes from
    Copernicus Marine's ARCO subset service in seconds, but the reflectance product has
    no subset service, so its 5.75 GB global daily file is range-read for the scene's
    bounding box — about 2 s per band against 17 GB of download avoided.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import requests
import xarray as xr
from _scenes import OUT_DIR, Scene, load_config, load_scenes, tanager_asset_url
from _verdicts import configure_stage_logging, upsert_verdict, verdict_row

logger = logging.getLogger(__name__)

#: Streaming chunk for the big Tanager files.
CHUNK_BYTES = 1 << 20

#: (connect, read) timeouts. The object store is fast; a stall means a dead connection.
TIMEOUTS = (10, 120)

#: A scene is useless below this many valid OLCI chlorophyll pixels in its box.
MIN_OLCI_PIXELS = 500


@dataclass
class StageConfig:
    """Stage 0 wiring."""

    out_dir: Path = OUT_DIR
    skip_pace: bool = False
    skip_olci: bool = False


def download(url: str, dest: Path) -> dict[str, object]:
    """Fetch ``url`` to ``dest``, resuming a partial file, skipping a complete one.

    Returns a small ledger entry. Writes through ``<dest>.part`` and renames on success,
    so a killed run can never leave a truncated file that looks finished.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    head = requests.head(url, allow_redirects=True, timeout=TIMEOUTS)
    head.raise_for_status()
    expected = int(head.headers.get("Content-Length", 0))

    if dest.exists() and expected and dest.stat().st_size == expected:
        logger.info("have %s (%.2f GB)", dest.name, expected / 1e9)
        return {"file": dest.name, "bytes": expected, "action": "skipped"}

    part = dest.with_suffix(dest.suffix + ".part")
    have = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    logger.info("fetching %s (%.2f GB)%s", dest.name, expected / 1e9,
                f", resuming at {have / 1e9:.2f} GB" if have else "")
    with requests.get(url, stream=True, headers=headers, timeout=TIMEOUTS) as response:
        response.raise_for_status()
        mode = "ab" if have and response.status_code == 206 else "wb"
        with open(part, mode) as handle:
            for block in response.iter_content(chunk_size=CHUNK_BYTES):
                handle.write(block)
    os.replace(part, dest)
    size = dest.stat().st_size
    if expected and size != expected:
        raise OSError(f"{dest.name}: got {size} bytes, server said {expected}")
    return {"file": dest.name, "bytes": size, "action": "downloaded"}


def fetch_tanager(scene: Scene, config: dict) -> dict[str, object]:
    """The scene's ortho surface-reflectance cube."""
    return download(tanager_asset_url(scene.tanager_id, "ortho_sr_hdf5", config), scene.sr_path)


def fetch_pace(scene: Scene) -> list[dict[str, object]]:
    """The coincident PACE OCI L2 apparent-optical-property and biogeochemistry granules.

    Granule names are pinned in ``scenes.yaml`` rather than re-searched, so a rerun
    cannot silently pick a different overpass than the one the design measured.
    """
    import earthaccess

    earthaccess.login(strategy="netrc")
    ledger: list[dict[str, object]] = []
    for name, dest in ((scene.pace_aop, scene.pace_aop_path),
                       (scene.pace_bgc, scene.pace_bgc_path)):
        if dest.exists() and dest.stat().st_size > 0:
            logger.info("have %s", dest.name)
            ledger.append({"file": dest.name, "bytes": dest.stat().st_size, "action": "skipped"})
            continue
        short_name = "PACE_OCI_L2_AOP" if "OC_AOP" in name else "PACE_OCI_L2_BGC"
        granules = earthaccess.search_data(short_name=short_name, temporal=(scene.date, scene.date),
                                           bounding_box=scene.search_bbox)
        match = [g for g in granules if g["umm"]["GranuleUR"] == name]
        if not match:
            raise FileNotFoundError(f"{name} not found in CMR for {scene.date} "
                                    f"{scene.search_bbox}; found "
                                    f"{[g['umm']['GranuleUR'] for g in granules]}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        earthaccess.download(match[:1], str(dest.parent))
        ledger.append({"file": dest.name, "bytes": dest.stat().st_size, "action": "downloaded"})
    return ledger


def fetch_olci(scene: Scene, config: dict) -> dict[str, object]:
    """OLCI chlorophyll (ARCO subset) and reflectance (range-read), cut to the scene box.

    Written as one netCDF per scene so no later stage ever touches the 5.75 GB global
    file again.
    """
    if scene.olci_path.exists() and scene.olci_path.stat().st_size > 0:
        logger.info("have %s", scene.olci_path.name)
        with xr.open_dataset(scene.olci_path) as existing:
            valid = int(np.isfinite(existing["CHL"]).sum())
        return {"file": scene.olci_path.name, "action": "skipped", "chl_valid_px": valid}

    import copernicusmarine as cm

    lon_min, lat_min, lon_max, lat_max = scene.bbox
    pad = 0.05
    logger.info("OLCI chlorophyll (ARCO subset) for %s", scene.key)
    chl = cm.open_dataset(dataset_id=config["olci"]["chl_dataset_id"],
                          minimum_longitude=lon_min - pad, maximum_longitude=lon_max + pad,
                          minimum_latitude=lat_min - pad, maximum_latitude=lat_max + pad,
                          start_datetime=scene.date, end_datetime=scene.date).load()

    url = f"{config['olci']['https_base']}/{scene.olci_rrs_file}"
    bands = [f"RRS{nm}" for nm in config["olci"]["rrs_bands_nm"]]
    logger.info("OLCI reflectance (range-read of the global daily file), %d bands", len(bands))
    import fsspec

    handle = fsspec.open(url, "rb", block_size=8 * 2**20).open()
    with xr.open_dataset(handle, engine="h5netcdf") as full:
        keep = [name for name in bands + [f"{b}_uncertainty" for b in bands] if name in full]
        rrs = full[keep].sel(lon=slice(lon_min - pad, lon_max + pad),
                             lat=slice(lat_max + pad, lat_min - pad)).load()

    merged = xr.merge([chl, rrs], compat="override", join="override")
    merged.attrs.update({"scene": scene.key, "date": scene.date,
                         "chl_dataset_id": config["olci"]["chl_dataset_id"],
                         "rrs_dataset_id": config["olci"]["rrs_dataset_id"],
                         "note": "CMEMS OLCI L3 300 m, cut to the Tanager scene bbox + 0.05 deg"})
    scene.olci_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = scene.olci_path.with_suffix(".nc.part")
    merged.to_netcdf(temporary)
    os.replace(temporary, scene.olci_path)
    valid = int(np.isfinite(merged["CHL"]).sum())
    return {"file": scene.olci_path.name, "action": "downloaded", "chl_valid_px": valid,
            "rrs_bands": len(keep)}


def gate(ledger: dict) -> list[str]:
    """G0: every scene has every file, each opens, and OLCI has usable coverage."""
    problems: list[str] = []
    for key, entry in ledger["scenes"].items():
        if not entry.get("tanager"):
            problems.append(f"{key}: no Tanager cube")
        if "olci" in entry:
            valid = entry["olci"].get("chl_valid_px", 0)
            if valid < MIN_OLCI_PIXELS:
                problems.append(f"{key}: only {valid} valid OLCI chlorophyll pixels "
                                f"(floor {MIN_OLCI_PIXELS})")
        for problem in entry.get("problems", []):
            problems.append(f"{key}: {problem}")
    return problems


def check_opens(scene: Scene, skip_pace: bool, skip_olci: bool) -> list[str]:
    """Open every downloaded file for this scene; return what failed."""
    problems: list[str] = []
    try:
        import tanager_io as tio
        dataset = tio.open_sr(scene.sr_path)
        n_bands = dataset.sizes["wavelength"]
        if n_bands != 426:
            problems.append(f"Tanager cube has {n_bands} bands, expected 426")
        wavelength = dataset["wavelength"].values
        if not np.all(np.diff(wavelength) > 0):
            problems.append("Tanager wavelengths are not strictly increasing")
        dataset.close()
    except Exception as error:  # noqa: BLE001 - the gate reports, it does not diagnose
        problems.append(f"Tanager cube will not open: {error}")
    if not skip_pace:
        for path in (scene.pace_aop_path, scene.pace_bgc_path):
            try:
                with xr.open_dataset(path, group="geophysical_data") as pace:
                    if not len(pace.data_vars):
                        problems.append(f"{path.name}: no geophysical_data variables")
            except Exception as error:  # noqa: BLE001
                problems.append(f"{path.name} will not open: {error}")
    if not skip_olci:
        try:
            with xr.open_dataset(scene.olci_path) as olci:
                if "CHL" not in olci:
                    problems.append("OLCI subset has no CHL")
        except Exception as error:  # noqa: BLE001
            problems.append(f"olci_l3.nc will not open: {error}")
    return problems


def fetch_all(scenes: list[Scene], config: StageConfig) -> dict:
    """Download everything for every scene and build the gate ledger."""
    scene_config = load_config()
    ledger: dict = {"scenes": {}}
    for scene in scenes:
        logger.info("=== %s (%s, %s) ===", scene.key, scene.tanager_id, scene.date)
        entry: dict = {"tanager_id": scene.tanager_id, "date": scene.date}
        entry["tanager"] = fetch_tanager(scene, scene_config)
        if not config.skip_pace:
            entry["pace"] = fetch_pace(scene)
        if not config.skip_olci:
            entry["olci"] = fetch_olci(scene, scene_config)
        entry["problems"] = check_opens(scene, config.skip_pace, config.skip_olci)
        ledger["scenes"][scene.key] = entry
    return ledger


def main(scenes: list[str] | None = None, out_dir: Path = OUT_DIR, skip_pace: bool = False,
         skip_olci: bool = False, log_level: str = "info") -> int:
    configure_stage_logging(out_dir / "fetch_scenes.log", log_level)
    config = StageConfig(out_dir=out_dir, skip_pace=skip_pace, skip_olci=skip_olci)
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    ledger = fetch_all(load_scenes(keys=scenes), config)
    problems = gate(ledger)
    ledger["problems"] = problems
    (metrics_dir / "g0_fetch.json").write_text(json.dumps(ledger, indent=2, default=str))

    upsert_verdict(metrics_dir, [verdict_row(
        "G0 every scene has a Tanager cube and coincident match-ups that open",
        "scenes complete; problems",
        f"{len(ledger['scenes'])}; {len(problems)}",
        "all scenes; 0 problems",
        "PASS" if not problems else "FAIL",
        "; ".join(problems) or f"scenes: {', '.join(ledger['scenes'])}")])

    print(json.dumps({"scenes": list(ledger["scenes"]), "problems": problems}, indent=2))
    if problems:
        raise RuntimeError("G0 failed: " + "; ".join(problems))
    return 0


def parse_opt() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exp 35 stage 0: download Tanager scenes and their PACE/OLCI match-ups",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenes", nargs="*", default=None,
                        help="scene keys from scenes.yaml (default: all)")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--skip-pace", action="store_true")
    parser.add_argument("--skip-olci", action="store_true")
    parser.add_argument("--log-level", choices=["debug", "info", "warning", "error"],
                        default="info")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main(**vars(parse_opt())))
