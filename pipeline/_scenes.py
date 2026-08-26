"""Exp 35 shared wiring: the scene table, and where each scene's files live on disk.

`scenes.yaml` is the single source of truth for which scenes the experiment uses and
which coincident PACE granule and OLCI file belong to each. Every stage reads it through
here rather than re-parsing it, so a scene added or dropped is one edit.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import tanager_io as tio
import yaml

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent
SCENES_YAML = HERE / "scenes.yaml"

#: Downloads and stage outputs. Both default to drift's layout and both are overridable,
#: because these same stage scripts are republished standalone (see `just
#: exp35-publish-study`) where there is no drift tree above them to be relative to.
DATA_DIR = Path(os.environ.get("EXP35_DATA_DIR") or PROJECT_ROOT / "data" / "tanager")
OUT_DIR = Path(os.environ.get("EXP35_OUT_DIR") or HERE / "outputs")

#: Static inputs, script-adjacent like ``scenes.yaml``: they are what the experiment was
#: run against, not something it derived, so they travel with the stages.
SCENE_INVENTORY = HERE / "tanager_open_scenes.json"
ERIE_INSITU = HERE / "erie_insitu.csv"


@dataclass(frozen=True)
class Scene:
    """One scene and everything coincident with it."""

    key: str
    tanager_id: str
    date: str
    tanager_utc: str
    bbox: tuple[float, float, float, float]
    search_bbox: tuple[float, float, float, float]
    pace_aop: str
    pace_bgc: str
    olci_rrs_file: str
    role: str = ""
    cloud_percent: float | None = None

    @property
    def dir(self) -> Path:
        """Where this scene's downloads live. Machine-local; ``data/`` is gitignored."""
        return DATA_DIR / self.key

    @property
    def sr_path(self) -> Path:
        return self.dir / f"{self.tanager_id}_ortho_sr_hdf5.h5"

    @property
    def pace_aop_path(self) -> Path:
        return self.dir / self.pace_aop

    @property
    def pace_bgc_path(self) -> Path:
        return self.dir / self.pace_bgc

    @property
    def olci_path(self) -> Path:
        return self.dir / "olci_l3.nc"


def load_config(path: Path = SCENES_YAML) -> dict:
    """The whole of ``scenes.yaml``."""
    return yaml.safe_load(path.read_text())


def _assert_ids_are_strings(table: dict) -> None:
    """Refuse a scene id that YAML turned into a number.

    Tanager ids look like ``20250514_193937_64_4001``, and YAML 1.1 reads underscores as
    digit separators — so an unquoted id silently becomes the integer
    ``20250514193937644001`` and every URL built from it 404s. This has already happened
    twice (once here, once reading the scene inventory as JSON), so it is a hard check
    rather than a comment.
    """
    numeric = [key for key, entry in table.items() if not isinstance(entry["tanager_id"], str)]
    if numeric:
        raise TypeError(
            f"scene id(s) {numeric} parsed as numbers, not strings — quote them in "
            f"scenes.yaml (YAML reads 20250514_193937_64_4001 as a digit-separated int)")


def load_scenes(path: Path = SCENES_YAML, keys: list[str] | None = None) -> list[Scene]:
    """The scenes, in file order, optionally restricted to ``keys``.

    Raises ``KeyError`` on an unknown key rather than silently returning fewer scenes —
    a typo in ``--scenes`` that quietly runs nothing is the failure this prevents.
    """
    config = load_config(path)
    table = config["scenes"]
    _assert_ids_are_strings(table)
    if keys is not None:
        unknown = [key for key in keys if key not in table]
        if unknown:
            raise KeyError(f"unknown scene(s) {unknown}; have {sorted(table)}")
    wanted = keys if keys is not None else list(table)
    return [Scene(key=key,
                  tanager_id=table[key]["tanager_id"],
                  date=str(table[key]["date"]),
                  tanager_utc=table[key]["tanager_utc"],
                  bbox=tuple(table[key]["bbox"]),
                  search_bbox=tuple(table[key]["search_bbox"]),
                  pace_aop=table[key]["pace_aop"],
                  pace_bgc=table[key]["pace_bgc"],
                  olci_rrs_file=table[key]["olci_rrs_file"],
                  role=" ".join(str(table[key].get("role", "")).split()),
                  cloud_percent=table[key].get("cloud_percent"))
            for key in wanted]


def tanager_asset_url(tanager_id: str, asset: str = "ortho_sr_hdf5",
                      config: dict | None = None) -> str:
    """Direct object-store URL for one Tanager asset, from the configured base.

    Thin wrapper over :func:`tanager_io.asset_url`; the layout rule lives in the library
    and only the base comes from ``scenes.yaml``.
    """
    return tio.asset_url(tanager_id, asset,
                         base=(config or load_config())["tanager_asset_base"])
