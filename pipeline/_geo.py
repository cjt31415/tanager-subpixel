"""
    _geo.py: the two geographic helpers stage 4 needs, by value.

    This is a deliberate copy, not an import. Stage 4 first imported
    ``drift.common.geo.haversine_km`` and ``drift.density.top_boxes``; the submission repo
    assembled from these stages (``just exp35-build-study``) is public and drift is not, so
    that import made stage 4 un-runnable for anyone but us. Libraries are imported by
    reference (``tanager_io``); a competition snapshot copies by value, and this file is
    the copy — drift commit 3b38629, ``src/drift/common/geo.py`` and
    ``src/drift/density.py``, specialised to what stage 4 uses: a global lon/lat grid
    that wraps at the antimeridian, ranked by distinct ids, and a geodesic box area from
    pyproj instead of drift's shapely-lune machinery (the same number to 0.1 km² on a
    Tanager-sized box). Nothing here may import drift.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from pyproj import Geod

logger = logging.getLogger(__name__)

EARTH_RADIUS_KM = 6371.0
WORLD_BOUNDS = (-180.0, -90.0, 180.0, 90.0)
_GEOD = Geod(ellps="WGS84")


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between arrays of points (NaN-safe)."""
    lat1, lon1, lat2, lon2 = (np.asarray(a, dtype=float) for a in (lat1, lon1, lat2, lon2))
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def box_area_km2(lon_west: float, lat_min: float, lon_east: float, lat_max: float,
                 step_deg: float = 0.01) -> float:
    """Geodesic area of a lon/lat box on WGS84, edges densified so they follow parallels.

    A box given as four corners is otherwise measured along geodesics that bow poleward;
    on a Tanager-sized box that is a ~0.1 % error, on a large one it is not.
    """
    if lon_west > lon_east:                       # crosses the antimeridian: two parts
        return (box_area_km2(lon_west, lat_min, 180.0, lat_max, step_deg)
                + box_area_km2(-180.0, lat_min, lon_east, lat_max, step_deg))
    lons = np.arange(lon_west, lon_east, step_deg).tolist() + [lon_east]
    lats = np.arange(lat_min, lat_max, step_deg).tolist() + [lat_max]
    ring_lon = lons + [lon_east] * (len(lats) - 1) + lons[::-1] + [lon_west] * (len(lats) - 1)
    ring_lat = ([lat_min] * len(lons) + lats[1:] + [lat_max] * len(lons) + lats[::-1][1:])
    area_m2, _ = _GEOD.polygon_area_perimeter(ring_lon, ring_lat)
    return abs(area_m2) / 1e6


def _distinct_id_windows(cell_y: np.ndarray, cell_x: np.ndarray, id_codes: np.ndarray, *,
                         n_anchor_y: int, n_lon: int, win_y: int, win_x: int,
                         wrap_lon: bool) -> np.ndarray:
    """``counts[anchor_y, anchor_x]`` = distinct ids with ≥1 point in that window.

    Exact, not a ranked shortlist: ranking windows by point count and then counting ids
    in the top few thousand is biased against exactly the regions this exists to find
    (the Irminger Sea box, 48 floats in 83 profiles, ranks 132,824th of 6.4 M windows by
    profile count). So each id is walked once and every window anchor covering any of
    *its* occupied cells is incremented, deduplicated per id.
    """
    counts = np.zeros(n_anchor_y * n_lon, dtype=np.int64)
    offset_y = np.arange(win_y)
    offset_x = np.arange(win_x)
    order = np.argsort(id_codes, kind="stable")
    codes_sorted = id_codes[order]
    y_sorted, x_sorted = cell_y[order], cell_x[order]
    starts = np.flatnonzero(np.r_[True, codes_sorted[1:] != codes_sorted[:-1]])
    for lo, hi in zip(starts, np.r_[starts[1:], len(codes_sorted)], strict=True):
        cy, cx = y_sorted[lo:hi], x_sorted[lo:hi]
        anchor_y = cy[:, None] - offset_y[None, :]
        anchor_x = cx[:, None] - offset_x[None, :]
        valid_y = (anchor_y >= 0) & (anchor_y < n_anchor_y)
        if wrap_lon:
            anchor_x = anchor_x % n_lon
            valid_x = np.ones_like(anchor_x, dtype=bool)
        else:
            valid_x = anchor_x >= 0
        flat = anchor_y[:, :, None] * n_lon + anchor_x[:, None, :]
        keep = valid_y[:, :, None] & valid_x[:, None, :]
        counts[np.unique(flat[keep])] += 1
    return counts.reshape(n_anchor_y, n_lon)


def _lon_cell_distance(a: int, b: int, n_lon: int, wrap_lon: bool) -> int:
    """Cells between two longitude anchors, the short way round when longitude wraps."""
    straight = abs(int(a) - int(b))
    return min(straight, n_lon - straight) if wrap_lon else straight


def _inside(df: pd.DataFrame, lon_west: float, lat_min: float, lon_east: float,
            lat_max: float, lon_col: str, lat_col: str) -> pd.DataFrame:
    """Rows inside the box, edges inclusive, in one or two parts when it wraps."""
    lon, lat = df[lon_col], df[lat_col]
    in_lat = (lat >= lat_min) & (lat <= lat_max)
    if lon_west <= lon_east:
        return df[in_lat & (lon >= lon_west) & (lon <= lon_east)]
    return df[in_lat & ((lon >= lon_west) | (lon <= lon_east))]


def top_boxes(df: pd.DataFrame, *, box_deg: tuple[float, float], top_n: int = 5,
              min_separation_deg: float = 10.0, cell_deg: float = 0.1,
              lon_col: str = "LONGITUDE", lat_col: str = "LATITUDE",
              id_col: str = "PLATFORM_NUMBER",
              bounds: tuple[float, float, float, float] = WORLD_BOUNDS) -> pd.DataFrame:
    """The ``top_n`` fixed-size boxes holding the most distinct ``id_col`` values.

    ``box_deg`` is ``(lon_extent, lat_extent)`` and must be a whole multiple of
    ``cell_deg``. ``min_separation_deg`` is a floor on how far apart two reported boxes'
    anchors can be — suppressing only overlap returns four tilings of the same blob.
    Returns one row per box, ranked by ``n_ids`` descending, with columns ``rank,
    lon_west, lon_east, lat_min, lat_max, wraps, n_points, n_ids, area_km2``.
    """
    lon_extent, lat_extent = box_deg
    for name, extent in (("lon", lon_extent), ("lat", lat_extent)):
        steps = extent / cell_deg
        if abs(steps - round(steps)) > 1e-9:
            raise ValueError(f"box_deg {name} extent {extent} is not a whole multiple of "
                             f"cell_deg {cell_deg}; the search window could not be that size")
    lon_min, lat_min, lon_max, lat_max = bounds
    n_lon = int(round((lon_max - lon_min) / cell_deg))
    n_lat = int(round((lat_max - lat_min) / cell_deg))
    win_x, win_y = int(round(lon_extent / cell_deg)), int(round(lat_extent / cell_deg))
    wrap_lon = abs((lon_max - lon_min) - 360.0) < 1e-9
    n_anchor_y = n_lat - win_y + 1
    if n_anchor_y < 1 or win_x > n_lon:
        raise ValueError(f"box_deg {box_deg} does not fit inside bounds {bounds}")

    lon_values = df[lon_col].to_numpy(dtype=float)
    lat_values = df[lat_col].to_numpy(dtype=float)
    cell_x = np.clip(((lon_values - lon_min) / cell_deg).astype(int), 0, n_lon - 1)
    cell_y = np.clip(((lat_values - lat_min) / cell_deg).astype(int), 0, n_lat - 1)
    id_codes = pd.factorize(df[id_col])[0]
    score = _distinct_id_windows(cell_y, cell_x, id_codes, n_anchor_y=n_anchor_y,
                                 n_lon=n_lon, win_y=win_y, win_x=win_x, wrap_lon=wrap_lon)

    separation = int(round(min_separation_deg / cell_deg))
    chosen: list[tuple[int, int]] = []
    rows = []
    for flat in np.argsort(score.ravel())[::-1]:
        if len(chosen) == top_n:
            break
        anchor_y, anchor_x = np.unravel_index(flat, score.shape)
        if score[anchor_y, anchor_x] == 0:
            break
        if any(abs(anchor_y - py) < separation
               and _lon_cell_distance(anchor_x, px, n_lon, wrap_lon) < separation
               for py, px in chosen):
            continue
        chosen.append((int(anchor_y), int(anchor_x)))
        west = lon_min + anchor_x * cell_deg
        east = west + lon_extent
        if wrap_lon and east > 180.0:
            east -= 360.0
        south = lat_min + anchor_y * cell_deg
        north = south + lat_extent
        inside = _inside(df, west, south, east, north, lon_col, lat_col)
        rows.append({
            "rank": len(chosen),
            "lon_west": round(west, 4), "lon_east": round(east, 4),
            "lat_min": round(south, 4), "lat_max": round(north, 4),
            "wraps": west > east,
            "n_points": int(len(inside)), "n_ids": int(inside[id_col].nunique()),
            "area_km2": round(box_area_km2(west, south, east, north), 1),
        })
    if len(rows) < top_n:
        logger.warning("top_boxes: asked for %d boxes, found %d — either the points are "
                       "sparse or min_separation_deg=%.1f leaves no room", top_n,
                       len(rows), min_separation_deg)
    return pd.DataFrame(rows)
