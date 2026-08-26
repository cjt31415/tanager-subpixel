# Tanager sub-pixel variance

**How much of the water inside a satellite ocean-colour pixel does the validation
protocol never see?** Measured with Planet Tanager (30 m, hyperspectral) against
same-day PACE OCI and Sentinel-3 OLCI over three scenes.

Submission to the **Planet Tanager Open Data Competition**, 2026.
The write-up is **[MEMO.md](MEMO.md)**, with methods, per-scene results and the full
limitations in **[METHODS.md](METHODS.md)**. This file is the repository's own
documentation.

---

## The finding, in three lines

- **~70 % of the 30 m variance in the water's *colour* is invisible to a 300 m sensor**
  and ~75 % to PACE — an exact law-of-total-variance split, consistent across a turbid
  estuary, a clear-water gulf and a Great Lake.
- Its *brightness* behaves the opposite way (14–18 % hidden at two of three scenes) —
  and brightness is what the standard homogeneity filter screens on, while what a
  match-up certifies is colour-derived. **The filter is screening the one quantity that
  passes its own test.**
- No open Tanager scene sits at any ocean-colour validation site with in-situ truth.
  The nearest open scene to the densest BGC-Argo box is 645 km away; two open scenes
  cover western Lake Erie and miss all 16 NOAA GLERL stations by 10 km.

The satellite-to-float join this bears on is an active one — see
[CHLA-Z](https://fish-pace.github.io/chla-z/), which pairs PACE OCI L3 Rrs with Bio-Argo
chlorophyll profiles, and the [GO-BGC Float Data and Science
Workshop](https://www.go-bgc.org/event/go-bgc-float-data-and-science-workshop).

## What is in here

```
MEMO.md                     the submission write-up
METHODS.md                  how it is measured, the validation gate, full results
docs/DESIGN.md              the pre-registration: gates V0-V6, written before the run
pipeline/                   the six numbered stages, with the inputs they read beside
                            them: scenes.yaml, the 153-scene inventory, the Erie stations
data/figures/               F1 homogeneity · F2 ladder · F3 pedestal · F4 tasking
data/metrics/               verdicts.csv (the gate ledger) + every statistic as JSON
data/*.parquet              per-coarse-pixel aggregates and per-box homogeneity tables
data/tanager/*/olci_l3.nc   the Copernicus OLCI subsets, ~330 kB each
```

Figures sit under `data/` rather than at the top level because that is where a re-run
writes them — so regenerating overwrites exactly what shipped, and a judge can diff.

`data/metrics/verdicts.csv` is the ledger every number in the memo comes from: one row
per gate, with its statistic, its threshold, its verdict, and a note recording which
gates were pre-registered and which were added on revision.

## Reproduce

Stages 3, 4 and 5 run from the data committed here — **no credentials, no downloads**:

```bash
pip install -r requirements.txt
export EXP35_DATA_DIR=$PWD/data/tanager EXP35_OUT_DIR=$PWD/data

python pipeline/compare.py        # stage 3 → gates V3/V4/V5, h_results.json
python pipeline/tasking.py        # stage 4 → tasking_boxes.parquet, g4_tasking.json
python pipeline/plot_figures.py   # stage 5 → figures F1, F3, F4
```

That regenerates every number in the memo's result tables and three of its four figures.
`tanager-io` is the only dependency not on PyPI; `requirements.txt` installs it from
GitHub. Nothing else of ours is imported: the two geographic helpers stage 4 needs
(`pipeline/_geo.py`) are copied in by value, and this exact sequence is run in a fresh
virtualenv before every push.

**What needs more than this repo.** F2 (the resolution ladder) and the variance split
need a Tanager cube — 1.3 GB, public, no account:

```bash
python pipeline/fetch_scenes.py --scenes lake_ontario   # anonymous HTTPS
python pipeline/plot_figures.py
```

Stages 1 and 2 additionally need the PACE granules, which do require a free NASA
Earthdata login in `~/.netrc`. The OLCI subsets are already committed; regenerating them
from source needs a Copernicus Marine account. Both are noted in `pipeline/scenes.yaml`.

Stage 4 ranks every Tanager-footprint-sized box on Earth against the BGC-Argo global
index. The five columns it uses (395,340 profiles; float id, position, time, parameter
list) are committed as `data/bgc_argo_index_subset.parquet` (8 MB), cut from the public
GDAC `argo_bio-profile_index.txt` as of 2026-08-05; its output, `data/tasking_boxes.parquet`,
is what figure F4 draws.

## The library

The reader is released separately as **[tanager-io](https://github.com/cjt31415/tanager-io)**
(MIT, 24 offline tests), because a naive `xr.open_dataset` on a Tanager scene returns an
empty `Dataset`. It absorbs the four things that cost a day: the data group's name
contains a space; the wavelengths live in an attribute, not a coordinate; the
georeferencing hides in a `StructMetadata.0` text blob and the grid is corner-registered;
and fill is −9999 on the floats but 255 on the `uint8` flags.

```python
import tanager_io as tio
ds    = tio.open_sr(tio.asset_url("20250919_170233_04_4001"))
water = tio.water_mask(ds)
rrs   = tio.to_rrs(ds, water=water)
ratio = tio.band_ratio(ds, 560, 490, source=rrs)
```

Every stage in `pipeline/` goes through those calls, so the snippet and the numbers
cannot drift apart. They did once — an earlier version subtracted the atmospheric
pedestal inline and never divided by π, so the gate that was meant to validate Tanager's
variance was silently regressing surface reflectance against Rrs. Nothing downstream
noticed, because every statistic here is a coefficient of variation or a correlation and
both are scale-free. That is why the pipeline now calls the library rather than
reimplementing it.

## Data and licences

| source | licence |
|---|---|
| Planet Tanager open scenes | CC-BY-4.0 |
| NASA PACE OCI L2 (AOP, BGC) | public domain |
| Copernicus Marine `OCEANCOLOUR_GLO_BGC_L3_MY_009_103` | Copernicus licence |
| BGC-Argo global index | CC-BY-4.0 |
| NOAA GLERL/CIGLR station positions | after Boegehold et al. 2023, ESSD 15:3853 |

Code in this repository is MIT (`LICENSE`). The memo and figures are CC-BY-4.0.

## Provenance

Developed as experiment 35 of [drift](https://github.com/cjt31415/drift), where the full
history lives — the brief, the design, and the commit trail. This repository is the
self-contained submission: the stages, the derived data, and the evidence.
