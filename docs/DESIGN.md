# Experiment 35 — DESIGN: Tanager ⇄ OLCI ⇄ PACE match-up and tasking proposal

*Written 2026-08-25 against `README.md` (the brief). Deadline 2026-08-31. Five working
days. Every stage below has a gate; a stage that fails its gate stops the run and is
reported, not worked around. Decisions superseded during the build go in §12.*

## 0. Deliverable and how we will know it worked

**Deliverable** (formats and scoring from the T&C, `docs/Planet-TermsConditions-TanagerCompetition.pdf`, read 2026-08-25 — §13.2):
1. `outputs/case_study.md` — the 1–3 page **technical memo**, with the two sections the
   rubric names explicitly: a *Project Summary* and an *Impact Statement*.
2. `outputs/tanager_matchup.ipynb` — runs top-to-bottom from the stage outputs, no
   network; pushed with the stage scripts to a **public repo** (reproducibility is 20 pts).
3. A pip-installable `tanager_io` (the `_io.py` reader + a `pystac` STAC walker) as the
   open-source component — the "reduce the technical barrier to using Tanager data"
   line in the rubric, and the +5 tie-breaker.
4. Derived parquets (agg tables, tasking boxes) on Zenodo — optional, listed as accepted.
Everything the memo and notebook show is computed by the numbered scripts; the notebook
only reads parquet and PNGs. Submission is one entry via the SurveyMonkey link, by
2026-08-31 23:59 PST.

**It worked if** all of these hold, in order, and the notebook says which did not:

| # | check | pass |
|---|---|---|
| V1 | Golden Gate: ≥ 200 PACE pixels and ≥ 2,000 OLCI pixels that are ≥ 90 % Tanager-water, cloud-free in all three | numbers printed in stage 2 |
| V2 | Spectral registration: aggregated-Tanager vs PACE Rrs, normalised at 555 nm, median per-band ratio within 0.85–1.15 over 443–670 nm | stage 3 table |
| V3 | Sign gate ([[gate-signs-not-just-magnitudes]]): in the Golden Gate plume, Tanager sub-pixel chl-a CV is *higher* inside the plume front than in the bay interior, and the Loreto CV is lower than both — direction asserted, not just magnitude | stage 3 assertion |
| V4 | H3 answered either way: the OLCI↔PACE chl-a agreement (RMSD, log space) and the Tanager↔each agreement are both reported with N | stage 3 table |
| V5 | Tasking map renders with the Ligurian and Erie boxes and their in-situ counts | stage 4 PNG |

## 1. Evaluation of the brief — what changes and why

1. **Tanager SR, not radiance.** The brief hedged on whether SR existed; it does
   (`ortho_sr_hdf5`, read 2026-08-25 by HTTP range request). Layout: HDF-EOS5 grid
   `HDFEOS/GRIDS/HYP/Data Fields/`, EPSG 32610 for Golden Gate (UTM zone in the group
   attrs per scene), `surface_reflectance (426, rows, cols)` float32 with `fwhm` per band
   (5.2–6.8 nm), `surface_reflectance_uncertainty`, `aerosol_optical_depth`,
   `column_water_vapour`, `beta_cloud_mask`, `beta_cirrus_mask`, `nodata_pixels`,
   `sensor_zenith/azimuth`, `sun_zenith/azimuth`, per-pixel `time`. Golden Gate is
   867 × 852 px ≈ 26 × 26 km. Fill is −9999.
2. **There is no water mask in the product.** Planet's SR is land-tuned and ships no
   land/water field. Stage 1 builds one from SWIR (water is ~0 at 1,600 nm) plus the
   cloud/cirrus/nodata masks. This is a **probe candidate** (surfaced once, here): the
   SWIR threshold is a by-eye choice over a plume that is bright in the visible; the
   stage prints the histogram and the chosen cut, and the notebook shows the mask.
3. **Tanager SR → Rrs is an approximation, and is said so.** SR is bottom-of-atmosphere
   reflectance including sky glint and no BRDF; Rrs ≈ SR/π after subtracting a
   near-infrared residual (mean SR over 860–880 nm) as a glint proxy. PACE `Rrs` and
   OLCI `ρw/π` are proper water-leaving quantities. H1 is therefore a *shape*
   comparison (normalised at 555 nm); absolute level is reported but not scored.
4. **One chl-a algorithm across all three sensors, at the coarse sensor's bands.** The
   sub-pixel CV (H2) must not be contaminated by algorithm differences, so Tanager is
   band-resampled to OLCI (Gaussian SRF from published OLCI centre/width) and to PACE
   (5 nm — direct linear interpolation), then the same OC4-form blue/green ratio is applied
   to Tanager-resampled, OLCI and PACE alike. The sensors' own products (`chlor_a`,
   `CHL_OC4ME`, `CHL_NN`) are compared *separately* as "what the user would get".
5. **Coarse pixel footprints from the L2 lat/lon, not from a reprojection.** PACE and OLCI
   L2 are swath products. Each coarse pixel becomes a polygon from the midpoints to its
   four neighbours (in UTM), and Tanager water pixels are binned into it. Edge rows lose
   their footprint and are dropped. Reported as N per pixel.
6. **Erie gets a cyanobacteria index, not just chl-a.** Tanager at 5 nm resolves the
   phycocyanin absorption at 620 nm; OLCI has 620 nm; PACE has it too. Stage 3 computes the
   Wynne cyanobacteria index (CI, 665/681/709) on all three and a 620 nm depth on the two
   hyperspectral ones. GLERL station values are **not** joined — the brief records that
   the week's sampling could not be read; the notebook states the stations' distance.
7. **Scope cut decided now, not on Friday.** Golden Gate is the full triple. Loreto and
   Erie are Tanager ↔ PACE ↔ OLCI only if stage 0 has all three files by end of day 1;
   otherwise Tanager ↔ PACE. The tasking section (stage 4) is independent of stages 1–3
   and is built on day 1 so it exists whatever happens.

## 2. Layout

```
experiments/35_tanager/
├── README.md, DESIGN.md
├── tanager_open_scenes.json        # 153-scene inventory (committed)
├── scenes.yaml                     # the three scenes: ids, dates, bboxes, PACE/OLCI granule ids
├── _io.py                          # open_tanager(), open_pace_l2(), open_olci_wfr(): xarray, one call each
├── _spectral.py                    # band resampling, Rrs approx, OC4-form ratio, CI, 620 depth
├── fetch_scenes.py                 # stage 0
├── build_water_mask.py             # stage 1
├── aggregate.py                    # stage 2
├── compare.py                      # stage 3
├── tasking.py                      # stage 4
├── plot_figures.py                 # stage 5
└── outputs/                        # parquet, png, notebook, case_study.md (gitignored except .ipynb, .md)
data/tanager/<scene_id>/            # SR h5, UDM, PACE nc, OLCI L3 subset nc — machine-local, gitignored (data/* already is)
```

No new credentials: OLCI comes from Copernicus Marine via the existing
`~/.copernicusmarine` login (see §13.1). `justfile` gains `exp35-fetch`, `exp35-mask`, `exp35-aggregate`,
`exp35-compare`, `exp35-tasking`, `exp35-figures`, `exp35` (all, in order).

Conda env `argo`: `h5py`, `fsspec`, `earthaccess`, `xarray`, `netCDF4`, `pyproj`,
`shapely`, `geopandas` are already present (checked by import on 2026-08-25). No new
dependencies.

## 3. Stage 0 — `fetch_scenes.py` (G0)

For each scene in `scenes.yaml`: Tanager `ortho_sr_hdf5` + `ortho_beta_udm` (public GCS,
`requests` streaming, resume on partial); PACE `PACE_OCI_L2_AOP` + `_L2_BGC` granule by
id via `earthaccess` (netrc); OLCI as a **CMEMS L3 300 m subset** (§13.1): CHL + flags by
`copernicusmarine.open_dataset` (ARCO), RRS400…RRS709 + uncertainties by HTTP range-read
of the global daily file, both cut to the scene bbox and written as one netCDF per scene.
Sizes logged; nothing re-downloaded if the file exists and its size
matches the header.

**G0:** every file opens (`_io.py`), Tanager has 426 bands with monotone wavelengths,
PACE `Rrs` has a `wavelength_3d` axis, the OLCI subset has the 11 `RRS*` bands and `CHL`
with > 1,000 valid pixels. Print the per-scene time offsets (Tanager − PACE) from the
files, not the STAC; for OLCI record the S3A/S3B pass times from the CDSE catalogue in
`scenes.yaml` (L3 merges the passes and keeps no time).
Fails loudly ([[loud-failures-over-quiet-noops]]) on any missing file; a scene can be
excluded with `--scenes`.

## 4. Stage 1 — `build_water_mask.py` (G1)

Per scene, on the Tanager grid: `valid = ~nodata & ~cloud & ~cirrus`; `water = valid &
SR(1600 nm) < t` with `t` chosen from the valley of the 1,600 nm histogram (Otsu), printed,
overridable with `--swir-threshold`. Morphological open (3 px) to drop boats and wakes.
Output `outputs/<scene>_water_mask.tif` (COG, same grid) and a PNG of RGB with the mask
outline.

**G1:** water fraction of valid pixels > 0.30 (Golden Gate, Loreto) or > 0.50 (Erie);
otherwise stop — the scene is not usable for a water comparison and the notebook says so.

## 5. Stage 2 — `aggregate.py` (G2)

For each scene and each coarse sensor: build coarse-pixel polygons (§1 item 5), rasterise
their ids onto the Tanager grid, and for each id compute `n_water`, `n_valid`,
`frac_water`, mean and std of Tanager SR across all 426 bands over water pixels, plus the
mean Tanager `time`. Write `outputs/<scene>_<sensor>_agg.parquet`, one row per coarse
pixel, columns `rrs_mean_<λ>` / `rrs_std_<λ>` wide (426 × 2 float32 — ~1 MB per 300 rows;
fine). Keep the coarse sensor's own values on the same row (`Rrs_<λ>` / `Oa##`, product
chl-a, flags).

**G2:** V1 counts. Also a registration check: the Tanager-aggregated 865 nm over water is
near zero and does not correlate with `frac_water` (if it does, the mask leaks land).

## 6. Stage 3 — `compare.py` (H1, H2, H3, Erie CI)

Reads the agg parquets only.

- **H1** — Tanager SR → Rrs (§1 item 3), band-resampled to PACE wavelengths and to OLCI
  bands; per band: median ratio, MAD, N, on pixels with `frac_water ≥ 0.9` and coarse
  flags clean. Normalised-at-555 shape comparison plus the absolute one, both tabulated.
  A wavelength-dependent bias that is smooth is atmospheric correction; one that is
  spiky at 760 nm (O₂ A) or 940 nm (H₂O) is a gas-correction residual — the table flags
  both by simple rules and the notebook says which it sees.
- **H2** — OC4-form ratio chl-a on (a) every Tanager water pixel, resampled to the coarse
  bands, (b) Tanager mean spectrum per coarse pixel, (c) the coarse sensor itself. Per
  coarse pixel: CV and skewness of (a); difference (b) − mean(a) (Jensen gap of the ratio
  algorithm — how wrong is "algorithm on the mean"). Maps of CV; V3 sign assertion by
  comparing the plume polygon (hand-drawn once in `scenes.yaml` from the RGB) to the
  bay-interior polygon, and Golden Gate to Loreto.
- **H3** — log-RMSD and bias: OLCI↔PACE, Tanager↔PACE, Tanager↔OLCI, on the common
  pixels (OLCI pixels aggregated into PACE footprints for the first pair). V4.
- **Erie** — CI and 620 nm depth on all three; map and histogram; CV of CI inside the
  scene's bloom polygon.

Outputs: `outputs/h1_spectral.parquet`, `outputs/h2_subpixel.parquet`,
`outputs/h3_basis.parquet`, `outputs/erie_ci.parquet`, and `outputs/verdicts.json`
(V1–V4 with the numbers that decided them).

## 7. Stage 4 — `tasking.py` (V5)

Ocean: run exp 20's `hotspots.py` machinery with `--box-deg 0.30 0.22`, top 10 by
distinct floats and top 10 by profiles, all years and last 24 months, adding a
`since_2023` column and a `chla_profiles` column. Freshwater: the Erie box centred at
41.80 N, −83.28 with the in-situ inventory from README §4b (a small table committed as
`erie_insitu.csv`, coordinates from Boegehold et al. 2023 Table 1 and the NDBC station
table). Output `outputs/tasking_boxes.parquet` and `outputs/tasking_map.png` (world
inset + Ligurian and Erie panels with station/float markers).

## 8. Stage 5 — `plot_figures.py` and the notebook

Figures: (F1) three-panel RGB + water mask; (F2) resolution ladder — the same 3 km of the
Golden Gate front at 30 m / 300 m / 1.2 km; (F3) H1 per-band ratio with MAD envelope;
(F4) CV maps and the plume-vs-interior box plot; (F5) H3 scatter triplet; (F6) Erie CI
triplet; (F7) tasking map. The notebook is assembled last, reads only `outputs/`, and
prints `verdicts.json` at the top.

## 9. Tests — `tests/test_exp35_tanager.py` (synthetic, no network)

- Band resampling: a synthetic 5 nm spectrum convolved with a Gaussian SRF returns the
  analytic value; flat spectrum returns itself for every sensor.
- Aggregation: a synthetic 30 m field of known block means, with a fake coarse lat/lon
  grid, recovers those means and the correct `n_water`; a masked block is excluded.
- Water mask: a synthetic SWIR image with two modes picks the valley; the override wins.
- OC4-form ratio: monotone in the ratio; Jensen gap is zero on a constant field and
  positive on a bimodal one (sign gate for H2's machinery).

## 10. Documentation and commits

`[PLAN]` this design. `[IMPL]` after stage 0–2 run on Golden Gate (day 2). `[IMPL]` after
stage 3–5 (day 4). `docs/WORKFLOWS.md` gains an "Exp 35" block (fetch needs CDSE creds,
~4 GB per scene, run on harpy). `docs/DECISIONS.md` entry on the day the SR→Rrs
approximation is fixed (§1 item 3) — it is the decision most likely to be questioned.
`docs/TODO.md` "Probe candidates" gets the SWIR water-mask threshold (§1 item 2) so it is
never re-raised. Submission itself happens outside the repo; the notebook and
`case_study.md` are the committed record.

## 11. Known risks, with the response already decided

| risk | response |
|---|---|
| CMEMS range-read of the 5.75 GB global reflectance file stalls | Measured 2026-08-25: 4.7 s open, 2.3 s per band for the Golden Gate box; if it regresses, `copernicusmarine get` the whole file once (×3 dates = 17 GB) and subset locally |
| Planet SR over water is unusable (negative or wildly glinty in the blue) | Fall back to `ortho_radiance_hdf5` and compare *top-of-atmosphere* Tanager against PACE `rhot` (L1B not fetched — use PACE L2 AOP's `Rrs` only if SR works). Documented as the finding it is |
| PACE 1.2 km pixels at scene edges have too few Tanager water pixels | `frac_water ≥ 0.9` filter already excludes them; V1 states the count |
| Golden Gate 2025-05-14 haze (`light_haze_percent` 29) | AOD field is in the file; H1 stratified by AOD tercile if the residual correlates with it |
| Sub-pixel CV is dominated by the water mask edge, not the plume | The 3 px open plus `frac_water ≥ 0.9`; V3 uses interior polygons only |
| Time: 5 days | Stage 4 on day 1; stages 0–2 on Golden Gate by day 2; if stage 3 slips, F2 + F4 + F7 alone are a submission |

## 11b. Scoring map (100 + 15)

| rubric | pts | where this submission earns it |
|---|---|---|
| Scientific integrity & innovation | 30 | V1–V4 with N; limitations section (SR→Rrs approximation, AC differences, L3 merge); novelty = Tanager as sub-pixel lens + float-dense tasking |
| Application / use case | 30 | Relevance: pixel-scale representativeness for float/mooring cal-val (GO-BGC, PACE validation); feasibility: the tasking map is an implementable decision for Planet's Open Committee |
| Workflow & tool | 20 | `tanager_io` (STAC walk with `pystac`, HDF-EOS reader, water mask), public repo, `just exp35` reproduces |
| Visualization & storytelling | 20 | F1–F7; memo Project Summary + Impact Statement |
| +5 strategic area | 5 | water quality (Erie HAB, Golden Gate plume) |
| +5 comparison vs public alternatives | 5 | PACE & OLCI same-day (H1–H3); EMIT sensor-level comparison if day 4 allows (§13.2) |
| +5 open-source / AI-ML | 5 | `tanager_io` |

## 12. Open decisions

- Whether to present chl-a at all for Golden Gate (turbid, CDOM-rich; OC4 is wrong there)
  or to present only the reflectance-shape and CV results and keep chl-a for Loreto. Decide
  after seeing F3 on day 2; default is to show both and say OC4 is a *consistency* metric,
  not a chl-a claim.
- Which PACE pixel geometry to trust: L2 `latitude/longitude` centres with midpoint
  polygons (chosen) vs the pixel-corner fields if the granule carries them (check in G0).
- Whether OLCI 300 m is aggregated into PACE for H3 (chosen) or PACE is nearest-neighbour
  sampled at OLCI pixels (more rows, wrong scale).

## 13. What the build changed — superseded decisions

### 13.8 The pedestal is a property of the scene, not the product (corrects §13.4)

§13.4 generalised from one scene and was wrong to. Measured over water on all three:

| scene | Planet `light_haze_percent` | 443 nm | 560 nm | 865 nm | 1600 nm |
|---|---|---|---|---|---|
| Lake Ontario | 0 | 0.0276 | 0.0235 | 0.0026 | **0.0025** |
| Loreto | – | 0.0612 | 0.0274 | 0.0113 | **0.0090** |
| San Francisco Bay | 29 | 0.1220 | 0.1458 | 0.1122 | **0.0971** |

Physical water is under 0.01 at 1600 nm. **Lake Ontario is physical; Loreto is close;
only San Francisco Bay carries the ~0.09 pedestal — and it is the scene Planet flags at
29 % light haze.** So the residual tracks the atmosphere on the day, not the product, and
the diagnostic is free: look at 1600 nm over water.

Three corrections follow:

1. The NDWI mask is still the right default, but the reason changes — not "the product
   always has a pedestal" but "NDWI survives both cases, and a SWIR threshold does not".
2. Absolute reflectance is usable in a clear scene. §1 item 3 and §13.4 said it never is.
3. **The lower-bound argument is scene-specific.** A pedestal deflates a CV, so San
   Francisco Bay's within-pixel CV is the most suppressed of the three — which is a
   confound on reading it as "the control that behaves the other way". Independent
   support that its water really is uniform at this scale comes from OLCI, which has no
   pedestal problem and also reports a low *between*-pixel CV there (0.029). The headline
   effect is cleanest at Lake Ontario and Loreto, the two near-physical scenes, and the
   memo must present it that way.

### 13.7 The CV needs a physical validity range, and V3 is per pairing (refines §6, §0)

A coefficient of variation diverges as its mean approaches zero. A handful of coarse
pixels whose mean green/blue ratio sat near zero — failed pedestal removal, not water —
dominated every upper percentile (lake_ontario/olci reported a p90 of 4.5). Coarse pixels
whose mean ratio falls outside [0.2, 5.0] are now excluded, in the headline table and the
structure score alike, so the two cannot diverge.

V3 is recorded per scene and sensor rather than as one worst-case gate, because it
licenses the interpretation of a given pairing. Result: structured on 4 of 5 assessable
pairings (0.17-0.44 against a shuffled field's 1.0), **loreto/pace marginal at 0.79**, and
sf_bay/pace not assessable at 19 adjacent finite pairs. The overall verdict is MIXED and
the memo says so. Loreto is the clear-water control, so less spatial organisation at PACE
scale is consistent with the science rather than a defect — but it is reported, not
explained away.

### 13.6 V1 tests containment, not equality (fixes §0 V1)

The ortho grid is a north-up array sized to hold a *rotated* imaging strip, so its extent
is the bounding box of the array while the STAC bbox is the bounding box of the imaged
data inside it. The grid is legitimately larger, by the nodata margin — 1.35 px on sf_bay,
2.94 px on loreto. Equality within a tolerance would have produced false failures on the
more strongly rotated scenes. The gate now requires the grid to *contain* the STAC bbox
(shortfall 0.00 px on all three) and reports the margin without failing on it.

### 13.5 The CMEMS "olci-300m" product is not 300 m (refines §1 item 5, §2)

`cmems_obs-oc_glo_bgc-*_my_l3-olci-300m_P1D` is mapped to a 1/180 deg grid: **488 x 618 m
at 38N**, about 335 Tanager pixels per cell, measured from the delivered file. The "300m"
in the identifier refers to OLCI's full-resolution source, not the mapped product. The
resolution ladder is therefore **30 m -> ~550 m -> ~1.2 km**, and the memo must say so.

### 13.4 Tanager surface reflectance over water carries a large additive pedestal (supersedes §1 items 2 and 3)

Measured, not assumed. A 16x16 box of deep open water in central San Francisco Bay reads
0.137 at 441 nm, 0.163 at 561 nm, 0.107 at 866 nm and **0.093 at 1598 nm**, where physical
water is under 0.01. What this is: real water-leaving signal — the green peak above the
blue is there — on an additive pedestal of ~0.08-0.10 rising towards the blue, i.e.
uncorrected Rayleigh and aerosol path reflectance that a land-tuned correction leaves over
a dark target. Three consequences:

1. **The SWIR water mask of §4 does not work.** The pedestal lifts water at 1598 nm to
   just under vegetation (0.110); Otsu on that histogram split bright land from dark land
   and returned a "44% water" mask whose clearest decile had a vegetation red edge. The
   mask is now **NDWI** (green vs NIR), where water sits at +0.21 and land at -0.18
   (urban) to -0.72 (vegetation) because the pedestal is common to both terms. Water
   fractions: sf_bay 33.9%, lake_ontario 48.4%, loreto 43.8%.
2. **Absolute Rrs is not recoverable**, so H1 is a shape and correlation comparison only —
   as §1 item 3 anticipated, now with the measurement behind it.
3. **Every CV is a lower bound.** An additive term inflates a mean and leaves a standard
   deviation alone, so it deflates every coefficient of variation. The headline finding is
   therefore conservative, which is the useful direction to be wrong in.

### 13.3 The scene set changed on the evidence of the thumbnails (supersedes §7, README §2)

Every planned scene was inspected before any download. Two were wrong:

* **"Golden Gate" contains no Golden Gate.** The bridge is at -122.478; the scene starts
  at -122.448, 2.7 km east. The scene is central San Francisco Bay and the
  Berkeley/Richmond shoals — still an excellent subject, but renamed `sf_bay`, and the
  memo must not describe a plume through the Gate.
* **The Lake Erie pair is unusable.** The north scene (`20250914_171527_18_4001`) is ~90%
  downtown Detroit with only the river as water; the south is a river mouth at ~25% water,
  and the 2025 blooms were ~40 km southwest that week. Erie is **dropped as a case study
  and promoted to the tasking argument**: the only two open scenes over the western basin
  are a city and a river mouth, both >= 12 km from the entire GLERL monitoring network,
  which is a sharper argument for tasking than a weak case study would have been.

**Lake Ontario / Rochester replaces it** and is the strongest scene in the set: ~55% water
(48.4% of imaged pixels measured), a sharp cross-shore sediment front, 0% cloud, PACE at
17:39, 2,879 valid OLCI pixels, and it falls inside the ROCX 2025 field campaign window.
It also carries the PACE analysis — SF Bay yields only **9** usable PACE pixels against
Lake Ontario's 189, exactly the land-adjacency problem §11 predicted for a narrow bay.

### 13.2 Rules read; three consequences (adds to §0, §2, §8; supersedes README §1 "unverified")

The T&C (`docs/Planet-TermsConditions-TanagerCompetition.pdf`) allow third-party data,
code and APIs; require Tanager Open STAC imagery; accept 1–3 page memos, annotated
notebooks/repos, slides, ≤ 3 min video; one submission per participant; English; deadline
2026-08-31 23:59 PST; ≥ 3 Planet judges, mean score, top 10 to committee. Consequences:
1. **Memo gets the rubric's own headings** — Project Summary, Impact Statement,
   Limitations — and the notebook is secondary.
2. **`tanager_io` is promoted from `_io.py` to a deliverable** (§0 item 3): `walk_stac()`
   over `tanager_open_scenes` with `pystac` 1.14 (installed), `open_sr()` returning an
   xarray with wavelengths/FWHM as coordinates and the UTM geotransform parsed from
   `StructMetadata.0`, `water_mask()`. Small, tested, in the public repo.
3. **EMIT is not a match-up but can be a performance comparison.** `earthaccess` search
   2026-08-25: no EMIT L2A within ±10 d of any scene; nearest Golden Gate 2025-05-31
   (+17 d, 19:59 UTC), Erie 2025-08-18 (−27 d) and 2025-10-04 (+20 d); Loreto none
   within ±30 d. Optional day-4 item: Tanager vs EMIT over the same *stable* water
   (Golden Gate outer coast) — spectral sampling, reported SR uncertainty, 760/940 nm
   residuals — framed as sensor characteristics, never as the same water. Dropped first
   if time is short; PACE/OLCI already satisfy the +5 wording.

### 13.1 OLCI comes from Copernicus Marine L3, not Copernicus Data Space L2 (supersedes §2, §3, §11 row 1)

Charlie has a `~/.copernicusmarine` login; drift already has a `copernicusmarine`
adapter (`cmems_l4`). The CMEMS product `OCEANCOLOUR_GLO_BGC_L3_MY_009_103` carries
OLCI at **300 m, daily, global, 2016-04-25 → 2026-08-16** (MY; NRT only from 2026-08-03,
so every 2025 date is MY): `…plankton_my_l3-olci-300m_P1D` (CHL, CHL_uncertainty, flags;
ARCO subset works — 3,303 / 4,002 / 3,264 valid pixels for Golden Gate / Loreto / Erie)
and `…reflectance_my_l3-olci-300m_P1D` (RRS400, 412, 443, 490, 510, 560, **620**, 665,
674, 681, 709 + uncertainties; original-files only, one 5.75 GB global netCDF per day,
chunked 256 × 256 so a bbox range-read costs seconds).

Gains: no new account; atmospheric correction is the CMEMS/POLYMER chain (documented,
consistent with the PACE-vs-OLCI comparisons in exps 23/26); the grid is already 1/336°
lat/lon so the coarse-pixel polygons of §1 item 5 are trivial for OLCI; **Lake Erie is
covered** (the GLO product does not mask the Great Lakes). Losses: 11 bands instead of
21 (no 400–412 pair beyond RRS400/412, no NIR beyond 709 — the OLCI SRF resampling in
`_spectral.py` targets these 11); the overpass time is not in the file; S3A and S3B are
merged where both saw the scene. The CDSE swath route stays as a fallback in `scenes.yaml`
comments only.

