# Methods, validation and full results

Companion to **[MEMO.md](MEMO.md)**. Nothing here is a summary of that document — this is
the material the memo points to: how the measurement is made, the gate that licenses it,
the per-scene numbers behind every figure quoted there, and the limitations in full.

Every number below is reproducible from this repository; see the README. The gate ledger
is `data/metrics/verdicts.csv` and the statistics are `data/metrics/h_results.json` and
`data/metrics/variance_split.json`.

---

## 1. How the measurement is made

**Water.** Tanager ships no water mask. We build one from the Normalised Difference Water
Index with an Otsu-chosen threshold, a binary morphological opening to remove boats, wakes
and single-pixel shoreline speckle, and the product's own cloud and cirrus flags
(`tanager_io.water_mask`). Fill is −9999 on the float fields but 255 on the `uint8` flags,
which is a trap worth naming.

**The 30 m scalar.** Planet ships a dimensionless surface reflectance ρ tuned for land. We
subtract a per-band additive pedestal — the 1st percentile over the scene's water pixels,
taken *per band* so the correction is a spectrum rather than the single number a
near-infrared dark-pixel correction assumes — and divide by π, since ρ = π·Rrs
(`tanager_io.to_rrs`). The primary scalar is then the **green/blue ratio, 560/490 nm**
(`tanager_io.band_ratio`), computed per 30 m pixel *before* any aggregation: the mean of a
ratio is not the ratio of the means, and that difference is what this experiment is about.
A ratio also divides out multiplicative error, which an uncorrected atmosphere largely is.

A second scalar, plain 560 nm reflectance, is carried through the identical path and
reported beside the first throughout — see §3.

**Coincidence.** The memo's table gives the same PACE time, 19:35:44 UTC, at San Francisco
Bay and at Loreto, 109 days apart. That is not a transcription slip: both L2 granules are
stamped `T193544` (`pipeline/scenes.yaml`, and the NASA CMR catalog lists both).

**Binning.** Each 30 m water pixel is assigned to the nearest coarse pixel center — a
Voronoi partition, via a KD-tree in the scene's UTM meters. For a smooth swath this is
equivalent to rasterising footprint polygons and has no topology to get wrong. A coarse
pixel enters the analysis only if **≥ 90 %** of it is Tanager water, which is what keeps
shoreline and swath-edge cells out. That leaves a median of 308–378 Tanager pixels inside
each OLCI cell and 1,590–3,075 inside each PACE pixel.

**The protocol, applied as the protocol applies it.** Bailey & Werdell (2006) take a 5×5
box of coarse pixels, require ≥ 50 % of it valid, apply an iterated 1.5σ rejection, and
accept the match-up when the coefficient of variation (CV) across the surviving pixels is
below 0.15. We reproduce that on **the coarse sensor's own retrieval** — OLCI's 560 nm
reflectance, or PACE's `chlor_a` — because that is what a real match-up filters on. That
statistic decides which boxes pass, and nothing else.

**The variance split** (the memo's headline) uses none of that. It applies the law of
total variance to a single 480 × 480 pixel window — 14.4 km, wider than a 5×5 PACE box —
block-averaged into square blocks of 18 × 18 Tanager pixels (540 m, standing in for the
OLCI cell) and 40 × 40 (1,200 m, for PACE). The blocks are not the coarse sensors' own
grid cells; they are a clean partition of the window at each sensor's scale, which is
what the identity needs. The window is trimmed to a whole number of blocks (468 px at the
OLCI scale, 480 at PACE's) and both the total and the between-block variance are
computed over the same trimmed pixels, each block weighted by its count of finite 30 m
pixels. The between-cell and within-cell parts then sum to the total exactly, which is
why this is the honest ladder statistic; a ratio of panel variances is not, and can
exceed one. (An earlier draft took the total over the untrimmed window; that leaked up
to one percentage point, and the numbers here are the exact ones.)

**A guard worth declaring.** A coefficient of variation diverges as its mean approaches
zero, so a coarse pixel whose mean green/blue ratio falls outside 0.2–5.0 is not reporting
a water spectrum but a failed pedestal removal. Those are excluded wherever the
within-pixel CV is used, identically in the headline table and in the structure score.
Without it a handful of near-zero denominators dominate every upper percentile. The
brightness variant is guarded only on a positive mean — a near-zero mean *inflates* its own
within-pixel CV, so a permissive guard can only push that variant toward the ratio's
answer, never away from it, and the gap reported in §3 is therefore a floor.

## 2. Validation: does Tanager's variance mean anything?

Before any 30 m claim, Tanager has to reproduce a sensor that has none of its problems at a
scale that sensor can see. Aggregating Tanager onto the OLCI grid and comparing, per band:

| scene | 443 nm | 490 nm | 560 nm |
|---|---|---|---|
| Lake Ontario (n = 952) | r 0.78 · ×0.90 | r 0.92 · **×1.13** | r 0.95 · ×1.09 |
| San Francisco Bay (n = 205) | r 0.67 · ×0.95 | r 0.82 · ×1.53 | r 0.88 · ×1.80 |
| Loreto (n = 685) | r 0.58 · ×1.03 | r 0.81 · ×1.49 | r 0.82 · ×1.44 |

Two numbers per band, because **a correlation cannot license a claim about variance**: it
is scale-invariant, and would read 1.0 on a Tanager field with ten times OLCI's spread. The
second number is the regression slope — the factor by which Tanager over- or under-states
OLCI's spatial variability at ~550 m. Both sides are remote-sensing reflectance in sr⁻¹,
so parity is a slope of 1.

Tanager reproduces OLCI's variability to within **13 %** at Lake Ontario (×1.13 at 490 nm,
×1.09 at 560 nm) and runs 44–80 % high at Loreto and San Francisco Bay, the two scenes with a heavier atmospheric residual
(§5). It is not a calibration; it is enough to license a claim about a factor of three to
five.

**443 nm is the weakest band everywhere and is reported rather than hidden.** Nothing
downstream is computed from it — the 30 m scalar is the 560/490 ratio, so 490 and 560 are
the bands the gate binds on, and the gate takes the worst of those two across the worst
scene rather than the best of nine.

## 3. Full results in the protocol's terms

Across 1,528 boxes that pass the standard filter, the variance inside the accepted pixels
against the variance across them — both measured on the same Tanager green/blue ratio over
the same 5×5 window. The box column below sums to 1,526: two Lake Ontario / OLCI boxes pass
the filter (675 do) but fall outside the 0.2–5.0 ratio guard of §1, so they count toward
the 1,528 and toward the last two columns, and not toward the CVs. San Francisco Bay / PACE
has no row: only 9 PACE pixels in that frame are ≥ 90 % water, so no 5×5 box can meet the
protocol's 50 % validity requirement, and nothing there passes or fails.

| scene / sensor | boxes passing | between-pixel CV | within-pixel CV | ratio | passing boxes whose *center* pixel has interior CV > 0.15 | … containing *any* such pixel |
|---|---|---|---|---|---|---|
| Lake Ontario / OLCI | 673 | 0.026 | **0.139** | **5.4×** | 48 % | 66 % |
| Lake Ontario / PACE | 100 | 0.023 | **0.106** | **4.6×** | 30 % | 69 % |
| Loreto / OLCI | 470 | 0.037 | **0.135** | **3.6×** | 38 % | 71 % |
| Loreto / PACE | 40 | 0.039 | **0.129** | **3.3×** | 35 % | 98 % |
| San Francisco Bay / OLCI | 243 | 0.017 | 0.022 | 1.3× | 1 % | 29 % |

Weighted by boxes, the two right-hand columns are the memo's 36 % and 63 %.

Both columns are the same scalar on the same instrument, differing only in scale, so the
ratio is a decomposition rather than a comparison of unlike things. The protocol's own
filter statistic — OLCI's 560 nm reflectance, or PACE's chlorophyll — is a different
quantity on a different instrument. It decides which boxes pass, it is recorded alongside
these numbers in `h_results.json`, and it is **not** the denominator of any ratio here.

**The brightness variant.** Repeating the table on 560 nm reflectance gives 0.6× (Lake
Ontario / OLCI), 0.4× (Lake Ontario / PACE), 1.9× and 1.2× (Loreto), and 0.9× (San
Francisco Bay) — the brightness row of the memo's split, restated in the protocol's terms.
Both variants come from the same code and both are in the verdict ledger; neither is shown
without the other. A reader who repeats this analysis on a single band will get the second
set of numbers, and that is the brightness/color result, not a discrepancy.

**San Francisco Bay is worth reading twice.** Its color ratio of 1.3× is the lowest here,
but that is a statement about *how much* structure the scene has, not where it lives: the
water is uniformly turbid at half a kilometer, so the within-pixel CV is 0.022 against
0.13–0.14 elsewhere. On the scale question it agrees with the others — 70 % of its color
variance is inside an OLCI cell.

## 4. The three supporting checks

- **It is water, not noise.** Within-pixel variability is spatially structured: adjacent
  coarse pixels differ by 0.17–0.44 of what the same field gives after being spatially
  shuffled (a shuffled field scores 1.0 by construction). Structured in 4 of 5 assessable
  scene/sensor pairings; above the 0.75 bar at 0.79 for Loreto/PACE — a fail, which we
  report rather than drop — and not assessable on San Francisco Bay/PACE, where only 19 adjacent finite pairs
  exist and the statistic is undefined rather than failed.
- **It clears the instrument's own floor.** Median within-pixel CV is 1.5× (San Francisco
  Bay), 5.1× (Lake Ontario) and 5.5× (Loreto) the relative uncertainty Planet's own
  `surface_reflectance_uncertainty` declares for those pixels.
- **The spectra agree band by band** — §2.

## 5. Limitations, in full

- **Tanager's atmospheric correction is tuned for land, and over water its quality is
  scene-dependent.** Median reflectance at 1600 nm, where physical water is under 0.01:
  Lake Ontario 0.0025 (physical), Loreto 0.0090 (close), San Francisco Bay 0.0971 — the
  scene Planet flags at 29 % light haze. This tracks the validation table in §2: Lake
  Ontario, the cleanest correction, is also the scene that reproduces OLCI's spread to
  within 13 %. The diagnostic is free and we recommend it as routine practice for anyone
  using Tanager surface reflectance over water. We removed the residual per band from the
  darkest water in each scene.
- **Every CV here is a lower bound.** An additive residual inflates a mean and leaves a
  standard deviation alone, so it deflates a coefficient of variation. San Francisco Bay
  is the most affected, which is a genuine confound on reading it as a clean control.
- **The hidden-variance fractions are upper bounds on hidden *water* variance**, because
  the pushbroom sensor's column striping counts toward the within-cell term. The
  cross-scene consistency (68–71 % / 75–78 %) is evidence against striping dominating —
  the three scenes are different strips on different dates — but it does not eliminate it.
- The CMEMS `olci-300m` product is **mapped to 1/180°, i.e. 488 × 618 m at 38 °N**, not
  300 m. The ladder is 30 m → ~550 m → ~1.2 km, and the memo says ~550 m for that reason.
- Variability is measured on a **band ratio**, not chlorophyll: a ratio survives
  multiplicative error, and chlorophyll algorithms are not valid in the turbid and
  optically shallow parts of these scenes. The ratio is the color signal a chlorophyll
  retrieval is built from, not the retrieval itself — and the step from one to the other
  is not measured here.
- **No in-situ data are used.** That is the point of the tasking section: this study
  demonstrates that the sub-pixel term is measurable and large, not that any particular
  retrieval is wrong.

## 6. Gates, and which were pre-registered

Seven gates were written down before the run and each is recorded in
`data/metrics/verdicts.csv` with its statistic, threshold, verdict and a note.

| gate | asks | verdict |
|---|---|---|
| G0 | every scene has a cube and coincident match-ups that open | PASS |
| V0 | surface reflectance is usable over water | PASS |
| V1 | the grid is where `StructMetadata.0` says it is | PASS |
| V2 | Tanager reproduces OLCI's variability at the scale OLCI can see | PASS |
| V3 | within-pixel variability is spatially structured, not noise | **MIXED** |
| V4 | it exceeds the declared noise floor | PASS |
| V5 | the filter accepts boxes with variable interiors | ANSWERED |
| V6 | the tasking proposal names boxes with in-situ truth and no scene | PASS |

V3 is recorded as MIXED, not PASS: 4 of 5 assessable pairings are structured and
Loreto/PACE sits above the bar at 0.79. V5 was pre-registered as "reported either way" —
it is the headline, and pre-committing to report it in whichever direction it came out is
the reason it can be believed.

**One gate was added on revision and is labelled as such in the ledger.** V2 originally
tested a correlation alone. A correlation is scale-invariant and therefore cannot license
a claim about variance, which is what everything downstream of V2 is. The spread half —
the regression slope in §2 — was added afterwards, and the ledger says so rather than
presenting it as foresight.

## 7. Reproducing this

Everything is computed by six numbered stages from the open scenes. Stages 3 and 5 re-run
from the data committed in this repository with no credentials and no downloads; see the
README for the exact commands and for what the remaining stages need.

The reader is released as a package:

```bash
pip install "git+https://github.com/cjt31415/tanager-io.git"
```

```python
import tanager_io as tio
ds    = tio.open_sr(tio.asset_url("20250919_170233_04_4001"))  # lazy, labelled, georeferenced
water = tio.water_mask(ds)                                     # NDWI + cloud flags
rrs   = tio.to_rrs(ds, water=water)                            # pedestal removed per band, / pi
ratio = tio.band_ratio(ds, 560, 490, source=rrs)               # the 30 m scalar measured here
```

Those are the calls the analysis itself makes — every stage goes through the package, so
the snippet above and the numbers in these documents cannot drift apart.

`tanager-io` exists because a naive `xr.open_dataset` on a Tanager scene returns an empty
`Dataset`. It absorbs the four things that cost us a day: the data group's name contains a
space; the wavelengths live in an attribute rather than a coordinate; the georeferencing
hides in a `StructMetadata.0` text blob and the grid is corner-registered (half a pixel,
silently, if missed); and fill is −9999 on the floats but 255 on the `uint8` flags. MIT
licensed, 24 tests, no network required to run them.
