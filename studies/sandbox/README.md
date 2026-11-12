# Sandbox: fast option tests

A small case (`prod_case("sandbox")`: Midwest part of the 21–22 October 2026 night, 28 radars, 6 frames 00:30–03:00 UTC, 4 km × 100 m, built by P1–P3 with no near-range or no-echo cuts) for trying inversion options quickly (10–250 s per run) and deciding later.

```matlab
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
S = sandbox_load("sandbox");                       % once, ~6 s
r = sandbox_run(S, struct("name", "my test", "prior", 0.5));
summary = sandbox_summary({r}, "my_test");          % table + comparison figure (list several runs)
run(fullfile(project_root, "studies", "sandbox", "SB_fast_solver.m"))   % an example test script
```

An option is a struct with a `name` and any of: `prior`, `smoothness`, `advection`, `coeff_horver`, `lambda_profile`, `prior_interpolation`, `prior_length_scale_m`, `prior_mean`, `shape`, `prior_correlation_m`, `row_weight` (a function `w = f(range_m, elev_deg, is_no_echo)` weighting each observation row), `nonnegative` (default false: plain PCG) and `video` (default true).

**Outputs of every run**
- `figures/runs/<name>.mp4`: one frame per time with the column map, the column / regional mean map (the radar imprint: blue discs or red rings around the radars), the KDMX section and the frame's radar-distance curve. `<name>.png` is the middle frame.
- `figures/sandbox_log.csv`: one row per run (time, name, scores, solve time), appended, so progress can be followed across sessions.
- `sandbox_summary(results, test_name)` writes `figures/<test_name>.csv` and `.png`: mean radar-distance curves, bias per frame and the bias / trusted-misfit trade-off.

**Scores**
- *Radar bias:* column density / regional mean (150 km Gaussian) against distance to the nearest radar, RMS of the log ratio within 200 km (0 = no imprint).
- *Misfit:* relative data misfit on all rows, and on trusted rows only (30–120 km from the radar), since the near-range rows are biased low.
- *Detail:* RMS of log(column / regional mean), the small-scale structure kept (lower = smoother map).
- *Negative share:* negative mass / total absolute mass (0 with `nonnegative`).
- *Success:* bias near 0 (a flat radar-distance curve at 1, no discs or rings in the relative map) without a higher trusted misfit or a much lower detail. Radar cross-validation (`sandbox_fold.m`) is the out-of-sample check.

**Solver options added for this** (`solve_advection_lean`, off by default):
- `row_weights`: soft data weights (rows and values scaled by √w).
- `shape`: separable penalty λ²‖x − p Σ_z x‖² per column, p the prior's vertical shape (column total free).
- `prior_correlation_m`: prior term (x − x̄)ᵀ(I + (ℓ/Δx)² L_h)(x − x̄), a first-order Matérn/SPDE Gaussian-Markov prior on the deviation from the prior, horizontally correlated over ℓ.

## Results (23 October 2026)

| Option | Bias | Misfit | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline (no cuts) | 0.132 | 0.243 | 0.226 | 0.37 | 25 s |
| hard cuts (drop < 25 km; no-echo rows only < 75 km) | 0.049 | 0.360 | 0.223 | 0.34 | 25 s |
| Gaussian tapers 15 km / 25 km | 0.122 / 0.106 | 0.248 / 0.259 | 0.226 / 0.227 | 0.36 / 0.35 | 26 s |
| logistic 20–75 km | 0.091 | 0.276 | 0.226 | 0.35 | 27 s |
| **logistic 28–75 km** | **0.056** | 0.324 | **0.226** | 0.34 | 50 s |
| logistic 28–75 km + GP 50 km | 0.028 | 0.419 | 0.307 | 0.30 | 101 s |
| shape 1 / shape 3 | 0.139 / 0.133 | 0.294 / 0.403 | — | — | 25–33 s |
| GP 50 / 100 / 150 km (prior 0.3) | 0.121 / 0.077 / 0.052 | 0.344 / 0.413 / 0.453 | — / 0.357 / — | — / 0.27 / — | 80–200 s |
| GP 100 km, prior 0.1 | 0.166 | 0.306 | 0.275 | 0.38 | 91 s |
| soft 25 + GP 100 km + shape 1 | 0.025 | 0.466 | 0.386 | 0.30 | 243 s |

Logistic weights: near range 1/(1 + exp(−(r − r0)/3 km)); no-echo rows 1/(1 + exp((r − 75 km)/10 km)). Gallery: `figures/gallery_key_options.png`.

**Reading**
- The imprint is mostly the suppressed near-range data. Down-weighting them removes most of it without changing the fit to the trusted rows: a steep but smooth (logistic) transition centred near 28 km is as good as the hard 25 km cut.
- The correlated (GP) prior flattens the remaining imprint but also spreads errors: alone, with a weak prior weight, it widens the near-radar deficit (bias 0.166); with good data weights it reaches 0.03 but costs trusted misfit (+35%) and detail (maps visibly smoother).
- The shape (separable) penalty as implemented does not reduce the imprint; it holds each column to the prior's vertical shape, which is the same shape near and far from radars, so it does not address a horizontal deficit.

## Logistic weights + prior strength and GP length (23 October 2026)

All runs with logistic weights (near range centred at 28 km, no-echo rows fading around 75 km). `figures/gp_tradeoff.png`, `figures/gp_tradeoff.csv`, maps in `figures/gallery_gp_candidates.png`.

| Prior | GP length | Bias | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.3 | — (white) | 0.056 | 0.226 | 0.34 | 23 s |
| 0.1 | 20 / 30 km | 0.090 / 0.085 | 0.228 / 0.235 | 0.40 / 0.39 | 50 s |
| 0.2 | 20 / 30 / 50 km | 0.062 / 0.056 / 0.044 | 0.245 / 0.257 / 0.279 | 0.35 / 0.35 / 0.33 | 38–63 s |
| 0.5 | — (white) | 0.042 | 0.245 | 0.30 | 17 s |
| 0.3 | 20 / 30 / 50 km | 0.047 / 0.040 / 0.028 | 0.262 / 0.279 / 0.307 | 0.32 / 0.32 / 0.30 | 36–67 s |
| 1.0 | — (white) | 0.026 | 0.321 | 0.22 | 14 s |

**Reading**
- What flattens the imprint is the overall stiffness of the prior: a weaker prior (0.1) lets the data imprint come back (bias 0.085–0.090), whatever the GP length.
- There is one trade-off curve between radar bias and fit to the trusted data. For moderate flattening (bias ≈ 0.04) the stronger white prior (0.5) is the cheapest: best trusted misfit (0.245), fastest, similar detail. For strong flattening (bias ≈ 0.027) the GP wins clearly: GP 50 km at prior 0.3 keeps detail 0.30 at misfit 0.307, while a white prior of 1 drops detail to 0.22 and dims the whole map.
- Visually: white 0.5 keeps the fine structure with faint marks left at a few radars; GP 30–50 km removes the marks but blurs the migration features; white 1 dims everything.
- Candidates: *logistic 28 + white prior 0.5* (sharp, cheap, faint residual marks) or *logistic 28 + GP 30 km, prior 0.3* (no visible marks, smoother).

## Plain PCG without non-negativity (24 October 2026)

`SB_fast_solver.m`, `figures/fast_solver.png`, `.csv`. The projected FISTA step after PCG makes the solver about 1.6× slower and changes nothing measurable: the scores agree to 3 decimals, and PCG alone leaves 0.2% negative mass (0.03% with prior 0.5). **Plain PCG is now the sandbox default** (`nonnegative=false`); non-negativity can be added back for production.

| Option | Bias PCG / FISTA | Trusted misfit PCG / FISTA | Detail | Negative mass (PCG) | Time PCG / FISTA |
| --- | ---: | ---: | ---: | ---: | ---: |
| no cuts | 0.132 / 0.132 | 0.226 / 0.226 | 0.37 | 0.17% | 15 / 24 s |
| hard cuts (25 km, no-echo < 75 km) | 0.049 / 0.049 | 0.223 / 0.223 | 0.34 | 0.18% | 15 / 23 s |
| logistic 28, prior 0.5 | 0.042 / 0.042 | 0.245 / 0.245 | 0.30 | 0.03% | 11 / 18 s |

Times are the solve only (PCG to 1e-4, 19–31 iterations); the density prior is cached in `S` after the first run.

## Near-radar artefact at fine range, per radar (24 October 2026)

`SB_near_range.m`: fine distance bins (1 km to 20 km, then 2.5, 5, 10 km; observations in 2 km gates from 4 km) and one line per radar, instead of the pooled 25 km bins of `radar_distance_bias`. Figures `near_range_columns.png`, `near_range_observations.png`, `near_range_vertical.png`.

- **Columns** (column / 150 km regional mean, own radar). Without cuts the median drops to 0.4–0.5 at 4–12 km and recovers only at 35 km. With hard cuts it is 0.9–0.98 inside 25 km, with a dip to 0.81 at 25–28 km, just outside the cut. Both runs fall to 0.84–0.88 at 100–150 km. The radars differ widely (0.3–1.1 at 5–15 km without cuts).
- **Held-out observations** (observed / predicted from the hard-cut solution, which never saw the rows within 25 km). The median deficit depends on elevation. At 0.5° it reaches 1 at 18 km, at 0.9–1.3° at 25 km, at 1.8–2.4° at 28–30 km, and at 3.1° and above it stays at 0.8–0.9 to 45 km. The 25 km cut is therefore not enough above 1.5°. Per radar, the spread at 0.5° and 10–20 km is 0.14–2.5 (10th–90th percentile): some radars are strongly suppressed, while others read up to 3× above the prediction.
- **The deficit is in the values, not in no-echo zeros.** No-echo gates are under 1% of near rows at 0.5–1.3° and 2–9% above. On the non-zero rows, the pooled ratio at 4 km is 0.17 at 0.5° and 0.03 at 3.1°.
- **Vertical.** Both solutions put a dense low layer (0–800 m AGL, about 45 birds km⁻³) at 15–60 km, where the 0.5–1.3° beams see it. Beyond about 60 km the beams pass above it, and the low layers come from the flatter prior profile (about 25). This matches the 0.85 column deficit at 100–150 km, so this part of the imprint is vertical sampling, not near-range suppression. Without cuts, the near-radar column is depleted at all heights inside 25 km.

## Gate estimator, no cuts, and a minimum range (1 November 2026)

`sandbox` now uses the production cleaning with the gate estimator (`studies/gate_density`) and no near-range or no-echo cuts; `sandbox_median` is the earlier sandbox build (median cleaning, no cuts). Minimum ranges as hard row weights. `figures/estimator_min_range.png`, `.csv`.

| Option | Bias | Misfit | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| median cleaning, no cut | 0.132 | 0.242 | 0.226 | 0.37 | 17 s |
| estimator, min 3 km (`minRange`) | 0.065 | 0.195 | 0.184 | 0.29 | 15 s |
| estimator, min 5 / 10 / 15 km | 0.065 / 0.066 / 0.068 | 0.201 / 0.221 / 0.262 | 0.184 | 0.29 | 15–16 s |

**Reading**
- The estimator removes the dark disc: within 25 km the column density is 1.02–1.06 of the regional mean, against 0.73 with the median cleaning, and the bias halves without any cut (0.065, against 0.049–0.056 for the old 25 km / 75 km cuts on the median cleaning).
- A minimum range beyond 3 km changes nothing (bias 0.065–0.068, same trusted misfit): no near-range cut is needed.
- What remains is at far range: a bright ring at 25–50 km (1.08) and a dip at 100–150 km (0.92), then 1.05 beyond 175 km. Misfits are against each cleaning's own data, so they compare options within a cleaning, not across.

## The far-range dip (2 November 2026)

`SB_far_range.m`. Same sandbox, estimator, no cuts.

| Distance to radar (km) | 12.5 | 37.5 | 62.5 | 87.5 | 112.5 | 137.5 | 175 | 250 | Bias |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| base | 1.02 | 1.08 | 1.01 | 0.95 | 0.92 | 0.94 | 1.03 | 1.05 | 0.065 |
| flat prior (one profile per time) | 1.03 | 1.09 | 1.02 | 0.96 | 0.94 | 0.94 | 1.00 | 1.02 | 0.062 |
| rows to 100 km only | 0.99 | 1.05 | 0.98 | 0.93 | 0.99 | 1.01 | 1.02 | 1.02 | 0.054 |

**Reading**
- The dip at 100–150 km comes from the observations beyond 100 km, not from the prior: without those rows it goes; with a flat prior it stays.
- Those rows observe less than the field from the nearer radars predicts (held out of the solve): 0.80 overall, 0.70–0.83 for rows that see 1.3–2.7 km, 0.96 for rows that see ~1.1 km.
- The deficit fades over the night: 0.43, 0.61, 0.71, 0.81, 0.91, 1.01 from 00:30 to 03:00 UTC. Near the radars the profile is steep early (30 → 8 birds km⁻³ from 950 to 1950 m at 00:30) and flatter later (49 → 32 at 03:00).
- Reading: far beams sit in the steep upper part of the layer, so a small error in their modelled height or width gives a large density error. A beam a few hundred metres higher than the 4/3-Earth model (sub-refraction) or narrower than modelled would read low early and fine later. Not tested: it needs kernels with another refraction or beam width.
- The bright ring at 25–50 km mostly follows from the dip (each curve is normalized to mean 1); 1.05 remains without the far rows.
- With the estimator, gates without echo get a small positive density, so no row is all zero: the `is_zero` argument of `row_weight` no longer selects anything.

**1 km bins** (`SB_far_range.m` section 4, `figures/far_range_fine.png`, `.csv`; grid columns are 4 km, so nothing finer shows)
- Near the radar the curve dips to 0.79–0.90 at 2–8 km, is back at 1 by 12 km, then 1.05–1.10 from 20 to 50 km. The 25 km bins hid the dip.
- The imprint map is mostly per-radar discs: each radar's coverage sits uniformly above or below its neighbours, with circular edges where the nearest radar changes. Differences in level between radars, not a distance effect, dominate the map.
- Vertical imprint: the density is high along the lowest beam (0–500 m near the radar, ~1.5 km at 80 km) and low above 2 km within 50 km; beyond 150 km, above 2 km, it is high.
- Beam lines on the vertical imprint (4/3 Earth, mean radar height): the high band near the radar sits inside the 0.5° beam (between its half-power edges) out to ~70 km; the high far-range region above 2 km also lies along the 0.5° beam; the air below the 0.5° beam beyond ~40 km is low. The solution puts birds where the lowest sweep looks.
- First 25 km in 0.5 km bins, mean over radars with a 95% bootstrap interval over radars: 0.8–0.95 at 2–8 km, with the interval below 1 only at 2.25, 4.25–4.75 and 6.25 km; back at 1 from 9 km.

**What to adjust for the gaps at 100–150 km** (quick runs, same sandbox; columns: curve at 87.5 / 112.5 / 137.5 km)

| Option | Bias | Trusted misfit | Detail | 87.5 | 112.5 | 137.5 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| base | 0.065 | 0.184 | 0.29 | 0.95 | 0.92 | 0.94 |
| row weight tapering from 1 at 75 km to 0.5 / 0.2 / 0.05 at 150 km | 0.062 / 0.060 / 0.058 | 0.184–0.185 | 0.28–0.29 | 0.95 | 0.92–0.93 | 0.94–0.97 |
| rows to 100 km only | 0.054 | 0.206 | 0.28 | 0.93 | 0.99 | 1.01 |
| shape 1 | 0.111 | 0.235 | 0.39 | 0.92 | 0.86 | 0.87 |
| prior 1 / 3 | 0.039 / 0.017 | 0.272 / 0.442 | 0.19 / 0.11 | 0.96 / 0.98 | 0.95 / 0.98 | 0.97 / 1.00 |

- Down-weighting far rows helps little; a stronger prior flattens the curve only by giving up the data (misfit ×1.5–2.4, detail ÷1.5–2.6); the regional shape makes it worse.
- Per-radar level within 100 km (column / regional mean, frame mean, geometric mean over the disc) ranges from 0.64 (KVNX) to 1.34 (KMPX), sd ×1.24. Radars with a high disc drop most at their edge (KHPX 1.18 → 0.81, KIWX 1.14 → 0.81, KOHX 1.22 → 0.85). The steps between radars, not the distance, make the gaps.

## Dry nights for a multi-night comparison (2 November 2026)

`SB_night_screen.m`, `figures/night_screen.png`, `.csv`. Dark Ecology daily series 2021–2025, night periods, the 28 sandbox radars. Dark Ecology rain shares are rarely below 1% (median 6% over radar-nights), so a night qualifies with rain below 5% at 24 or more radars, at least 26 radars reporting (<10% missing), and a median night-mean bird reflectivity of at least 300 cm² km⁻². 11 nights qualify (3 spring, 8 autumn); the 3 spring and the 6 strongest autumn nights are kept: 2021-05-12, 2023-05-17, 2023-05-21, 2025-10-07, 2024-10-06, 2021-09-08, 2025-09-25, 2025-09-27, 2021-10-21.

First look at radar levels, from the same series (all dry radar-nights, rain < 5%, reflectivity > 100 cm² km⁻²): each radar's night-mean reflectivity relative to the median of its neighbours within 300 km.
- The long-run level ranges from 0.37 (KFSD) and 0.44 (KDDC) to 1.65 (KOAX) and 1.69 (KMPX); this mixes calibration and real geography (the low ones are on the western edge).
- From night to night the level varies by a factor of 1.3–1.9 (sd), so a single night's disc level of 0.64–1.34 is within normal variation.
- It correlates only weakly with the sandbox disc levels (r = 0.22 on log scale, Spearman 0.32, 27 radars). A fixed per-radar offset explains little of the sandbox discs; an overlap test on the volumes is needed to separate calibration from birds.

**Lowest-cut mosaics** (`SB_night_maps.m`, `figures/night_maps.png`): one Level II volume per radar and night, closest to 3 h after sunset at the domain centre (251 of 252 found, 2.3 GB in `local_data_dir/sandbox_nights/raw`; KGRR 2021-09-08 missing, KVNX 2024-10-06 unreadable). Raw reflectivity, each pixel from its nearest radar within 150 km.
- Every radar shows a bright disc fading with range, as the 0.5° beam rises above the bird layer.
- Some radars sit clearly below their neighbours on a given night (e.g. KILX on 2021-10-21, KUEX on 2021-05-12 and 2023-05-17), with sharp edges where the nearest radar changes.
- Rain or other echoes remain on some nights (bands on 2021-09-08 and in the south-east on 2025-09-25), so these two are not fully dry.

## Overlap test: neighbouring radars on the same air (2 November 2026)

`SB_overlap.m`, `figures/overlap.png`, `overlap_pairs.csv`, `overlap_offsets.csv`. Volumes of `SB_night_maps.m` (9 nights, same UTC time, volumes within 5 min of it). Each pair of radars less than 300 km apart is compared on gates where both beams are at the same height (within 100 m) and range (within 10 km), i.e. near the midline: same air, same height, same beam size. Rain (REF ≥ 20 dBZ and RHO ≥ 0.97), echoes above 35 dBZ and missing gates are dropped. Median difference in dB per pair and night (median 2440 gates, at ~1950 m on the 0.5° cut, ~2600 m on the 0.9° cut); per night, one offset per radar (mean zero) by least squares.

- Within a night, pairs differ by 2.1–3.2 dB (sd), and one offset per radar explains most of it (residual 0.5–1.5 dB): the differences are radar-specific.
- Across nights they change: the mean offsets spread 1.6 dB (KFSD −5.5, KABR −2.3, KSGF −1.7 … KLSX +2.1, KIWX +2.5, KPAH +2.5), the night-to-night sd is 2.45 dB; 30% of the variance repeats. The mean offsets follow the Dark Ecology long-run levels (r = 0.61).
- On the 0.9° cut the offsets are a third as large (slope 0.36, r = 0.48, 199 radar-nights). A calibration error would be the same on both cuts, so most of the 0.5° offset depends on the beam height: propagation (or pointing), per radar and per night. A per-radar scale factor would fix the repeating part only.

## Does beam propagation explain the overlap offsets? (2 November 2026)

`SB_propagation.m`, `figures/propagation.png`, `propagation_pairs.csv`. Refractivity from radiosondes (IEM archive, 18 sites, 00 and 12 UTC; 303 of 324 usable), inverse-distance to each radar and linear in time to the volume time; beam centre by ray tracing (checked: −3 m from 4/3 Earth at 150 km in a standard atmosphere); bird profile from each radar's own volume at 5–35 km; predicted pair difference = beam-weighted profile at A's true height minus at B's. No fitting.

| Cut | Pairs | Measured sd | Predicted sd | r | Variance explained |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.5° | 394 | 2.66 dB | 0.65 dB | 0.06 | none |
| 0.9° | 321 | 2.27 dB | 0.63 dB | 0.20 | 3% |

- With the soundings, the beams at the matched gates are within −74 to +40 m of the 4/3-Earth height (10–90%); the mean gradients over 0–3 km are −40 to −80 N km⁻¹, no ducting. That moves the pair difference by well under 1 dB, and not in the measured direction.
- So propagation, as far as radiosondes resolve it, does not explain the 2–3 dB radar offsets. Not excluded: local inversions or ducts between the sounding sites or hours (soundings are ~300 km apart, 2–4 h from the volumes).
- Sensitivity: at the matched gates (median range 111 km, beam sd 540 m) a 100 m shift changes the beam-weighted profile by 0.9 dB (10–90%: 0.5–1.5 dB), so shifts of tens of metres are resolved (heights are continuous, the 100 m profile is interpolated). The two radars of a pair differ by 42 m (sd; max 325 m) with the soundings; explaining the measured differences would need a median 170 m (10–90%: 0–590 m), i.e. a refractivity-gradient difference of ~30 N km⁻¹ over the whole 0–2 km layer between neighbouring radars.

**Censoring or calibration?** (`SB_overlap.m` with each sweep's detection threshold, c + 20 log10(range), c from the 1st percentile of reported values at 20–150 km)
- The matched gates are mostly close to the noise: median margin above threshold 4.3 dB (0.5°) and 3.3 dB (0.9°); only pairs less than 150 km apart (matched at ~77 km, ~1.2 km height) have a median margin of 11 dB. The thresholds differ by 5.5 dB (sd) between radars.
- The disagreement shrinks with the margin: mean |pair difference| on 0.5° is 2.3 dB below noise, 1.8–2.0 dB at 0–9 dB, 1.2–1.3 dB above 9 dB; on 0.9° from 2.0 to 0.7 dB. Pairs less than 150 km apart: median 0.0 dB, sd 1.75 dB.
- The pair difference follows the difference in detection threshold (r = 0.36): the radar with the higher noise reads higher, as expected when only gates detected by both are compared.
- Selecting strong gates by their own values (both ≥ 5 dB above threshold) still gives a 0.9°/0.5° slope of 0.27, but that selection is biased towards the noisier radar.
- Reading: most of the radar-to-radar disagreement at the overlaps comes from echo near the detection limit; with strong echo the radars agree within ~1 dB (calibration-scale). Beam propagation (soundings) adds < 0.5 dB.

## Where the column goes with distance (2 November 2026)

Base sandbox solve (estimator, no cuts), mean over frames: mean column contribution (birds km⁻²) per height band, by distance to the nearest radar.

| Distance (km) | 0–500 m | 500–1000 m | 1000–1500 m | 1500–2000 m | 2000–3000 m | Column |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0–25 | 16.8 | 21.8 | 12.5 | 9.8 | 11.9 | 72.8 |
| 25–50 | 12.4 | 24.3 | 17.3 | 11.7 | 12.3 | 78.0 |
| 50–100 | 9.4 | 17.5 | 15.0 | 12.6 | 15.1 | 69.5 |
| 100–150 | 7.8 | 16.2 | 14.4 | 11.7 | 16.2 | 66.3 |
| 150–250 | 7.0 | 16.2 | 15.3 | 13.0 | 20.0 | 71.5 |

- The column changes little with distance, but its profile does: below 1 km it falls from 38.6 to 23–24 birds km⁻² beyond 100 km, while above 2 km it rises from 11.9 to 16–20. Birds are low near the radars and high between them.
- Below 1 km between radars no beam looks: that air is filled by the solver from the beams above and by smoothness, not from the well-observed profile near the radars. Without the rows beyond 100 km only the 2–3 km band at 100–150 km changes (16.2 → 20.0).

The prior mean (`build_density_prior_fast`: one profile per radar and time fitted to all of the radar's rows, Gaussian-interpolated) has the same shape at every distance, and it is the shape seen between radars, not near them:

| Distance (km) | 0–500 m | 500–1000 m | 1000–1500 m | 1500–2000 m | 2000–3000 m | Column |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0–25 | 10.5 | 16.9 | 14.6 | 11.8 | 17.8 | 71.6 |
| 100–150 | 7.8 | 16.4 | 14.9 | 12.7 | 20.2 | 71.9 |

- Near the radars the data pull the low layer up (16.8 + 21.8 against the prior's 10.5 + 16.9) and the upper layer down; between radars the solution falls back to the prior (7.8 / 16.2 at 100–150 km, same as the prior).
- The profile fit is dominated by the many far rows, which see only the upper part of the layer and the weakest echo, so the prior profile is too low near the ground and too high aloft. The vertical imprint is the solve going from data-driven near the radars to this prior between them.

## Profile prior and correlated prior on the estimator data (2 November 2026)

**Profile prior** (`SB_profile_prior.m`, `figures/profile_prior.png`, `.csv`). The radar profiles of the prior mean were smoothed in height with `lambda_profile` = 10, which flattens the low-level peak: too few birds below 1 km and too many above 2 km, and between radars the solution follows it. Weighting the profile fit by range (range^−q, q = 1–3, all rows kept) did not fix it and was removed.

| lambda_profile | Bias | Misfit | Trusted misfit | Detail | Prior below 1 km, 0–25 / 100–150 km | Solution below 1 km, 0–25 / 100–150 km |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 (current) | 0.065 | 0.195 | 0.184 | 0.29 | 27.4 / 24.2 | 38.6 / 24.0 |
| 1 | 0.043 | 0.194 | 0.183 | 0.26 | 41.9 / 33.7 | 41.7 / 33.5 |
| 0.1 | 0.045 | 0.194 | 0.183 | 0.26 | 42.3 / 34.3 | 41.4 / 34.1 |

- With `lambda_profile` = 1 the prior matches the solution near the radars, the dip at 100–150 km goes (0.92 → 0.97) and the ring at 25–50 km shrinks, at the same misfit. The low-layer share still falls with distance (0.48 near radars to 0.31 at 250 km, prior and solution alike).

**Correlated prior** (`SB_gp_prior.m`, `figures/gp_prior.png`, `gp_tradeoff_normalized.png`, `.csv`; all with `lambda_profile` = 1). The Gaussian-Markov prior was 10–28× stiffer (sd) than the white prior of the same weight; its precision is now normalized to the same marginal variance (MIGRATION.md).

| Prior | Weight | Bias | Trusted misfit | Detail |
| --- | ---: | ---: | ---: | ---: |
| white | 0.3 / 0.5 / 1 / 2 | 0.043 / 0.035 / 0.026 / 0.016 | 0.183 / 0.199 / 0.260 / 0.354 | 0.26 / 0.23 / 0.18 / 0.13 |
| GP 50 km (normalized) | 0.5 / 1 / 2 / 4 | 0.084 / 0.073 / 0.056 / 0.037 | 0.182 / 0.193 / 0.213 / 0.250 | 0.35 / 0.32 / 0.29 / 0.25 |
| GP 30 / 50 / 100 km, 50 km + vertical 300 m (normalized) | 0.3 | 0.088 / 0.089 / 0.090 / 0.092 | 0.177 / 0.177 / 0.178 / 0.176 | 0.36–0.39 |
| GP 30 / 50 / 100 km, old (unnormalized) | 0.3 | 0.041 / 0.034 / 0.020 | 0.231 / 0.259 / 0.310 | 0.25 / 0.24 / 0.21 |

- At equal trusted misfit the white prior leaves less imprint than the correlated one (e.g. 0.035 at 0.199 against 0.073 at 0.193). The flattening seen with the old GP came from its stiffness, not from the correlation.
- A correlated prior spreads each radar's data level over its neighbourhood: the imprint is a level difference between radars in the data, which a smooth prior on the field cannot tell from real spatial variation. That needs a term for it in the observation model (a per-radar level, as the overlap test measures), not a different field prior.

## Gaussian covariance without a radar prior (research, 2 November 2026)

`SB_gp_model.m`, `figures/gp_model.png` (shared layout `sandbox_case_figure.m`: per run, prior mean and reconstruction as maps and KDMX sections at 02:00, plus the radar-distance curve and the share below 1 km for all runs). `SB_prior_cases.m` (`figures/prior_cases.png`) uses the same layout.

Model: data on every row, + w_p² (x − m)ᵀ Q (x − m), + advection; smoothness off. m is one domain-mean profile per frame (no horizontal structure); Q = c (I + (l_h/dx)² L_h + (l_v/dz)² L_v), normalized to the white prior's variance. Assumptions: stationary deviations in birds km⁻³ (same variance everywhere, not relative to the local level), the same error for every row.

| Case | Bias | Trusted misfit | Detail | Negative mass | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| reference: radar prior (lambda_profile 1), white 0.3, smoothness 0.3 | 0.043 | 0.183 | 0.26 | 0.2% | 15 s |
| GP 150 km × 300 m, flat mean, no smoothness | 0.084 | 0.176 | 0.37 | 2.8% | 114 s |
| GP 50 km × 300 m, flat mean, no smoothness | 0.084 | 0.176 | 0.36 | 1.9% | 88 s |

- The Gaussian covariance fits the data slightly better and keeps more of the low layer between radars (share below 1 km about 0.47 at 100–150 km, against 0.40), but the maps carry more radar imprint and contrast: a dip at 100–150 km (0.92) and a strong patch at 2 km near −98.5°, 41.5–42°N where only far beams look.
- In all three the KDMX section has a dark cone above the radar (near-range rows at the higher elevations reading low); the covariance does not remove it.
- With a stationary covariance a few far, high observations spread over l_h; relative (level-dependent) variance and an observation error for per-radar levels and near-noise rows are the next ingredients to test.

## Two full nights with the production defaults (2 November 2026)

`SB_nights.m`: cases `night_20261021` (sandbox night, 23 frames, 00:30–11:30 UTC) and `night_20210512` (dry spring night, 16 frames, 02:30–10:00 UTC), every 30 min from 1 h after sunset to 30 min before sunrise, production pipeline (gate estimator, no cuts, `lambda_profile` = 1, profiles clipped at zero). Videos and mid-night frames: `<figure_dir>/night_YYYYMMDD/night_review.mp4`, `.png` (column, radar imprint, radar-distance curve, sections through KDMX and KLSX, mean profile near and between radars); scores `figures/nights.csv`.

| Night | Frames | Bias | Column below 1 km near / between radars (birds km⁻²) | Solve |
| --- | ---: | ---: | ---: | ---: |
| 21 Oct 2026 | 23 | 0.047 | 56 / 52 | 29 it, 85 s |
| 12 May 2021 | 16 | 0.024 | 147 / 167 | 29 it, 64 s |

- 12 May with `lambda_profile` 10 (sandbox_run on the same case): bias 0.137, trusted misfit 0.172; with 1: 0.022, 0.169.
- Spring night: towers of high density in the lowest 1 km right above several radars (KDMX, KDVN, KOAX, KLSX; red discs in the imprint map). Autumn night: a dark cone aloft above the radars, a dip of the curve within 10 km, and bright vertical streaks along beams (KLSX section).

**Review of the two nights** (night review now has W–E sections through KDMX and KLSX and S–N sections through KDMX and KILX; every frame in `night_review_frames/`).
1. *Ground steps.* The field is layered in height above sea level (100 m) and the lowest active cell is about as dense as the one above (median ratio 1.01), so where the ground rises by one layer the column loses a layer of birds: between neighbouring columns (4 km) a one-layer ground step changes the column by 8% (21 Oct) and 14% (12 May, dense low layer), against < 0.3% without a step. The maps and imprint maps show blocky contours that follow terrain.
2. *Near-radar elevations disagree.* At 500–1500 m the density within 15 km of a radar is 25–30% below that at 15–60 km on both nights (27.6–29.1 against 39.5–40.2 birds km⁻³ on 21 Oct; 75–79 against 92–94 on 12 May): dark cones above the radars in every section, with narrow bright spikes along single beams inside them (KLSX, KDMX) and, on 12 May, dense towers in the lowest kilometre right above KOAX, KDMX, KDVN, KLOT, KLSX.
3. *Advection.* Radar-centred anomalies (red discs on 12 May at KOAX, KDMX, KDVN, KLOT from 06:00 to 09:30; the KDMX disc at 10:00 on 21 Oct) stay on the radar from frame to frame instead of moving with the wind: each radar's data set its local level at every frame, and the advection term does not carry it across.
4. *Between radars.* Areas dense at two neighbouring radars are lower in between (12 May 06:00: KDMX–KDVN, KDVN–KLOT); the S–N section through KDMX on 21 Oct 05:00 shows a dense layer at 37–39°N and 44–45°N with less in between.

## Column steps at ground-layer boundaries (2 November 2026)

`SB_ground.m`, `figures/ground_sandbox.png`, `ground_night_20210512.png`, `ground.csv`. Column change between neighbouring 4 km columns (mean log2 over frames) where the ground is one 100 m layer higher, and without a step:

| Case | Whole cells: no step / one layer | Air share of the lowest cell: no step / one layer |
| --- | ---: | ---: |
| sandbox (21 Oct 2026) | −0.001 / −0.110 | 0.000 / −0.026 |
| 12 May 2021 | −0.004 / −0.222 | 0.000 / −0.059 |

- The steps come from counting the lowest cell whole although most of it can be below the ground. Counting its share above the ground (`bmm4d.grid.air_fraction`, now in `export_maps`, `radar_distance_bias`, `sandbox_run`, `sandbox_case_figure` and the night review) removes three quarters of the step; the rest is the air really lost as the ground rises.
- Tested and not kept: the share in the forward operator (bias and misfit unchanged to three decimals) and the prior's radar profiles in height above the ground, raw DEM or smoothed 30 km (12 May: step −0.22 → −0.20 with whole cells, −0.059 → −0.056 with the share; bias 0.024 → 0.028).
- Scores from here on count the air share; the sandbox bias reads 0.046 instead of 0.043.

## Near-radar elevations (2 November 2026)

`SB_near_elevation.m`, `figures/near_elevation.png`, `near_elevation.csv`. Solve without the rows within 30 km of their radar, then compare those held-out rows with the prediction (farther beams of the same radar and the neighbours, same heights).

- At the same range the elevations disagree: at 10–20 km, 0.5–1.3° read 1.1–1.4 of the prediction and 1.8°+ read 0.6–0.8 (both nights); at 5 km the highest elevations read 0.2–0.5. Same range means same gate size, so gate size alone does not explain it.
- Against the expected number of birds per gate (λ = predicted density × counting volume) the elevations ≥ 1.5° fall on one curve on both nights: 0.15–0.5 for λ < 1, 0.5–0.75 at 2–4, 0.7–0.9 at 7–13, 0.9–1.1 at 24–42, 1.1 above 75. The gate estimator underestimates when gates hold few birds. High elevations near the radar sample sparse air aloft in small gates (lowest λ): the dark cones.
- Elevations below 1.5° sit 0.2–0.4 higher at the same λ: the low-level towers. Part of it may be the prediction itself near the ground, which few far beams see.
- Reading: at low λ the estimate rests on how many gates are empty. If birds are clumped (flocks: more empty gates than Poisson for the same mean) or the effective counting volume is smaller than modelled, the empty fraction reads too high and the density too low, by a constant factor at low λ, while at high λ the mean echo sets the estimate and it is unbiased. That is the shape seen.

## Clumped birds in the gate estimator (2 November 2026)

Birds per gate negative binomial with dispersion k instead of Poisson; each case cleaned again and solved with and without the rows within 30 km. The option, its script and its products were removed after the test.

| Night | Count model | Held-out near rows < 1.5° / ≥ 1.5° | Bias | Trusted misfit |
| --- | --- | ---: | ---: | ---: |
| sandbox | Poisson / k = 3 / k = 1 | 1.15 / 0.76, 1.15 / 0.75, 1.17 / 0.79 | 0.046 / 0.046 / 0.039 | 0.183 / 0.221 / 0.245 |
| 12 May 2021 | Poisson / k = 3 / k = 1 | 1.13 / 0.82, 1.14 / 0.83, 1.15 / 0.88 | 0.024 / 0.024 / 0.024 | 0.169 / 0.240 / 0.285 |

- Rejected: clumping leaves the low-λ deficit as it was (0.15–0.5 below one bird per gate) and worsens the fit to the 30–120 km rows. A dispersion of 1–3 changes the empty-gate share by a few percent; the deficit is a factor 2–5. 

## Raw reflectivity and estimator by elevation near the radar (3 November 2026)

`SB_near_raw.m`, `figures/near_raw.png`, `near_raw.csv`, `near_raw_fit.csv`. Each volume of the two nights (night_20261021 and night_20210512, 28 radars, one per radar and hour: 571 volumes) cleaned again at full gate resolution (production read settings, rhoHV rain filter instead of MistNet). Per sweep and 100 m of height above the antenna, within 3–60 km: mean raw reflectivity as density (gates without echo count as zero; this mean does not depend on how many birds a gate holds) and mean estimator density, on the same gates. Adjacent sweeps of a volume at the same height (11 977 pairs) differ in range (the upper one is closer) and possibly in waveform (lowest cuts long-PRT surveillance, ~320 Hz, 0.5° radials; upper cuts short-PRT Doppler or batch, ~925 Hz, 1° radials). Fitted jointly: a factor per range bin (1 at 30–60 km) and a factor for a short-PRT upper sweep over a long-PRT lower one.

| | 3–5 km | 5–7.5 | 7.5–10 | 10–15 | 15–20 | 20–30 | short-PRT, VCP 212/215/35 | short-PRT, VCP 31/32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| raw | 0.16 | 0.43 | 0.54 | 0.67 | 0.76 | 0.87 | 0.85 | 1.01 |
| estimator | 0.55 | 1.07 | 1.11 | 1.11 | 1.07 | 1.04 | 0.84 | 0.96 |

- The elevation disagreement is in the radar data. In raw reflectivity the upper sweep reads 0.5–0.75 of the lower one at the same height when it is within 10 km, 0.8–0.9 beyond. Fitted as a range factor it is the radar output's range bias of GD_I (0.18 / 0.47 / 0.63 / 0.74 at 3.5 / 9.5 / 13 / 17 km there). It is not vertical beam smoothing: the upper sweep reads low whether density rises or falls with height.
- The estimator removes the range part from 5 km on (1.04–1.11, a slight over-correction at 5–15 km against 30–60 km; GD_I's estimator curve is lower there, with a different method: a vertical profile through the beam weighting). Within 5 km it reads 0.55.
- A further factor for short-PRT upper sweeps over long-PRT lower ones: 0.84 on VCP 212/215/35, 0.96 on VCP 31/32, varying by radar and night (0.65–1.15). **Correction (3 November, `SB_sweep_offset.m`): this is not the waveform.** On the same beam and gates the short-PRT Doppler pass of the split sweeps reads 0.96–1.02 of the long-PRT pass (median gate difference 0.0 dB). The factor only stands for whatever differs between the lowest and the upper sweeps beyond a range-only curve.
- Link with the held-out test (`SB_near_elevation`): the rows ≥ 1.5° are mostly these short-PRT sweeps, and the inversion has no term for a 16% offset between sweeps of one volume, so near the radar, where the upper sweeps carry the data aloft, it settles between them. The low-λ shape of the held-out ratio is not reproduced by this factor, which does not depend on density.

## Why sweeps of one volume disagree at the same height (3 November 2026)

`SB_sweep_offset.m` (274 VCP 212/215/35 volumes of the two nights, 02–07 UTC; `figures/sweep_offset_example.png`, `sweep_offset_split.png`, `sweep_offset_azimuth.png`) and `SB_processing_era.m` (224 rain-free gate-study volumes of 2015–2025; `figures/processing_era.png`, `processing_era.csv`). Radar output (REF) only, no cleaning.

- What it looks like (KUEX 22 Oct 2026 04 UTC): above about 1 km and beyond 25 km all sweeps agree. Closer in, every sweep reads lower the closer it gets; at 1250 m, 1.3–3.1° at 23–50 km read 79–88 birds km⁻³, 4.0 / 5.1 / 6.4° at 17 / 14 / 11 km read 64 / 44 / 24.
- Waveform (short against long PRT): ruled out. The two passes of the split sweeps (same beam and gates, about 30 s apart) agree: mean ratio 0.96–1.02 at 3–60 km, median gate difference 0.0 dB.
- Point clutter filter (RDA Build 19, 2019–2020): ruled out as the main cause. The range factor from same-height pairs is the same before and after it: 0.20 / 0.31 / 0.40 / 0.50 / 0.63 at 4 / 9 / 12.5 / 17.5 / 25 km in 2015–2018, 0.17 / 0.30 / 0.41 / 0.53 / 0.67 in 2020–2025. Point-flagged gates (CFP code 1) keep their reflectivity, and leaving them out changes nothing.
- Vertical beam smoothing (far beams wider): minor. Five profile shapes through the 0.95° beam and the same pair fit give 0.86–1.0 (0.69 at 4 km for an empty layer below 500 m).
- Noise and detection: not near the radar. The 1st percentile of REF is about −20 dBZ at 10 km on every sweep type, far below one bird in a gate there.
- Clutter filter notch: part of it within about 15 km. The share of gates with clutter power removed rises toward the radar at the same height (KUEX at 1250 m: 0% at 50 km, 26% at 11 km); adding the removed power back lifts 6.4° at 11 km from 24 to 46. On the Doppler and batch sweeps the filtered share peaks where birds fly across the beam (zero radial velocity); on the surveillance pass it is flat at about 60% (aliased velocities). Adding power back also adds ground clutter, so this is an upper bound.
- Bird aspect: reflectivity around the ring is 20–25% above the mean where the beam sees birds from the side (±90–120° from the heading) and 20–25% below head- or tail-on, the same on all sweep types; it does not differ between sweeps.
- Not explained: the range factor keeps rising past 30 km (0.84–0.87 at 37 km, 1.15–1.19 at 52 km against 1 at 30–60 km), which none of the above produces. The bright single-bird gates (> 20 dBZ) become rarer close to the radar, not more frequent as one bird in a small gate should make them.

## Post-estimator filters of clean_volume (3 November 2026)

`SB_post_filters.m`. `clean_volume` caps densities at 3000 birds km⁻³ and removes patches of `data > 0` smaller than 25 km², both designed for the paper's chain where no-echo gates were 0. Checked on 4 mid-night volumes of 12 May 2021 (KDMX, KLSX, KMPX, KICT), cleaned at full resolution with the production options.

- `data > 0` equals `~isnan(data)` on every sweep: the estimator gives every gate with information a positive density (minimum 0.056–0.075 birds km⁻³, below the 0.1 grid floor after the shape correction). The patch filter therefore no longer separates echo from no echo; it can only remove valid islands between missing gates (rain mask, blocked beam, range limits). It removed 0 gates on all 4 volumes.
- Density cap: 33 gates above 3000 at KDMX (0.13% of the bird sum, maximum 8373, so the shape correction can lift values past the 5000 grid top); none on the other 3 (maxima 1288–2744). Not checked whether the KDMX gates are clutter or birds.
- Recommendation: remove the patch filter (or redefine echo as `data >` a threshold if isolated-echo removal is wanted); the cap is rare and harmless.
- Decision: both removed from production (`MIGRATION.md`, "No density cap or patch filter after the estimator"); `clean_volume` only thins.

## Observation sections and KUEX on 21 October 2026 (3 November 2026)

`SB_obs_sections.m` (`figures/obs_sections.png`, `obs_sections_score.csv`): the cleaned, thinned observations as range-height sections (median over azimuth, and over the night for the two full nights), 28 radars × 11 dry nights (the nine screened nights with one volume per radar, cleaned here with the production settings). Score: median log2 of an upper sweep within 20 km over the next lower sweep at the same height. All 28 radars are negative on average (upper sweeps read low): from −0.11 (KSGF) and −0.13 (KDMX) to −0.48 (KIWX) and −0.52 (KABR); it varies more between nights than between radars (KABR −1.45 to +0.05).

`SB_kuex.m` (`figures/kuex_maps.png`, `kuex_section.png`, `kuex_beam.png`): KUEX, 22 Oct 2026 02/04/06 UTC, lowest 1000 m within 30 km.
- Two layers moving in different directions (velocity profile of the volume): north (340–50°) at 4–7 m/s below about 300 m, turning through east at 500–700 m to south (185–190°) at 7–10 m/s above 1 km. The southward layer aloft is the autumn migration; the low northward layer drifts like a southerly near-surface wind (insects or other drifting scatterers are likely; not yet checked).
- The 0.5–1.3° sweeps skim the low layer over 0–30 km and read 15–30 birds km⁻³ there; the profile from the narrow upper beams (3 at 200 m, 100 at 1200 m) passed through their beams predicts 3–6. So the low sweeps carry the low layer, which the upper sweeps cross only within a few km.
- A ground-clutter fan to the south-west on 0.5–1.3° (radial streaks and a strip at 5–15 km), with heavy clutter filtering everywhere within 30 km on 0.5 and 0.9°: adding the removed power back gives 10⁴–10⁵ birds km⁻³ equivalent below 300 m.

## Disagreement cases and the density dependence (3 November 2026)

`SB_disagreement_cases.m` (`figures/disagreement_cases.png`, `disagreement_cases.csv`): the SB_kuex diagnostics on eight radar-nights with a strongly negative score and two with a score ≥ 0. `SB_obs_sections.m` adds `figures/obs_sections_density.png`.

- Across all 301 radar-nights the score follows the density aloft (Spearman 0.62): median score −0.78 (×0.58) below 3 birds km⁻³ at 500–1500 m, −0.52 at 3–10, −0.29 at 10–30, −0.23 at 30–100, −0.10 (×0.93) above 100. The sweeps agree on busy nights and disagree on quiet ones.
- The eight bad cases have 0.3–3 birds km⁻³ at 500–1500 m (KLOT 26-10-21: 18); the two good ones 150–400. In the bad cases the narrow-beam profile is concentrated in the lowest 100–300 m (2–30 times the density at 500–1500 m in six of eight), and the low sweeps read 2–7 times what that profile predicts through their beam (good cases 1.0–1.5).
- In five of eight bad cases the low layer moves in another direction than the air above (85–173° apart; good cases 18–35°). ZDR does not separate the low layer from the air above (medians 0.4–2.9 dB in both).
- Reading: on quiet nights what is left near the ground (local or drifting birds, insects, residual clutter) is a large share of the echo; the low sweeps sample it over 0–40 km, the upper sweeps only within a few km. Gates near the radar on the upper sweeps also hold well below one bird, where the gate estimator reads low (SB_near_elevation).

## Can a shallow dense layer alone explain the score? (3 November 2026)

Synthetic profiles (100 birds km⁻³ at the ground, 0.5 aloft) through the VCP 212 sweeps (two-way Gaussian beam, 0.925°; the part of the beam below ground reads empty) and the same score as `SB_obs_sections` (median log2 upper / lower sweep at the same height, upper within 20 km): exponential layers of scale 100 / 150 / 250 / 400 / 800 m, uniform 0–200 m and 0–500 m, a layer at 1000 m and a deep layer all give −0.03 to 0. Refraction not in the model (beams bent down by a surface inversion, effective earth radius factor 2, 4 and ∞ instead of 4/3; up to 94 m lower at 40 km) gives at most −0.05. The observed −0.3 to −1.5 on quiet nights is therefore not what the geometry gives for any profile: at the same place the sweeps disagree beyond what the forward model can reconcile. (Inline test, not a script.)

## Simulated sparse birds through the observation chain (3 November 2026)

`SB_sparse_sim.m` (`figures/sparse_sim.png`, `sparse_sim.csv`): Poisson birds in every sweep of a VCP 212 volume (0.5–6.4°, 0.5° radials up to 1.3°, 1° above, 2–45 km) with the density of a vertical profile at each gate's beam centre (consistent between sweeps by construction), lognormal cross-sections and the two-way beam weighting of the gate model (`gd_simulate_sweep`), 0.5 dB coding and a detection floor (`gd_observe`), then the production gate estimator; scored as `SB_obs_sections`.

| Profile | Density 500–1500 m | Birds per gate, upper sweeps ≤ 20 km | Score, raw mean | Score, estimator |
| --- | ---: | ---: | ---: | ---: |
| uniform 2 / 10 / 50 / 200 | 2 / 10 / 50 / 200 | 0.19 / 0.94 / 4.7 / 18.8 | −0.03 / −0.03 / −0.01 / −0.01 | −0.06 / −0.01 / −0.01 / 0 |
| quiet / medium / busy, dense low layer | 2.9 / 22 / 117 | 0.38 / 2.4 / 11.5 | 0 / 0 / 0 | −0.02 / −0.02 / 0 |

- Sparse birds alone do not make the sweeps disagree: at 2 birds km⁻³ the simulated score is −0.06, the real radar-nights at that density −0.5 to −1.5.
- Under its own model the estimator reads 0.6–0.9 of the truth at 0.01–0.3 birds per gate and 0.9–1.0 above 1 (azimuth medians); on held-out real rows it reads 0.15–0.5 below one bird per gate (`SB_near_elevation`). Adjacent sweeps at the same height have similar birds per gate, so even that mild bias barely changes the score.
- So the disagreement is in how real radar data differ from the gate model: the echo of sparse, isolated birds in small gates is weaker or missing in the real output (the raw mean shows it too, `SB_near_raw`), which neither the geometry, refraction, nor Poisson sampling produces.

## Real gates against simulated birds, and refraction on quiet nights (3 November 2026)

`SB_gate_compare.m` (`figures/gate_compare.png`, `gate_compare.csv`): density profile from the lower sweeps at 25–45 km (radar output, no echo = 0), birds simulated with it through the real 2.4° (2.5°) sweep, compared with the real gates at 3–25 km. Two quiet (KGRB 2023-05-20, KABR 2025-09-27) and two busy radar-nights (KSGF 2025-10-07, KDMX 2023-05-16).

- Busy: real and simulated distributions overlap (mode 15–25 dBZ, no empty gates); real means 0.8–2.5 of the far profile (the real field is not uniform), with no sign of loss.
- Quiet: the real gates are not what sparse birds look like. At KABR (0.01–0.4 expected birds per gate) the simulation leaves 79–98% of gates empty with a few bright single-bird gates; the real sweep has echo in 28–96% of gates, no gate above 15 dBZ, and a smooth weak mode at −5 to −15 dBZ that falls with range. At KGRB the far sweeps give 60 birds km⁻³ (simulated mode +10 dBZ); the real upper sweep reads 2–5% of it as a weak mode near −10 dBZ.
- The radar's clutter flags on those quiet near gates: point clutter filter 9–21% (busy 2–17%), clutter power removed ≥ 3 dB 34–72% (busy 4–34%).

Refraction (inline test, 12 UTC soundings at the eight radars that are sounding sites, 68 radar-nights): refractivity gradient in the lowest 200 m −5 to −162 N km⁻¹ (ducting below −157: two cases), 0.5° beam at 40 km 351–460 m against 443 m standard. The score is weakly better, not worse, with stronger refraction (Spearman of score and dN/dh −0.29). Refraction and ducting do not explain the quiet-night disagreement.

## Does the estimator's own uncertainty cover the sweep disagreement? (3 November 2026)

Inline test: the gate estimator rerun on every sweep of one quiet volume (KABR 2025-09-28 02:58) and one busy one (KSGF 2025-10-08 02:42), keeping `log_sd` (sd of log density from the curvature of the window likelihood, computed by `estimate_sweep` but not saved by `estimate_density`). Per 1 km range bin: medians over azimuth of density, `log_sd`, the spread of log density over azimuth (IQR / 1.349) and, for upper-sweep bins within 20 km, the offset ln(upper / lower sweep at the same height).

| Radar | Expected birds per gate λ | Rows | log_sd | Spread over azimuth | Offset, ln upper / lower |
| --- | ---: | ---: | ---: | ---: | ---: |
| KABR | < 0.1 | 106 | 1.51 | 0.52 | 0 |
| KABR | 0.1–0.3 | 19 | 0.33 | 0.72 | −1.07 |
| KABR | 0.3–3 | 28 | 0.12–0.15 | 0.36–0.41 | −0.38 to −0.41 |
| KSGF | 3–10 | 40 | 0.08 | 0.56 | +0.08 |
| KSGF | > 10 | 112 | 0.05 | 0.72 | +0.09 |

- Where the gates hold almost no bird (λ < 0.1) the estimator already says it does not know (sd 1.5 in log, a factor 4.5).
- At λ 0.1–3 it is confident (sd 0.1–0.3) while the sweeps disagree by 0.4–1.1 in log, 3–7 times its sd: the near-radar deficit is a bias the likelihood does not contain, not noise.
- Busy: sd 0.05–0.08 and the sweeps agree within 0.1.

## Observation errors from the expected birds per gate (3 November 2026)

`SB_sigma_weights.m` (`figures/sigma_weights.png`, `sigma_weights.csv`, `sigma_weights_validation.csv`, `sigma_pairs.csv`). The solver now takes one weight per observation row and time (`solve_advection_lean` `row_weights` [n_rows n_time], data term Σ w (Fx − y)²); `sandbox_run` builds it from an error model σ(λ) in log density, with λ = (F × radar prior) × counting volume.

Error model, from same-height sweep pairs (gate estimator on every sweep, one volume of each of the ten radar-nights of `SB_disagreement_cases`, 1 157 upper-sweep bins within 20 km): RMS of ln(upper / lower) per half-decade of λ.

| λ | 0.06 | 0.18 | 0.56 | 1.8 | 5.6 | 18 | > 50 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RMS offset | 1.07 | 1.06 | 0.89 | 0.56 | 0.48 | 0.36 | 0.4–0.8 |
| Estimator log sd | 1.01 | 0.57 | 0.33 | 0.14 | 0.09 | 0.06 | 0.03–0.05 |
| σ = (0.4² + 0.5²/λ)^½ | 2.15 | 1.25 | 0.78 | 0.55 | 0.45 | 0.42 | 0.40 |

Solves (weights normalized to mean 1 over the observations):

| Case | Weights | Radar-distance bias | Misfit (trusted) | Detail |
| --- | --- | ---: | ---: | ---: |
| sandbox | none | 0.046 | 0.183 | 0.271 |
| sandbox | (0.4 / σ(λ))² | 0.044 | 0.182 | 0.271 |
| sandbox | 1 / (prediction × σ(λ))² | 0.030 | 0.372 | 0.184 |
| 12 May 2021 | none | 0.024 | 0.169 | 0.258 |
| 12 May 2021 | (0.4 / σ(λ))² | 0.024 | 0.165 | 0.259 |
| 12 May 2021 | 1 / (prediction × σ(λ))² | 0.009 | 0.524 | 0.153 |

Radar cross-validation (sandbox, 7 folds of 4 radars), held-out log bias at 3–15 / 15–30 / 30–60 / 60–120 km: none −0.79 / −0.32 / −0.13 / −0.12; σ(λ) the same to 0.002; relative −0.80 / −0.34 / −0.15 / −0.13.

- The λ weights change almost nothing on these two nights: only 1.2% (21 October) and 3.9% (12 May) of the observations have λ < 1 (6% and 11% of those within 30 km). The near-radar deficit seen in cross-validation (−0.79 at 3–15 km) is at λ above 1, where σ(λ) is already near its floor.
- Relative errors lower the radar-distance bias by weakening the data overall (misfit ×2–3, detail −30%), i.e. a smoother solution; cross-validation is not better. Not a fix.
- The weights should matter on quiet nights (λ < 1 in most near rows), but neither full-night case is quiet.

### Known truth, and range (3 November 2026)

Range: σ depends on range through λ = density × V_sample(r) (V_sample ∝ r²). On the sweep pairs, a robust regression of ln(offset²) on ln λ and ln range gives −0.37 ± 0.03 for λ and +0.05 ± 0.13 for range (3–20 km): at a fixed number of birds per gate, range adds nothing.

`SB_sigma_synthetic.m` (`figures/sigma_synthetic.png`, `sigma_synthetic_maps.png`, `sigma_real.png`, `sigma_synthetic.csv`): truth = the unweighted sandbox solution scaled ×1 / ×0.1 / ×0.02 (busy / medium / quiet); observations on the sandbox rows y = F ρ_true × exp(e), e ~ N(μ(λ), s(λ)²) with the deficit μ = −0.8 / (1 + (λ/1.2)^1.5) and spread s = min((0.2² + 0.2²/λ)^½, 0.8) of the real sweep pairs; the radar prior rebuilt from each synthetic data set.

| Scenario | Rows λ < 1 | Method | Near-radar bias log2(< 30 / 30–150 km) | Total bias log2 | Column log RMSE | Field rel. L2 | Solve (s) | Iterations |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| busy | 3% | perfect / unweighted / σ(λ) | 0.003 / −0.035 / −0.027 | 0.006 / 0.007 / 0.008 | 0.028 / 0.027 / 0.026 | 0.07 / 0.10 / 0.10 | 15.7 / 17.3 / 17.4 | 28 / 29 / 29 |
| medium | 17% | perfect / unweighted / σ(λ) | 0.006 / −0.224 / −0.158 | 0.007 / −0.129 / −0.122 | 0.029 / 0.106 / 0.096 | 0.08 / 0.17 / 0.17 | 18.9 / 17.7 / 19.2 | 29 / 30 / 31 |
| quiet | 49% | perfect / unweighted / σ(λ) | 0.008 / −0.300 / −0.196 | 0.008 / −0.391 / −0.369 | 0.031 / 0.270 / 0.254 | 0.09 / 0.30 / 0.29 | 16.4 / 16.8 / 18.2 | 30 / 32 / 34 |

- The weights remove about a third of the near-radar deficit (−0.22 → −0.16, −0.30 → −0.20), mostly within 15 km where λ is lowest: the dark centres of the radar discs fill in, a ring at 15–30 km stays (λ 1–3, weight 0.6–0.9).
- They do not fix the overall level: when half of the rows read low (quiet), the whole field and the radar prior built from them read low (−0.37 in log2); the far rows with λ < 1 are down-weighted too, so there is nothing unbiased to fall back on.
- Cost: weights from F × prior take 1–2 s; the solve takes 0–1.5 s longer (0–2 more PCG iterations, +0–9%).
- Real nights (both busy, 6% and 11% of rows within 30 km below λ = 1): the weighted column differs from the unweighted by at most 3.5% (sandbox, within 10 km) and 1.2% (12 May); radar-distance bias 0.046 → 0.044 and 0.024 → 0.024, solve time unchanged.

## The estimator's own uncertainty against the misfit of the inversion (3 November 2026)

`SB_estimator_sd.m` (`figures/estimator_sd.png`, `estimator_sd.csv`, `estimator_sd_solves.csv`, `estimator_sd_validation.csv`). The gate estimator now returns `sd`, the spread of each window's likelihood in birds km⁻³ (MIGRATION.md). The 168 sandbox volumes were cleaned again with it (`prod/volumes/sandbox_sd`) and the observation rows built a second time from sd (`prod/sandbox_sd`), so every observation has its sd. On simulated uniform sweeps the sd is 0.3–0.6 of the actual scatter of the estimate over azimuth (overlapping radials and windows): the likelihood is too narrow by about 2.

Sandbox night, unweighted solution, by expected birds per gate λ = prediction × counting volume (8.5 million observations, 97% with λ > 1):

| λ | Share of rows | RMS residual y − Fρ | RMS sd | sd² / residual² | Median sd / y | Median \|residual\| / y |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.02–0.1 | 0.3% | 1.1–1.9 | 0.9–1.0 | 0.28–0.69 | 1.2–3.3 | 1.0–2.2 |
| 0.1–1 | 2.4% | 3.1–4.9 | 1.2–1.5 | 0.10–0.16 | 0.37–0.66 | 0.50–0.78 |
| 1–10 | 15% | 5.0–6.1 | 1.8–1.9 | 0.09–0.14 | 0.15–0.21 | 0.20–0.30 |
| 10–100 | 52% | 5.4–6.1 | 2.1–2.3 | 0.14–0.16 | 0.09–0.12 | 0.14–0.17 |
| > 100 | 31% | 8.3–16 | 2.4–3.0 | 0.04–0.08 | 0.01–0.05 | 0.13 |

- The sd is not negligible: it explains 10% of the residual variance overall, about 40–60% after the ×2 for overlap, and more where windows are nearly empty. The rest (a relative misfit of about 13% on busy rows) is what the reconstruction cannot represent.
- In absolute units, which the least squares uses, both the residual (5–8 birds km⁻³) and the sd (1–3) change little with λ: the equal weights of the current inversion are already close to what the sd implies on this night.

Correction (same day): the first sd collapsed to near 0 where a window's likelihood is narrower than the density grid step (2.8% of rows, sd < 0.1 at a median density of 33); those rows got up to 100 times the median weight and made the first weighted solves look like they filled the near-radar cones. `estimate_sweep` now takes at least density × log_sd there. The residual table above used the first sd; with the corrected one the RMS sd is 1.6–3.8 birds km⁻³ (λ > 0.5) and sd²/residual² 0.10–0.20 between λ 0.5 and 60 (all rows: RMS sd 2.7 against RMS residual 7.3, 14% of the variance). The solves below use the corrected sd.

Solves weighted by 1/σ², weights normalized to mean 1 (so a common factor on sd changes nothing): "sd" σ = sd; "2 sd + 13% model" σ² = (2 sd)² + (0.13 × prior prediction)².

| Weights | Radar-distance bias | Misfit (all / 30–120 km) | Detail | Total birds | Solve (s) | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 0.046 | 0.194 / 0.183 | 0.271 | 1 | 19.1 | 30 |
| sd | 0.051 | 0.440 / 0.405 | 0.210 | 0.93 | 36.7 | 67 |
| 2 sd + 13% model | 0.045 | 0.429 / 0.391 | 0.192 | 0.95 | 35.8 | 68 |

`figures/estimator_sd_frames.mp4` / `.gif` / `.png`: column map and sections through KDMX, KLSX (west-east) and KDVN (south-north) for the three solves, every frame.

- The sd in birds km⁻³ grows with the density read (about 10% of it), so 1/sd² weights low readings far more than high ones. Dense structures (the plumes west of KLSX, the low layer south of KDVN) are smoothed away, the dark cones above the radars stay, and the night loses 7% of its birds: low readings, including the low near-radar windows, now pull hardest.
- So the estimator's sd, used as it is, makes the reconstruction less realistic, not more. Two things are missing: the error must not follow the value read (evaluate the likelihood's spread at the density the inversion predicts), and the prior, smoothness and advection terms were tuned for equal weights.
- The cross-validation above was run with the first (wrong) sd; it is not repeated.

Then with the sd at the prediction (same day): the reading-dependence removed by one reweighting step. From all observations of the night, the estimator's median sd / density by expected birds per gate λ and range (5–20% at λ 1–10, 3–5% above 30; larger at long range for the same λ) is read at λ = (unweighted prediction) × counting volume and multiplied by that prediction. Variants now: unweighted, sd at the reading, sd at the prediction (the "2 sd + 13% model" variant dropped).

| Weights | Radar-distance bias | Misfit (all / 30–120 km) | Detail | Total birds | Solve (s) | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 0.046 | 0.194 / 0.183 | 0.271 | 1 | 16.1 | 30 |
| sd at the reading | 0.051 | 0.440 / 0.405 | 0.210 | 0.93 | 31.0 | 67 |
| sd at the prediction | 0.047 | 0.427 / 0.388 | 0.202 | 0.95 | 33.1 | 69 |

- Taking the sd at the prediction gives back 2 of the 7% of birds lost, but the picture is the same: the KLSX plumes and the dense low patch south of KDVN are smoothed away and the cones above KDMX and KDVN stay.
- The cause is that the estimator's error is relative (a few % of the density where gates hold many birds): in birds km⁻³ dense observations get weights 10–100 times smaller than sparse ones, and the prior, smoothness and advection, tuned for equal weights, take over where the birds are.

Then the scaling of the weights (same day). The prior, smoothness and advection strengths were tuned with every weight 1. Scaled to mean weight 1, a weight-1 observation stands for σ = 1/√mean(1/σ²) = 0.47 birds km⁻³ (set by the near-empty windows, which then carried 65% of the summed weight), while the unweighted solve misses by 7.3 and the RMS sd is 2.7: typical observations (10–30 birds km⁻³, sd 1.5) got weight 0.1 against the regularization, hence the smoothing. Scaled to mean variance, w = mean(σ²)/σ², a weight-1 observation stands for the RMS σ. All with σ at the prediction; "floor" adds an absolute error, σ² = sd² + floor².

| Weights | Radar-distance bias | Misfit (all / 30–120 km) | Detail | Total birds | Solve (s) | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 0.046 | 0.194 / 0.183 | 0.271 | 1 | 14.9 | 30 |
| sd, scaled to mean weight (old) | 0.047 | 0.427 / 0.388 | 0.202 | 0.95 | 29.1 | 69 |
| sd, scaled to mean variance | 0.051 | 0.295 / 0.250 | 0.282 | 0.96 | 77.7 | 187 |
| sd + floor 1, mean variance | 0.053 | 0.274 / 0.234 | 0.278 | 0.98 | 22.3 | 47 |
| sd + floor 3, mean variance | 0.050 | 0.244 / 0.211 | 0.271 | 0.99 | 16.7 | 33 |

- With the variance scaling the structures come back (detail as unweighted): the KLSX plumes and the low patch south of KDVN are there, with lower peaks, and the low layers are smoother. The floor moves the solution towards the unweighted one and makes the solve as fast as unweighted (floor 3: 33 iterations, against 187 without floor, where near-empty windows get weights up to 230).

### Speed, and the streaks west of KLSX at 03:00 (3 November 2026)

Speed. PCG needs more iterations the wider the spread of the weights (the system is worse conditioned): mean-variance scaling without floor gives near-empty windows weights up to 230 and takes 187 iterations; the floor caps them. Floor 1 / 1.5 / 2 (σ from the unweighted solution, warm start): 41 / 35 / 32 iterations, 20.5 / 16.2 / 14.4 s; warm start saves only 5 iterations (floor 1.5 cold: 40). Taking σ from F × the radar prior instead of the unweighted solution needs one solve instead of two: floor 1.5, 37 iterations, 17.7 s (unweighted 15 s, 30), the field 6% different from the two-solve version. `SB_estimator_sd.m` now solves unweighted, "sd + floor 1.5, prior" and "sd + floor 1.5, solution" (the animation shows these three):

| Weights | Radar-distance bias | Misfit (all / 30–120 km) | Detail | Total birds | Solve (s) | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 0.046 | 0.194 / 0.183 | 0.271 | 1 | 15.5 | 30 |
| sd + floor 1.5, σ from the prior | 0.045 | 0.230 / 0.205 | 0.274 | 0.99 | 18.2 | 37 |
| sd + floor 1.5, σ from the unweighted solution | 0.053 | 0.267 / 0.228 | 0.275 | 0.98 | 19.5 (+15.5) | 40 |

Streaks (`figures/klsx_streak_0300.png`, inline). The bright vertical stripes at 1.5–3 km near −92.3° at 03:00 come from KLSX's lowest sweep (0.48°) at 136–138 km, azimuth about 270°: observations of 200–370 birds km⁻³ against 30–80 around, rising through the night (6 at 00:30, 103 at 02:30, 225 at 03:00). In the radar data the box reads 20 dBZ (16.5–17.5 in the whole western sector at 60–145 km), RHO 0.70–0.76, ZDR 0.4–1.3 dB, smooth radial velocity; the 0.8° sweep sees it at 15 dBZ (beam centre about 3.3 km). No gate of the western sector reaches RHO > 0.95 from 02:30 to 07:00, so it is not rain; MistNet and the estimator keep it, as they should for a slightly denser part of the layer. The stripes are made by the inversion: at 136 km the 0.48° beam is about 2.2 km deep, nothing else sees that place (edge of KLSX coverage, the neighbours beyond 150 km), and the thinned observations 2 km apart each become a column of their own.

### Why the sections follow the lowest beam: batch cuts read low

`SB_cut_type_deficit.m` (`figures/cut_type_deficit.png`, `cut_type_deficit.csv`); details in the gate-density report, Part III. The V-shaped sections (a dark cone above each radar, density along the lowest beams) come from the observations, not from the weights. At 10–30 km the batch cuts (≥ 1.8°) read about 0.7 of the split cuts at the same height in the radar output, beyond the near-range deficit that both share. The batch cuts are clutter-filtered at nearly every gate, the split cuts mostly not; putting the removed power back closes the gap. The estimator keeps the offset (batch/split 0.62–0.73 at 14–22 km). Diagnostic solves on the sandbox night (01:30 sections in `figures/near_range_bias_sections_*.png`, `range_test_sections.png`): without rows within 30 km the cone fills; dividing the batch cuts by their deficit or adding it to the observation error does not remove it.

### The clutter filter's removed power as an observation error

`SB_cfp_error.m`, dev only (needs `dev/cfp_error_toolbox.patch`; `figures/cfp_error_sections.gif`, `cfp_error_sections_0130.png`, `_0300.png`, `cfp_error_weights.csv`). After `SB_cut_type_deficit.m` and gate_density GD_L showed that the removed power (CFP) is neither all birds nor all clutter, the estimate is kept and the filter only lowers weights. `sd_cfp` is the estimate with filtered gates as [output, output + CFP] and point-clutter gates as "at least the output", minus the production estimate; it is added in quadrature to the sd-based sigma (floor 1.5), with the weight scale taken without it.
- **Volumes:** re-cleaned in `prod/volumes/sandbox_cfp`, rows in `prod/sandbox_cfp`. The density is identical to the production cache; 40% of observations have `sd_cfp` > 0.
- **Weight against production (median):**

  | Cuts | 0–10 km | 10–20 km | 20–30 km | 30–40 km | 40–60 km |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | batch ≥ 1.8° | 0.07 | 0.15 | 0.47 | 0.69 | 0.89 |
  | 0.5° | 0.001 | 0.004 | 0.02 | 0.10 | 0.95 |
  | 0.9–1.3° | 0.01 | 0.07 | 0.60 | 0.96 | 1 |

  The 0.5° cut loses its weight out to 40 km, because ground clutter is filtered there.
- **Solve:** 37 iterations and 16 s, as production. Birds ×1.004.
- **Sections:** the cone above each radar is filled by ×1.4–2 below about 1.5 km, and KLSX's bright column at the radar goes. The dark wedge above KLSX is lighter but still there, and the bands at 30–90 km (in the data) stay.
- **Column / regional mean** within 25 km of a radar: 0.93 → 0.95. Detail 0.274 → 0.260; trusted misfit 0.205 → 0.223.

### Levels per radar, cut type and distance in the 3D inversion (dev, 4 November 2026)

`SB_cut_levels.m` (about 4 min), the 3D version of gate_density GD_M. Figures: `figures/cut_levels_sections_0130.png`, `_0300.png`, `cut_levels_sections.gif`; table `figures/cut_levels.csv`. Every observation is read as y = c(radar, cut type, distance) × β(λ) × (Fρ).
- **c** is linear in distance between knots at 6 / 15 / 25 / 37.5 / 62.5 km. It is kept within [P_out/P_all, 1.1] of that radar's cuts of that type over the night (`gd_cfp_power`) and pulled toward 1 by one sweep's worth of observations. Split cuts at 62.5 km are the anchor.
- **β** is the sampling curve learned on the 124 reference scans (GD_M), fixed here.

It is solved in rounds: the field with y/f and weights w f², with the prior rebuilt from y/f; then the levels from the new field. Round 0 is production, followed by three rounds. The rounds take 32 / 21 / 16 / 15 s (warm start).
- **Batch levels**, median over radars at 6 / 15 / 25 km: 0.63 / 0.88 / 0.95. KMPX (22% filter loss) 0.94 / 1.07 / 1.06; KDDC (86%) 0.79 / 0.99 / 1.02, mostly clutter; KVNX (70%) 0.51 / 0.47 / 0.66.
- **Against production weights:**
  - column within 25 km of a radar 0.93 → 0.98 of its regional mean;
  - radar-distance bias 0.045 → 0.035;
  - misfit at 30–120 km 0.205 → 0.201;
  - detail 0.274 → 0.264;
  - total birds × 0.995.
- **Sections:** the cone above KLSX, KDVN, KDMX and KVNX rises × 1.5–2 below about 1.5 km and is lighter, though not gone. KMPX barely changes. KVNX gets speckle at 2–3 km.

**Full night** (`case_name = "night_20261021"`: 23 frames, 00:30–11:30 UTC, run as a batch job with 4 workers in 11 min; `figures/cut_levels_night_20261021.csv`, `cut_levels_sections_night_20261021_0600.png`, `_1130.png`, `.gif`). This case was cleaned before the estimator returned its sd, so round 0 uses equal weights.
- **Solves:** 93 / 119 / 107 / 91 s for 29 / 19 / 16 / 14 iterations, 410 s in total against 93 s for one solve. The client ran the morph study at the same time, so the per-round times are probably inflated.
- **Batch levels**, median at 6 / 15 / 25 km: 0.81 / 0.93 / 0.99 after round 1, 0.72 / 0.90 / 0.98 after round 2. They were still moving at 6 km, so three rounds is not converged.
- **Against round 0:**
  - column within 25 km of a radar 0.91 → 0.95;
  - radar-distance bias 0.047 → 0.034;
  - misfit at 30–120 km 0.182 → 0.181;
  - detail 0.263 → 0.255.
- **Sections:** the cones above KLSX, KDVN, KDMX and KVNX rise × 1.5–2 below about 1.5 km, as on the sandbox night.

### Smoother animations: advection morph between frames (dev, 4 November 2026)

`SB_morph_animation.m` (`figures/morph_night_20261021.mp4`, `.png`). The frame at fraction τ of an interval is (1 − τ) A(ρ_t, +τΔt) + τ A(ρ_t+1, −(1 − τ)Δt):
- A moves each layer along the interval's mean layer wind, interpolated from the radar profiles as in the advection term (semi-Lagrangian, two-stage departure point, bilinear).
- Cells below ground take the lowest active value of their column first, so that terrain steps do not pull in zeros.
- It equals the frames at τ = 0 and 1. Structures travel with the wind instead of cross-fading, and change the wind does not explain fades in.

It takes 0.24 s per intermediate frame on the 317 × 282 × 30 grid (132 frames at 5 min for the full night). It moves the radar-fixed artefacts (cones) with the wind between frames.

**us_east** (`SB_us_east_cut_levels.m`: production cleaning of the night, then `SB_cut_levels` with `case_name = "us_east"`; 101 of 107 radars, 26 frames, 5.8 million rows; figures `cut_levels_sections_us_east_0530.png`, `_1200.png`, `.gif`, table `cut_levels_us_east.csv`). Run in the agent's MATLAB session, 4 November 2026, restarted once after the session was closed (cleaned volumes kept in the cache).
- **Time:** cleaning 63 min (776 volumes; 1616 cached), static clutter 9 min, observations 7 min, clutter-filter power 15 min (4 workers). Solves 805 / 488 / 494 / 388 s for 47 / 35 / 28 / 24 iterations (round 0, then three rounds of levels), about 36 min against 13 min for one solve.
- **Batch levels**, median over 101 radars at 6 / 15 / 25 km: 0.60 / 0.81 / 0.92, as on the sandbox and the full Midwest night.
- **Against round 0** (production weights):
  - column within 25 km of a radar 0.91 → 0.98 of its regional mean;
  - radar-distance bias 0.044 → 0.020;
  - misfit at 30–120 km 0.333 → 0.322;
  - detail 0.359 → 0.343;
  - total birds × 0.996.
