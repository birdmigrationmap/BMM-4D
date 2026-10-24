# Production test (dev, exploratory)

Purpose: find out what it takes to turn the method into an operational product, meaning animated 3D bird-density maps for a radar network that are robust to rain and anomalies, with realistic regularization, and fast and light enough to run routinely. It reuses the shared `functions/` and keeps its own code in `prod/functions/`. The paper workflows (`01_`–`11_`) stay in the paper repository `BMM-within` with their own frozen copy of the code (tag `paper-v1` here).

## Layout

| Path | Content |
| --- | --- |
| `functions/prod_case.m` | All case settings (radars, times, grid, cleaning, inversion). Edit here to try settings. |
| `scripts/P1_get_volumes.m` | Find the raw volumes, or download them in parallel over HTTPS into the local raw cache. |
| `scripts/P2_clean_volumes.m` | Clean every volume once (MistNet by default) into a compact volume cache, in parallel. |
| `scripts/P3_build_case.m` | Grid, then forward-operator rows per radar (in parallel), stacked into `case.mat`. |
| `scripts/P4_invert.m` | Density prior, advection mappings, joint inversion (paper full-grid solver for now). |
| `scripts/P5_plot.m` | Nine-time column-density map and MP4. |
| `sandbox/` | Small Midwest case (28 radars, 6 frames) for trying solver options in seconds to minutes; `sandbox_run` scores each option. |
| `experiments/two_level/` | The two-level separable model (23–24 October 2026), set aside: code, scripts, figures and results. |
| `gate_density/` | Bird density from the distribution of gate values under local homogeneity (instead of outlier removal and median smoothing); tests and the near-range deficit by method. |
| `cleaning_study/` | Screening of ρHV vs MistNet on 300 random volumes and the 24-case review gallery. |
| `legacy_us_night/` | First full-US scripts (superseded by P1–P5; kept for reference). |

Each script reads a `case_name`: `keywest` (fast test), `midwest` (17-radar network), `us_east` (107 radars, one night). Add `_rho` or `_none` to change the rain removal, e.g. `keywest_rho`. A `case_name` already set in the workspace is kept, so a batch can loop over cases.

```matlab
project_root = setup_paths();
case_name = "keywest";
for f = ["P1_get_volumes" "P2_clean_volumes" "P3_build_case" "P4_invert" "P5_plot"]
    run(fullfile(project_root, "prod", "scripts", f + ".m"))
end
```

## Data

- Raw volumes: `nexrad_raw_dir()` (local cache outside OneDrive), plus the old project folder `data/raw_scans/nexrad`.
- Cleaned volumes: `local_data_dir()/prod/volumes/<hash>/`. The hash comes from the cleaning settings (`settings.json` in each folder), so the volumes are shared by all cases with the same cleaning.
- Case products: `local_data_dir()/prod/<case>/` holds `manifest.mat`, `case.mat`, `inversion.mat` and `stage_log.csv` (time, duration and client memory of each stage).
- Figures: `prod/figures/<case>/`.
- Worker memory was sampled every 2 s by an external `ps` loop into `local_data_dir()/prod/memory_samples.csv`, because a client snapshot misses parfor peaks.

## Approach

The goal is a coherent, plausible density field, produced quickly. There is no ground truth, so matching the paper results exactly is not a goal. Removing aberrant echoes matters most; the reconstruction can be approximate.

Current defaults (`prod_case.m`):
- **Cleaning:** MistNet within 150 km (its Cartesian grid ends there), with the rain mask grown by 2 km. Gates without echo are kept as zero-density observations (`noEchoRange`), not as missing data; blocked and rain gates stay missing. `read_rd` runs at full gate resolution, then `clean_volume` removes densities above 3000 birds km⁻³ (about 40 dBZ) and connected echo patches smaller than 25 km² (point targets, rain fragments, speckle), and thins the gates to 1 km (2 km for `us_east`).
- **Active domain:** columns within 250 km of a radar (more than half the radar spacing, so there are no gaps between radars), cells from the ground to 3 km above it.
- **Solver (`P4`, `solver_path = "lean"`):** `solve_advection_lean`. It uses the paper's objective and term weights (data, prior, 3D smoothness with the land–water cut, advection), on active cells, in single precision, with all times batched, and an approximate tolerance of 1e-4. Each advection link is weakened where the density prior changes along it (`advection_relax_exponent` 2), so takeoff and landing are not forced to conserve birds. `solver_path = "paper"` keeps the full-grid double-precision path.
- **Night window:** `prod_case` keeps the frames where the sun (`sun_elevation`) is below −6° (civil twilight) at one radar at least, and `radar_night` drops each radar's volumes while the sun is above −6° there (P1 does not download them). `P5` hides cells where the sun is up, fading over ±1°. Midwest: 19 frames (01:30–10:30 UTC) and 285 of 323 volumes; Key West unchanged; `us_east`: 2556 of 2782 volumes.
- **Solve grid:** `solve_aggregate` [x y z] groups kernel-grid cells (2 km × 100 m) into blocks. `construct_F(..., aggregate=...)` sums the kernel weights of each block into one column, so the fine grid is never built as columns; `coarsen_grid` gives the matching grid (block centres, majority water mask, mean DEM), and the velocity profiles are averaged to the coarse layer spacing. Default [1 1 1]; `us_east` and the `_coarse` suffix use [2 2 2] (4 km × 200 m, 8× fewer unknowns).
- **Maps (`P5`):** the density fades to 40% opacity between 150 km (end of the radar data) and 250 km (edge of the solved domain), where it comes only from the prior and advection.
- **Density prior:** `build_density_prior_fast`. Each radar gets a vertical profile per time (100 m layers), fitted to its own observations through the beam kernels with vertical and temporal smoothness; the profiles are fitted on land and water cells. Each column then takes a Gaussian-weighted average of all radar profiles within 3 L (L = 150 km). The paper interpolation (`prior_interpolation = "idw"`: 3 nearest radars, inverse distance, exact at each radar, land-only fit) reproduces the paper function to 1e-15; it leaves seams where the nearest-radar set changes and a peak at each radar, which become visible map features away from the radars (`figures/us_east/prior_review*.png`).
- **Domain (`us_east`):** a fixed geographic domain (`build_domain_mask`): contiguous US land east of 104°W and water (sea, Great Lakes) within 50 km of it. Cells far from radars are solved too and shown without fading. `create_grid` now classifies cells north of the land–water raster (47.02°N) from the Mapping Toolbox land and lake polygons.

## Findings (22 October 2026)

**The pipeline is equivalent to the paper path.** With ρHV cleaning and the paper solver, prod reproduces the published Key West result (median column correlation 0.992, same mean). The fast prior matches the paper function to 1e-15, the batched full-grid solver matches exactly, and the lean solver with every cell active and in double precision matches to 2e-5.

**MistNet needs `maxRange` = 150 km.** MistNet gives no label beyond its ±150 km grid, so rain cells at 150–160 km east of KBYX stayed in and became a large blob on the grid edge. Within 150 km, MistNet removes the rain speckle that ρHV leaves over Florida and keeps the Cuba pulse. On the Midwest it agrees with ρHV until 07:30 UTC and removes more of the north-west rain cells afterwards.

**The patch filter removes point targets.** Moving point targets south of KBYX (05:30–07:00, possibly ships) are gone. Two compact cells near KBYX at 01:00 (above 25 km², probably showers) remain.

**Paper Key West solve.** The paper's Key West solve stops at the 300-iteration cap (relative residual 3e-4, flag 1); Appendix D states 10⁻⁶. This is flagged for the paper, not changed here.

**Where the time goes** (36 GB MacBook; memory from an external `ps` sampler):

| Case (radars × times, unknowns) | Clean (MistNet) | Rows | Prior | Solve | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Key West (1 × 25, 20 M), paper solver | 75 s, 6 workers | 13 s | 1.5 s | 200 s (300 it, cap) | ~5 min |
| Key West, lean | 129 s, 4 workers | 30 s | 0.6 s | 105 s (300 it, cap) | ~4.5 min |
| Midwest (17 × 22, 131 M), paper solver | 947 s, 4 workers | 92 s | 246 s | 1206 s (161 it) | ~42 min |
| Midwest (121 M active), lean | 1000 s, 4 workers | 81 s | 64 s | 436 s (90 it) | ~26 min |

- **Cleaning memory:** about 1.6 GB per worker with ρHV and 2.6–3.7 GB with MistNet, so MistNet runs on 4 workers. Idle pools keep their memory, so P2 and P3 shut the pool down when they finish.
- **Solver per iteration (Key West):** paper 0.73 s; batched full grid 0.52 s; lean double, all cells 0.54 s; lean double, active cells 0.35 s; lean single, active cells 0.27 s. Single-precision sparse products work in R2026a (1.5× faster than double, relative difference 7e-8).
- **Advection term:** it is now the largest part of each iteration (0.20 s of 0.5 s). A stacked sparse operator is 1.5× faster but needs GB of memory, so `accumarray` is kept.
- **Iterations still matter:** at 100 iterations the Key West map differs from 300 (correlation 0.987, mean +4%). The weak Key West prior (0.003) converges slowly.

## Aberrant echoes (22 October 2026, cleaning_study SC_D)

`SC_D_filter_variants` cleans the volumes behind each aberration and the 24 gallery cases with the current settings and candidates (figures in `cleaning_study/figures/variants/`).

- **Rings around KILX and dark discs around radars after sunset were not echoes.** Two causes in the method:
  - *No-echo gates were missing data.* `read_rd` returned NaN for gates below the detection threshold, so beams above the bird layer gave no information and coverage ended where the echo ended (e.g. 1.5° stopped at 45 km at KILX at 01:30). They are now zero-density observations up to the altitude cap. Rain, blocked and removed gates stay missing. Level II flags no-echo and range-folded gates with the same code, but the surveillance cuts do not fold within 3.5 km altitude.
  - *Advection forced conservation during takeoff.* At 01:30 the radars and the prior both show almost no birds (prior 1–20 birds km⁻²), but the solved field rose to 290 birds km⁻² at 180 km from KDMX: the advection term pulled the first frames toward the dense later frames wherever the beams did not reach. With the relaxed links, the 01:30 column is 4–27 birds km⁻² at all ranges and the rings are gone; mid-night maps change by less than 5% (domain mean 239 → 233 birds km⁻²); dawn frames decline more cleanly.
- **The KBYX "shower" cells look biological.** They have noisy ZDR (−2 to 4 dB) and mixed ρHV like the surrounding birds; the real rain cells east of them (ρHV ≈ 1, ZDR ≈ 0) are removed by MistNet. They are kept.
- **A ρHV cap on top of MistNet is not safe.** Removing gates with median ρHV above 0.90 or 0.95 removes rain remnants at KBUF, but also most of the low-signal bird echo at KBYX and KILX after sunset, where bird ρHV is 0.85–0.95. Kept as the opt-in `removeRain = "mistnet_rho"` in `read_rd`, not used.
- **A local-contrast filter** (gates more than 12 dB above the 20 km median) flags 0.01–0.09% of gates in clean volumes and speckle inside rain; not adopted (`cleaning_study` scratch test only).
- **KFSD after 09:30 looks biological.** From 04:40 to 10:30 the echo next to the radar keeps the same signature: ρHV peaking at 0.75–0.85 for dz above 15 dBZ, positive and noisy ZDR, and the azimuthal ρHV pattern of birds (0.3–0.5 to the north-east, 0.7–0.9 to the south). The widespread rain to the west (ρHV ≈ 1) is removed by MistNet. Density near the radar grows as the front approaches (up to 1000 birds km⁻³ below 500 m), consistent with birds concentrating ahead of it. Kept; the post-sunrise frames (11:30) are now outside the night window.
- **Robustness:** `read_rd` failed on volumes almost entirely removed as rain (the shape fit had too few azimuths). It now skips the shape correction with a warning.

Shared `read_rd` options added (defaults leave the paper results unchanged): `noEchoRange`, `rainDilation`, `removeRain = "mistnet_rho"`, and the shape-fit guard. Tests: 52/52 pass.

Timings with these settings: Key West ~3.5 min (clean 78 s, rows 24 s, solve 116 s); Midwest ~22 min (clean 822 s, rows 55 s, prior 56 s, solve 315 s, 73 iterations). The zeros add rows (Midwest 3.1 M) but not solve time.

## Coarse solve grid (22 October 2026)

Midwest, fine (2 km × 100 m) vs `midwest_coarse` (4 km × 200 m), same cleaned volumes and night window, compared on 4 km columns:

| | Fine | Coarse |
| --- | ---: | ---: |
| Unknowns | 105 M | 13 M |
| Rows (P3) | 98 s | 44 s |
| Prior | 84 s | 36 s |
| Solve | 344 s (94 it) | 72 s (65 it) |

Column correlation per frame 0.99–0.997 (0.98 at the first and last frames), mean column density 2–3% lower on the coarse grid (8–12% at the first and last frames, where the layer below the birds is thicker).

## First full eastern-US night (21–22 October 2026)

`us_east`: 107 radars east of −104°, 23:30–12:00 UTC every 30 min, night window at −6° sun, MistNet cleaning with 2 km gates, solve on 4 km × 200 m.

| Stage | Time | Notes |
| --- | ---: | --- |
| P1 download | 131 s | 2393 of 2556 night volumes (968 new, ~7 GB); no data for KHGX, KILN, KJAN, KLIX, KOUN, KRMX |
| P2 clean (MistNet, 4 workers) | 5588 s | 2392 cleaned, 1 failed |
| P3 rows | 227 s | grid 900 × 745 × 15 (fine 1799 × 1489 × 30); 101 radars; 5.9 M rows, nnz(F) 4.9e7; 145 sweeps without kernel |
| P4 prior + steps | 106 s | |
| P4 solve | 399 s | 154 M unknowns (5.9 M active cells × 26), 68 iterations to 1e-4 |
| Total | ~1 h 45 min | cleaning is 89% of it; client peak ~8 GB; `case.mat` 1.2 GB |

The map (`figures/us_east/`) is coherent: takeoff after dusk moving west across the frames, the main flow over Nebraska–Iowa–Missouri and the lower Mississippi valley, and landing before dawn moving east to west. Visible issues: small bright spots at single radars (KBRO, west Texas speckle), scalloped edges where the fade around each radar's 150 km data range overlaps, and gaps filled from neighbours where radars had no data.

## Eastern US at 4 km × 100 m with the geographic domain (22 October 2026)

Cleaning with the MistNet screen (`mistnetScreen`: MistNet only on volumes with rain-like echo, about 42% of volumes skipped in the study set): 3981 s instead of 5588 s for the night (1.4×). The continuity filter in `read_rd` now computes medians only where there is echo (identical output, about 10% faster single-threaded).

| Stage | Time |
| --- | ---: |
| P3 grid / rows / save | 72 s / 286 s / 76 s (900 × 745 × 30, nnz(F) 9.2e7) |
| P4 prior (Gaussian, L = 150 km) | 260 s |
| P4 solve | 1062 s, 88 iterations, 238 M unknowns (9.2 M active cells × 26) |

## Solver review (22 October 2026)

**Objective** (`solve_advection_lean`, same terms and scaling as the paper solver). For frames t = 1…T with unknown densities x_t on the active cells:

J = Σ_t s_t² ‖M_t (F x_t − d_t)‖² + Σ_t p_t² ‖x_t − x̄_t‖² + Σ_t q_t² ‖L x_t‖² + Σ_t a_t² ‖W_t (x_{t+1} − A_t x_t)‖²

- Data: M_t keeps the observed rows (no-echo gates are observations of zero), s_t = 1/√n_t so every frame has the same total data weight whatever its number of observations.
- Prior: x̄_t the density prior, p_t = λ_prior · β_t / √N (β_t = ‖s_t F‖_F, the data block norm; N active cells).
- Smoothness: L the 3D Laplacian on active cells, cut at land–water links, neighbour weights (1 + h)/6 horizontally and (1 − 2h)/6 vertically (h = `coeff_horver`); q_t = λ_smooth · β_t / ‖L‖_F.
- Advection: A_t moves each cell along the radar wind profile to the next frame, W_t the link weights relaxed where the prior changes (takeoff, landing); a_t = λ_adv · mean(β_t, β_{t+1}) / ‖A_t‖_F.
- Solved by preconditioned conjugate gradients (diagonal), all frames at once, single precision, to a relative residual of 1e-4. No non-negativity constraint: negative densities are clipped to zero only in the maps.

Because every term is scaled by the data norm, the λ are relative strengths, but the scaling depends on the grid (N, ‖L‖_F) and on the number of observations: the same λ is a different balance at 2 km and 4 km, or with and without no-echo zeros.

**What the artefacts are** (us_east, first full run):
- *Negative densities:* 13% of cells, down to −1464 birds km⁻³, 1.8% of the positive mass; 21% of cells in the first frame, about 9% at midnight. Clipping them in the maps leaves holes and rings.
- *Beam imprint ("V" under a radar):* the strongest case is KDGX (Jackson, Mississippi), with a dense echo (up to ~800 birds km⁻³ at 300–700 m). In a vertical section the density follows the path of the lowest beams, with zeros above (no-echo gates of the higher beams) and the prior below (unobserved), so the column map shows a disc or ring around the radar. 74% of the cells below −200 birds km⁻³ are around KDGX.
- *West Texas rings and specks:* thin, stationary vertical columns near KMAF, KDYX and KSJT, present at all times, are most likely wind-turbine clutter (a cleaning issue: a static-clutter mask per radar, e.g. gates with echo in most volumes of the night and near-zero radial velocity).

**Weight tests** (`P4b_weight_tests`, 6 frames 02:00–04:30, `figures/us_east/weight_tests.png`):

| Configuration | Time | Iterations | Negative cells | Negative mass | Data misfit |
| --- | ---: | ---: | ---: | ---: | ---: |
| base (prior 0.1, smoothness 0.3, h 0) | 168 s | 80 | 10.8% | 1.6% | 0.284 |
| smoothness 1 | 176 s | 81 | 9.0% | 0.9% | 0.291 |
| prior 0.3 | 81 s | 32 | 5.5% | 0.3% | 0.294 |
| vertical smoothing (h −0.5) | 164 s | 75 | 10.6% | 1.5% | 0.284 |
| smoothness 1 + vertical | 141 s | 64 | 8.8% | 0.8% | 0.289 |

The maps barely differ; the beam imprint under KDGX stays in all of them. Prior 0.3 halves the negative cells, divides the negative mass by five and converges in less than half the iterations, for a data misfit 1% higher. The vertical weighting has no visible effect.

**Non-negativity** (`nonnegative`, now on in prod with prior 0.3 as default): after CG, projected FISTA continues from the clipped solution with the diagonal preconditioner and a step 1/λmax(D⁻¹H) from 15 power iterations, until the relative step is below 1e-3. Same 6-frame window (`figures/us_east/nonnegative_tests.png`):

| Solver | Time | Negative cells | Data misfit |
| --- | ---: | ---: | ---: |
| CG only (prior 0.3) | 82 s | 5.5% | 0.2936 |
| CG + projected FISTA (13 iterations) | 128 s | 0 | 0.2944 |

With a tolerance of 1e-4 (39 iterations) the field changes by 0.7%. The non-negative solution differs from the clipped CG solution by only 2% (mass +0.25%): the constraint makes the densities valid but does not remove the beam imprint under KDGX, which comes from the data geometry (lowest beam in a dense layer, zeros above, unobserved below), not from negative lobes.

**Open issues, 23 October 2026:**
- *KDGX echo is birds.* The depolarization ratio (Kilambi et al. 2018, from ZDR and ρHV) is −5 to −2 dB at KDGX, with 2–3% of gates below the −12 dB insect limit (1–10% at bird-dominated KDMX volumes); the high ZDR fits birds seen side-on. Kept.
- *Hole above each radar (part of the beam imprint).* At the lowest sweep the density relative to 30–50 km is 0.06 / 0.27 / 0.89 at 0–5 / 5–10 / 10–15 km (median of 59 radars with birds, 04:00), although the beam there is inside the bird layer: near-range bird echo is suppressed (clutter filtering). `c.near_range_m` = 15 km drops these gates in P3. Under KDGX the empty centre of the V below 1 km is gone and the column map no longer shows a dark disc with a ring. Prior or smoothing weights could not do this: prior 1 or smoothness 3 with h −0.5 (`figures/us_east/imprint_tests.png`) keep the V and raise the data misfit by 7–37%.
- *Wind farms and towers.* `P2b_static_clutter` builds one mask per radar from the night's lowest sweeps: gates more than 8 dB above their 20 km local median in more than half of the volumes (birds and rain move; turbines do not), grown by one gate; P3 removes them below 2° elevation. Median 0.5% of low-sweep gates flagged, up to 10% (KFDR, KMBX, KBIS, KCXX: wind-farm areas). Most west Texas specks are gone; a few remain at 12:00 (e.g. 102°W 31.4°N), probably turbines seen in only part of the night or above 2°.

## Radar imprint (23 October 2026)

**Metric** (`radar_distance_bias`): each column density divided by its regional mean (Gaussian blur, σ 150 km, over the domain), averaged in bins of distance to the nearest radar and normalized to mean 1. A map without radar imprint gives 1 in every bin; the bias index is the RMS of the log ratio. The 200–300 km bin is mostly domain edges (104°W line, northern Maine, offshore) where fewer birds are expected, so the index within 200 km is also reported.

| Setting (us_east) | Bias index | Within 200 km | Shape |
| --- | ---: | ---: | --- |
| near range 15 km, no-echo zeros to 150 km | 0.067 | — | dark centre 0–25 km (0.87–0.95), bright ring 25–50 km (up to 1.14 early), dip 75–125 km, brighter between radars |
| near range 25 km, no-echo zeros within 75 km (now default) | 0.056 | 0.036 | within ±4% up to 200 km |

**Inversion weights do not remove it.** On 6 early frames (00:30–03:00): smoothness 1 or 3, horizontal smoothing (h 0.3), advection 3 and advection 3 + smoothness 1 all give a bias index of 0.060–0.063 (current 0.060); advection 3 shrinks the ring slightly but raises the data misfit by 24% (`figures/us_east/radar_bias_tests*.png`).

**What helps is data-side:** dropping near-range gates to 25 km (the lowest sweeps are suppressed close to the radar) and keeping no-echo zeros only within 75 km (bias 0.060 → 0.046 on the window, misfit unchanged).

**What remains** (e.g. `figures/us_east/zoom_midwest_0130.png`: darker discs about 75 km across in low-density areas early in the night): the beams observe the whole column near a radar but only above about 1 km far from it, so the low layers there come from the prior. The prior itself is flat with distance (±3% within 200 km), but any mismatch between the prior profile and the near-radar data becomes a ring: with `lambda_profile` 10 the far-field profile is too flat (too much density aloft), with 1 or 0.1 the far-field shape matches the near field but the bias index rises to 0.062–0.067 (dark ring at 25–50 km, bright 100–150 km).

**Option study:** soft data weights, the separable shape penalty and a correlated (GP) prior were compared in the sandbox (`prod/sandbox/README.md`).

## Best setup so far (us_east, 23 October 2026)

MistNet cleaning with screen, no-echo zeros (kept within 75 km), 2 km rain dilation; static clutter masks and 25 km near-range cut; 4 km × 100 m, US domain with 50 km offshore; Gaussian L = 150 km prior; prior 0.3, smoothness 0.3, advection 1 with takeoff/landing relaxation; CG then projected FISTA (non-negative). Night of 21–22 October: P2b 676 s, P3 856 s, prior 445 s, solve 1263 s (40 CG + projected iterations), client peak about 16 GB. Video: `figures/us_east/column_density.mp4`.

## What is needed to go live (working list)

1. **Aberrant echoes.** Remaining: rain-adjacent echo that MistNet misses, short radial "petals" near KDMX in its first night frame (02:30), and insects. Candidates: a temporal consistency check per gate and insect screening.
2. **Cleaning cost.** This is now the largest stage by far: 93 min of the 1 h 45 min `us_east` night (about 9.3 s per volume per worker with MistNet on 4 workers). Options: run MistNet only on volumes with rain-like echo (cheap screening first), load fewer moments, or reuse MistNet predictions between neighbouring times.
3. **Full US.** The 2 km × 100 m grid does not fit the solver in 36 GB; `us_east` now solves on 4 km × 200 m blocks (see below for the first run). Next: time chunks with warm starts for longer runs.
4. **Solver convergence.** A better preconditioner, warm starts between chunks or nights, and fewer time steps.
5. **Terrain.** The prior uses sea-level layers; a terrain-following prior is needed west of the Mississippi.
6. **Output.** Animated 3D view (column density plus vertical structure), and tiles.
