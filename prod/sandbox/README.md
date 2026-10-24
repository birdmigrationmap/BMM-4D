# Sandbox: fast option tests

A small case (`prod_case("sandbox")`: Midwest part of the 21–22 October 2026 night, 28 radars, 6 frames 00:30–03:00 UTC, 4 km × 100 m, built by P1–P3 with no near-range or no-echo cuts) for trying inversion options quickly (10–250 s per run) and deciding later.

```matlab
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"), fullfile(project_root, "prod", "sandbox"))
S = sandbox_load("sandbox");                       % once, ~6 s
r = sandbox_run(S, struct("name", "my test", "prior", 0.5));
summary = sandbox_summary({r}, "my_test");          % table + comparison figure (list several runs)
run(fullfile(project_root, "prod", "sandbox", "SB_fast_solver.m"))   % an example test script
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
