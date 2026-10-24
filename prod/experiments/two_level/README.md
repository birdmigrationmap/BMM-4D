# Two-level separable model (set aside, 24 October 2026)

An alternative to the production solver, built 23–24 October 2026 to remove the radar imprint. It is kept here, with its code and results, but is not part of the production pipeline. Work went back to the production (lean) model: make it fast, then improve its observation model and prior (see `prod/HANDOFF.md`).

**The model**
- Density = column amount × vertical shape, in height above a smoothed ground (30 km), per 4 km column and frame.
- Vertical shape: 3 EOF patterns of the night's radar profiles (`scripts/SB_double_layer.m`, saved as `<data>/prod/us_east/shape_patterns.mat`), with coefficients on a 40 km grid, started from the radars' regional shape (`functions/regional_shape.m`).
- Column amount: a smooth large-scale level V_L on the 40 km grid (correlated over 600 km, estimated from all radars, no radar-profile prior) plus a fine deviation on the 4 km columns (correlated over 30 km).
- Observation model instead of data cuts: the prediction near the radar is multiplied by a recovery curve a(r) = 1 − exp(−(r/r0)^k), one per elevation group, fitted on held-out rows within 50 km. One-sided no-echo and multiplicative error are implemented but off (no gain).
- Advection on the fine deviation only; the level's own change acts as takeoff and landing.
- Solved in 3 alternating rounds (amounts by PCG then FISTA with a penalty on negative density, then shapes per frame). 8.9M unknowns for the us_east night against 238M for production.

**Verdict**
- Sandbox (Midwest, 6 takeoff frames): radar-distance bias 0.043 against 0.056 for production with logistic weights; in 7-fold radar cross-validation it predicts withheld radars better (trusted L2 0.596 against 0.623, 19 of 28 radars) and without the 20% over-prediction of the current solver.
- Full us_east night: worse than production (bias 0.098 against 0.036, trusted misfit 0.395 against 0.271, 37 min and 12 GB against 13 min and 7–8 GB). The near-range curve fitted on the whole night under-corrects the rows close to each radar: a dark spot at every radar (0–25 km bin 0.79), worst at landing.
- Reusable lessons: near-range suppression depends strongly on elevation (deficit to 30–40 km at 2° and above); 3 shape patterns are enough (double layers aloft are under 1% of radar-times); a weak prior on the amount lets data imprints through, a strong one dims the map.

## Layout

| Path | Content |
| --- | --- |
| `functions/solve_two_level.m` | The two-level solver (level, deviation, shape, observation model). |
| `functions/regional_shape.m` | Radar profiles' shape coefficients interpolated to the columns. |
| `functions/solve_separable.m` | First prototype (step 2): pattern amounts on the 4 km columns, no coarse level. |
| `functions/two_level_run.m` | `sandbox_run` with the separable and two-level solvers (`solver="separable"` / `"two_level"`). |
| `scripts/SB_double_layer.m` | Shape EOFs, misplaced share for K = 1–6, double-layer count; writes `shape_patterns.mat`. |
| `scripts/SB_separable.m`, `SB_two_level.m` | Steps 2 and 3 in the sandbox. |
| `scripts/SB_observation_model.m` | Step 4: near-range curve, one-sided no-echo, multiplicative error. |
| `scripts/SB_validation.m` | 7-fold radar cross-validation (uses `prod/sandbox/sandbox_fold.m`); predictions cached in `<data>/prod/sandbox/validation/`. |
| `scripts/SB_far_range.m` | Far-range deficit diagnosis and fix. |
| `scripts/P4c_invert_two_level.m` | Full-night inversion on the uncut case `us_east_near0_zero150`; saves `inversion_two_level.mat` there. |
| `scripts/P5c_compare_two_level.m` | Two-level against the production run (maps, radar bias, totals, misfit). |
| `figures/sandbox/`, `figures/us_east/` | Outputs of the scripts above. |

The scripts add `prod/functions`, `prod/sandbox` and `experiments/two_level/functions` to the path. Data stay in `local_data_dir()/prod/`: `us_east/shape_patterns.mat`, the uncut case `us_east_near0_zero150/` (case.mat 1.5 GB, inversion_two_level.mat 0.85 GB) and `sandbox/validation/`.

## Design and decisions (from MODEL_REVIEW, 23 October 2026)

### 5. Proposed model

#### Field
Density = column amount × vertical shape, in height above ground:

ρ(h, z, t) = V(h, t) · Σ_k c_k(h, t) φ_k(z − ground(h)),

- **φ_1…φ_3:** the mean shape and the first two shape patterns (EOFs of the radar profiles above ground, recomputed per night or per season), with c_1 = 1 and c_2, c_3 small.
- **V:** the column amount (birds km⁻²).

Linear form for the solver: x_k = V·c_k, so ρ = Σ_k x_k φ_k. That gives 3 2D fields per frame instead of 30 layers (about 10× fewer unknowns). The shape terms are regularized relative to x_1 (|x_k| ≲ a few tenths of x_1).

#### Prior structure, all relative to the local level (heteroscedastic)
- **Large-scale level V_L:** on a coarse grid (about 40 km), log-scale smooth, correlated over about 500–800 km and about 2–3 h. Starting values from section 3 and BMM-US October.
  - Its mean can be the night-time template: a sun-angle curve for take-off and landing.
  - This replaces the hand-built radar prior. It's estimated from all radars at once, so it has no seams or bullseyes.
- **Fine-scale deviation:** V = V_L · (1 + ε), with ε correlated over about 20–50 km. Its variance (about 0.4² from the neighbour statistics) is relative, so where V_L is near zero the fine field is near zero too.
- **Shape fields c_2, c_3:** correlated over about 500–1000 km (the shape semivariogram), with a radar-level nugget.
- **Density weighting.** All structure terms are weighted by the local level (or a threshold around 10–20 birds km⁻²), so empty areas contribute nothing to the structure. The same weighting goes into the scores.

#### Observation model (replaces hard cuts)
- **Multiplicative error.** The observation error is proportional to the predicted value: σ = 0.3–0.5 × predicted, plus a small floor. Solved by iterative reweighting (re-solve with weights from the previous solution, 2–3 rounds), so the operator stays linear.
- **No-echo gates as inequalities.** "Predicted ≤ the detection limit at that range" (−19.5 dBZe at 50 km clear-air, +20 log10(r/50 km)), as a one-sided penalty active only when the prediction exceeds the limit. This removes the false-zero problem at long range without cutting data.
- **Near-range suppression as part of the forward model.** Predicted = a(r) × beam average, with a(r) a smooth recovery curve (e.g. 1 − exp(−(r/r0)^k)), r0 about 10–15 km, estimated from the data (globally, then possibly per radar). The near-range gates are kept, but no longer read as true densities.
- **Station calibration.** One multiplicative factor per radar (log-scale prior sd about 0.25, i.e. about ±1 dB), constrained to average 1. This matches the VID and shape nuggets and the Nussbaumer 2019 cross-validation.

#### Dynamics
- **Advection on the 3D field reconstructed from the shapes.**
  - It's already a weak constraint, with relaxation where the prior changes (take-off and landing). That's the source term W of Nussbaumer 2021.
  - Next step: an explicit W with a night-time structure (take-off after dusk, landing before dawn) instead of the prior-based relaxation.

#### Solving
Still one linear system per outer iteration: the same preconditioned conjugate gradient, non-negativity by projection, and an outer loop for the reweighting, the one-sided penalties and the level-dependent variances (2–4 rounds). The coarse level V_L and the shape coefficients are small; the main cost is the fine field at 4 km.

### 6. Decisions (23 October 2026, with Raphaël)
1. **Height reference: hybrid.** Shapes are defined above a *smoothed* ground surface (DEM low-pass filtered over about 20–50 km), where profiles are known to be more uniform. The full-resolution DEM still defines which cells are below ground (active mask) and the beam blockage in the forward model. To be tested in the prototype, since it costs little.
2. **Shape patterns: few by default, room for rare structure.** 3 patterns by default, plus up to 6 with a penalty growing with pattern order, so double layers (Dokter et al. 2011: about 0.8 and 1.9 km) are possible when the data demand them. A small free 3D residual is the fallback if rare profiles are still not representable. To check first: how often double layers occur in the radar profiles and how well K patterns represent them.
3. **Density weighting:** whatever is faster and smooth (continuous weighting by level).
4. **Near-range:** a recovery curve a(r) multiplying the prediction near the radar, estimated from held-out near rows. *Updated 23 October 2026:* one curve per elevation group (< 1.2, 1.2–2, 2–3, ≥ 3°), since the deficit at 2° and above extends to 30–40 km (section 7b). This is the default in `solve_two_level`.
5. **Detection limit:** per scan mode if the volume metadata (VCP, pulse) is readily available, otherwise the clear-air value.
6. **Station calibration factors: no.** All stations are the same WSR-88D system; revisit only if the residuals show persistent per-radar offsets.
7. **Time: advection only,** as a strong constraint. There is no radar-profile prior in the new model; the large-scale level is part of the GP estimated from all radars. No separate temporal correlation term for now.
8. **Hyperparameters: fixed** from BMM-US (October: about 600 km, a few hours) and section 3, tuned later on the multi-night test set (section 8).

### 6b. Original list of choices
1. **Height reference:** above ground (literature default, and simpler for the shape), above sea level (current), or a mix.
2. **Number of shape patterns:** 2 (91%) or 3 (96%); fixed per night or per season.
3. **Density weighting:** a threshold (BMM-US: 10 birds km⁻²) or continuous weighting by level.
4. **Near-range model:** one global recovery curve or per radar; estimated or fixed.
5. **Detection limit:** clear-air value for all scans, or per scan mode (VCP and pulse length are in the files).
6. **Station calibration factors:** yes or no for a first version.
7. **Time:** advection only (current), or also a temporal correlation in the prior (BMM-US style).
8. **Hyperparameters:** fixed from section 3 and BMM-US, or fitted per night (empirical covariances weighted by density, as in BMM-US).

### 7. Prototype plan (sandbox)
1. **Shape basis.** EOFs above ground from the sandbox radars. Check the reconstruction error of the radar profiles with 2 and 3 patterns.
2. **Separable solver.** ρ = Σ x_k φ_k, the current data term and advection, a fine Laplacian on x_1 weighted by the level, regional smoothness on x_2 and x_3. Compare with the current solver on radar bias, trusted misfit, detail and time.
3. **Two-level, heteroscedastic.** Add the coarse V_L and make the fine-scale weights relative to it.
4. **Observation model.** Near-range recovery curve, no-echo inequalities, multiplicative error by reweighting, station factors. Each one switched on separately, so its effect is measurable.
5. **Validation.** Leave-one-radar-out on the sandbox (predict a withheld radar's observations), standardized errors, and the radar-distance bias weighted by density. Compare with the current model and a Kranstauber-style column × fixed-shape baseline.

### 7b. Prototype results (23 October 2026)
Details and tables are in the sandbox sections below.

**Shape patterns, us_east night**
- 1901 radar-times, shifted above the ground smoothed over 30 km.
- 1/2/3/4/6 uncentered patterns explain 85/96/99.1/99.8/100% of the shape variance.
- With K = 3 the median misplaced share is 3% (90%: 8%, 99%: 14%), and only 0.1% of profiles have more than 20% misplaced.
- Double layers aloft are rare: 7–18 radar-times out of 1901 (0.4–0.9%). Most profiles flagged as double are a ground layer plus the migrating layer, at dusk and dawn.
- K = 4 halves the misplaced share on double layers (5% → 2.8%).
- The coefficient spread halves with each order (sd of x_k/x_1: 0.52, 0.28, 0.14, 0.065, 0.03).

**Recommendation from the profiles**
- K = 3 by default, with room up to 6.
- A Gaussian penalty on c_k with these sd, i.e. a weight growing about 4× per order.
- No free 3D residual is needed to represent the radar profiles.

**Separable solver, first sandbox comparison** (column × shape, local corrections, advection on 3D, density-weighted Laplacian)
- Without an amount anchor it fails: overshoot 75–150 km from radars, collapse beyond 150 km, poor convergence.
  - Far from radars only the upper profile is seen, so the column depends on the assumed shape.
  - A regional shape from the radar profiles reduces this, but the level must be constrained. The large-scale level V_L of step 3 is therefore needed, not optional.
- With a weak level prior (a stand-in for V_L) the imprint is as low as with the stiffest current options: bias 0.039 against 0.056 (current prior 0.3) and 0.042 (prior 0.5).
- The fit is worse: trusted misfit 0.30 at K = 3 and 0.28 at K = 6, against 0.23–0.25.
  - The extra misfit sits at high elevations near the radars.
  - The current field needs about 5–6 patterns to keep its fit, so it carries vertical detail beyond what the radar profiles show. Whether that detail is real or overfit can only be decided by leave-one-radar-out validation (step 5).
- Not faster yet: K ≥ 3 takes 2–4 min, because the shape corrections need a coarse grid.

**Implications for the plan**
- Bring V_L (step 3) forward and use it as the amount anchor.
- Put the shape coefficients on the coarse grid with it.
- Run the leave-one-radar-out test before deciding K, since in-sample misfit favours the free 3D field.
- A small free 3D residual is the fallback if the validation shows the extra vertical detail is real.

**Two-level solver (step 3)** (`solve_two_level.m`)
- Setup: coarse level V_L and coarse shape coefficients on a 40 km grid, a fine deviation relative to V_L, no radar-profile prior; the shape prior is centred on the radar profiles' regional shape.
- It works once the level GMRF is loose enough. With weight 1, V_L is flat and the deviation carries all the structure. With weight 0.03, V_L is a smooth regional field and the deviation has 0.24 rms.
- It matches the step-2 runs that used the radar prior as the anchor (bias about 0.05–0.06 at trusted misfit 0.30, K = 3), in about 60 s.
- Loosening the level brings the imprint back into V_L (bias up to 0.16), so the imprint is a data problem.

**Implication:** the observation model (step 4: near-range recovery a(r), no-echo inequalities, multiplicative error) is the next lever, together with the leave-one-radar-out test before choosing K.

**Observation model (step 4)** (`solve_two_level.m`, `SB_observation_model.m`)
- **Near-range recovery.** The deficit is elevation dependent. Low sweeps recover by about 10 km; sweeps at 2° and above are still at 0.5–0.8 between 20 and 40 km. A global a(r) fitted on held-out near rows removes little (bias 0.144 → 0.093). One curve per elevation group (r0 7.5 / 8.5 / 20.5 / 29.5 km for < 1.2 / 1.2–2 / 2–3 / ≥ 3°) reaches 0.073 with no cut and the same trusted fit; the logistic cut gives 0.061.
  - Accepted by Raphaël (decision 4 updated); now the default. The cause, perhaps cleaning or near-range sensitivity, is still unknown.
- **No-echo gates as one-sided penalties.** No effect on this night: the empirical detection limit (1% quantile of echoes, about 0.2–0.4 birds km⁻³) is far below the densities. The clear-air limit is 10× lower again. Scan-mode limits are not available in the case.
- **Multiplicative error by reweighting.** Worse in fit and bias at every floor tried. Keep it off until it can be judged by leave-one-radar-out with standardized errors.
- **Remaining imprint.** It is a far-range deficit (0.91–0.95 at 150–250 km) and a 25–100 km excess. It also appears with near rows held out, so it is the vertical extrapolation at far range, not near-range data.

**Radar cross-validation (step 5)** (`SB_validation.m`; 7 folds of 4 withheld radars, each model's own observation operator)
- **Out-of-sample errors are similar and large for all models** (relative L2 0.62–0.65 at 30–120 km). The in-sample advantage of the current solver (trusted misfit 0.226 against 0.30) does not carry over, so its extra vertical detail is mostly fitted to each radar's own data.
- **The current solver over-predicts withheld radars** (log bias −0.17 to −0.29). The two-level model is nearly unbiased (±0.06) with a lower log RMS (0.81 against 0.84–0.86), and is better on 61% of radars (largest gains at the network edge). The pooled L2 slightly favours the current solver with prior 0.5 (0.623 against 0.642).
- **Near-range recovery per elevation is confirmed:** near-range log bias −1.05 → −0.50.
- **Choices:** K = 6 adds nothing over K = 3, and multiplicative reweighting adds nothing. Decision proposal: K = 3, near-range per elevation, no reweighting.
- **Next:** the far-range deficit, and the multi-night test library (section 8) before replacing the production solver.

**Far-range deficit** (`SB_far_range.m`)
- **Two causes.** The relative (1/level²) amplitude of the fine deviation pulled V_L below dense patches (+12% deviation on observed columns). Advection of the whole density kept earlier, lower densities in the data-free areas while the level rose (−5 birds km⁻² by the last frame).
- **Fix, now the default:** an unweighted amplitude term (the smoothness part stays relative) and advection on the 3D deviation only, so the level's change acts as the source term.
- **In sample:** bias 0.073 → 0.043 at the same fit, with far bins at 0.98–1.01.
- **Cross-validated:** best model on every L2 score (trusted 0.596 against 0.623 for the current solver), better on 19 of 28 radars, with a small over-prediction (log bias −0.06).
- **Advection adds nothing measurable on the 6-frame sandbox;** to be judged on a full night.

**Full night (us_east, 26 frames, 101 radars)** (`prod/scripts/P4c_invert_two_level.m`, `P5c_compare_two_level.m`, `prod/README.md`)
- **Runs at acceptable cost:** 37 min and 12 GB peak (production 13 min, 7–8 GB). It has 8.9M unknowns instead of 238M.
- **Worse than production on this night:** bias 0.098 against 0.036, trusted misfit 0.395 against 0.271.
- **The cause is near-range under-correction.** With all near rows in, the fitted recovery curve is weaker (r0 5.5–26.5 km), columns at every radar drop (0–25 km bin 0.79), and the 25–75 km fit degrades. It is worst at landing.
- **The sandbox (6 takeoff frames) did not reveal it.** The near-range model needs more work (per frame, or a soft cut below 25 km kept on top) before the two-level model can replace production.

## Sandbox results

### Vertical shape patterns and double layers (23 October 2026)

`SB_double_layer.m` on the full us_east night (101 radars × 26 frames). The radar profiles come from `build_density_prior_fast(lambda_profile=1)`. They are shifted to height above the ground smoothed over 30 km and resampled at 50–2450 m AGL. The 12 high-ground radars, where the grid top leaves less than 2.5 km, are dropped. That leaves 1901 radar-times with VID ≥ 5, normalized to a unit column. Patterns are uncentered EOFs (the solver fits ρ = Σ x_k φ_k with free x_k). The *misplaced share* is the fraction of the column that the K-pattern fit puts at the wrong height: half the L1 distance between the two unit-column shapes. Outputs: `figures/sandbox/double_layer_eofs.png`, `figures/sandbox/double_layer_eof_summary.csv`, `figures/sandbox/double_layer_profiles.csv`, and the patterns in `<data>/prod/us_east/shape_patterns.mat`.

| K | Variance explained | Misplaced, median | 90% | 99% | VID-weighted | Share > 20% |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.852 | 0.159 | 0.412 | 0.584 | 0.210 | 38% |
| 2 | 0.958 | 0.082 | 0.169 | 0.285 | 0.086 | 5% |
| 3 | 0.991 | 0.032 | 0.079 | 0.139 | 0.041 | 0.1% |
| 4 | 0.998 | 0.016 | 0.040 | 0.071 | 0.020 | 0 |
| 5 | 0.999 | 0.008 | 0.020 | 0.034 | 0.010 | 0 |
| 6 | 1.000 | 0.005 | 0.011 | 0.019 | 0.005 | 0 |

(Centered: mean + 1/2/3 EOFs explain 78/95/99% of the anomaly variance.)

**Double layers** are two peaks of prominence ≥ 10% of the maximum, ≥ 500 m apart, the secondary ≥ 50% of the primary. They occur in 149 of 1901 radar-times (7.8%, 7.2% of the birds), on 22 radars. Most are not two layers aloft:
- 127 are a peak at the ground (50–150 m AGL) plus a layer at 1.1–2.4 km: a low layer plus the migrating layer, mostly at dusk and dawn. Examples are KTWX and KOAX at 00:30, where takeoff is in progress, and KBRO, KFDR, KGRK, KINX and KSRX through the night.
- 15 more only have the upper peak at the top of the 2.45 km window: an edge effect.
- Only 7 (0.4%) are two layers aloft. With the 150 m peaks included there are 18 (0.9%), on KDLH (4), KGRK (4), KFDR, KGRB (2 each), KARX, KDFX, KEWX, KFWS, KOAX and KTWX.

On the double layers the median misplaced share is 5.0% with K = 3, 2.8% with K = 4 and 0.7% with K = 6. The worst case at K = 3 is 15%.

**Shape coefficients.** The spread of c_k = x_k / x_1 roughly halves with each order: sd 0.52, 0.28, 0.14, 0.065, 0.030 for k = 2…6. A Gaussian prior on c_k with these sd is the penalty growing with pattern order (weight ∝ 1/sd², about 4^(k−2)). They are saved as `shape_sd`.

### Separable solver, first prototype (23 October 2026)

`solve_separable.m` models the density as ρ(h, z, t) = x_1(h, t) s(h, z, t) + Σ_{k=2..K} x_k(h, t) φ_k(z − smoothed ground).
- s is the mean shape φ_1. With `regional_patterns` J it is bent by the radar profiles' coefficients on φ_2…φ_{J+1}, interpolated with Gaussian weights weighted by VID (`regional_shape.m`).
- The full DEM still sets the active cells and beam blockage.
- The terms are:
  - the current data term;
  - advection on the 3D field Bx;
  - a fine squared Laplacian on x_1;
  - a GMRF (I + (ℓ/Δx)² L_h) on the corrections x_k with sd from `shape_sd`;
  - optionally a level prior ‖W(x_1 − level)‖², a stand-in for the step-3 V_L;
  - x_1 ≥ 0 plus a penalty on negative 3D density (projected FISTA).
- All structure terms are weighted by W = median(level)/max(level, 10 birds km⁻²), with the level taken from the radar prior's column.
- It runs from `two_level_run` with `solver="separable"`; the comparison is `SB_separable.m`.

All runs below use logistic 28/75 km data weights. A level weight of 0.12 is about the current prior 0.3 per column.

| Option | Bias | Misfit | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| current, prior 0.3 | 0.056 | 0.324 | **0.226** | 0.34 | 57 s |
| current, prior 0.5 | 0.042 | 0.351 | 0.245 | 0.30 | 36 s |
| sep K3, night mean shape, no level | 0.209 | 0.447 | 0.368 | 0.88 | 237 s |
| sep K1, regional shape, no level | 0.137 | 0.430 | 0.340 | 1.01 | 184 s |
| sep K1, regional, level 0.1 | 0.051 | 0.432 | 0.341 | 0.43 | 61 s |
| sep K1, regional, level 0.3 | 0.039 | 0.449 | 0.353 | 0.28 | 36 s |
| sep K3, regional (J = 5, 75 km), corrections 50 km, level 0.3 | **0.039** | 0.401 | 0.297 | 0.27 | 139 s |
| sep K3, same, weak corrections, level 0.12 | 0.057 | 0.366 | 0.263 | 0.34 | 118 s |
| sep K6, regional, level 0.3 | 0.040 | 0.390 | 0.278 | 0.25 | 264 s |

Projection of the current solution (prior 0.3) on K patterns, column by column, then its trusted misfit: 0.580 / 0.353 / 0.283 / 0.255 / 0.244 / 0.237 for K = 1…6 (current: 0.226). Files: `figures/sandbox/separable_projection.csv`, `figures/sandbox/separable_comparison.csv`, `figures/sandbox/separable_curves.png`, `figures/sandbox/gallery_separable.png`, plus one panel per run.

**Reading**
- **A shape is not enough without an amount anchor.**
  - With the night's mean shape, the columns overshoot 75–150 km from the radars (ratio 1.3), then collapse beyond 150 km. Far from a radar the beams only see the upper part of the profile, so the column depends on the assumed shape. The regional radar shape halves this, but a fine Laplacian alone extrapolates badly beyond the data.
  - It also converges slowly: PCG reaches 500 iterations, and stronger smoothness *raises* detail through unconverged modes.
  - A weak level prior fixes both: about 36 iterations and the expected behaviour.
- **With the level prior, the separable model removes the imprint as well as the stiffest current options** (bias 0.039–0.040, detail 0.25–0.28). It pays in fit: trusted misfit 0.30 at K = 3 against 0.245 for the current solver at the same bias (prior 0.5).
- **The misfit cost is vertical.**
  - The separable residuals are largest at 2–5° elevation within 60 km (0.55–0.73 against 0.26–0.31).
  - The current field itself needs about 5–6 night-wide patterns to keep its fit (K = 3: 0.283).
  - It is not yet known whether this vertical detail is signal or overfit. Only the leave-one-radar-out test (step 5) can tell.
  - Advection weight (0 / 0.3 / 1) and the top taper of the patterns change nothing.
- **It is not faster yet.** K = 1 runs in 36 s (same as the current solver). K ≥ 3 takes 2–4 min because the GMRF on x_k at 4 km with ℓ = 50–100 km is ill-conditioned for the Jacobi preconditioner. The x_k (or the c_k) should live on a coarse grid (about 40 km), as planned for V_L.

### Two-level solver: coarse level V_L and coarse shape (23 October 2026)

`solve_two_level.m` (step 3), run from `two_level_run(..., solver="two_level", two_level=struct(...))` and compared by `SB_two_level.m`.
- **Model:** ρ = x_1 · (φ_1 + Σ_k (P c_k) φ_k), with x_1 = P v + e.
  - v is the large-scale level V_L and c_k are the shape coefficients, both on a 40 km node grid; P is the bilinear interpolation to the 4 km columns.
  - e is the fine deviation on the columns.
- **Level term:** an intrinsic GMRF (smoothness over 600 km) with no mean. There is no radar-profile prior: the level is estimated from all radars at once.
- **Fine deviation term:** relative to the level, W = median/max(P v, 10), with a GMRF correlated over 30 km. In the first round W is uniform; later rounds use the previous level.
- **Shape term:** coefficients are penalized by (c − c0)'(I + (ℓ/Δ)² L_c)(c − c0)/sd_k². c0 is the radar profiles' regional shape (`regional_shape.m`), or 0, the night's mean shape.
- **Solving:** ρ is bilinear in (x_1, c), so each of 3 rounds solves the amount step (PCG + FISTA penalty on negative density, advection on 3D) and then the shape step (one direct solve per time, advection ignored).
- **Remaining radar dependence:** the advection relaxation still uses the radar prior.

All runs use the logistic 28/75 km data weights. "Fine" is the weight of the deviation prior; "level" is the weight of the level GMRF (default 1).

| Option | Bias | Misfit | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| current, prior 0.3 | 0.056 | 0.324 | **0.226** | 0.34 | 24 s |
| current, prior 0.5 | 0.042 | 0.351 | 0.245 | 0.30 | 18 s |
| step 2: sep K3, radar-prior level 0.3 | **0.039** | 0.401 | 0.297 | 0.27 | 61 s |
| 2L K3, mean-shape start, fine 0.1 | 0.121 | 0.420 | 0.318 | 0.34 | 204 s |
| 2L K3, fine 0.2 | 0.046 | 0.421 | 0.319 | 0.26 | 184 s |
| 2L K3, fine 0.1 | 0.066 | 0.395 | 0.296 | 0.37 | 158 s |
| 2L K6, fine 0.2 | 0.049 | 0.414 | 0.300 | 0.25 | 208 s |
| 2L K6, fine 0.1 | 0.068 | 0.385 | 0.278 | 0.37 | 179 s |
| 2L K3, fine 0.2, level 0.03 | 0.061 | 0.401 | 0.299 | 0.31 | 63 s |
| 2L K3, fine 0.2, level 0.01 | 0.087 | 0.396 | 0.294 | 0.44 | 63 s |
| 2L K6, fine 0.2, level 0.03 | 0.061 | 0.391 | 0.281 | 0.31 | 84 s |

Figures: `figures/sandbox/two_level_tradeoff.png` (bias against trusted misfit), `figures/sandbox/two_level_fields.png` (amount, V_L, relative deviation, c_2, c_3 and the domain-mean shape for "2L K3, fine 0.2, level 0.03"), `figures/sandbox/two_level_comparison.csv` and one panel per run.

**Reading**
- **The two-level structure works as intended once the level is loose enough.**
  - With the default level weight 1, V_L is almost constant (about 60 birds km⁻²) and the "deviation" carries all the large-scale structure (±100%).
  - At level 0.03, V_L is a smooth regional field (5–95%: 18–88 birds km⁻²) and the deviation has 0.24 rms, below the 0.4 coefficient of variation between neighbouring radars.
  - The PCG then converges 3× faster (63 s for K = 3).
- **The amount anchor no longer needs the radar prior.** It gives the same trade-off as the step-2 runs that used the radar-prior column as level: bias 0.046–0.061 at trusted misfit 0.30–0.32.
- **The radar shape matters most for the bias.** Starting from the night's mean shape (c0 = 0) gives bias 0.121, because the sandbox profiles sit higher than the night's mean (`figures/sandbox/two_level_fields.png`). The data alone don't pull the coarse shape there in 3 rounds. The shape prior must be centred on a regional shape.
- **The same stiffness trade-off as everywhere else.** Loosening the level lets V_L pick up the near-radar deficit (a faint bullseye at KDMX): bias 0.061 → 0.087 → 0.161 for level weight 0.03 → 0.01 → 0.003.
  - The imprint is in the data, not in how the level is regularized.
  - It has to be removed by the observation model (step 4: near-range recovery a(r), no-echo gates as inequalities), not by smoothing.
- **Fit is limited by the vertical basis.** K = 3: 0.296–0.319; K = 6: 0.278–0.300. The ceilings from projecting the current field are 0.283 and 0.237; coarse shapes at 40 km lose a bit more.
- In-sample, every separable variant is dominated by the current solver (prior 0.5: 0.042 / 0.245). Only the leave-one-radar-out test can say whether the current solver's extra vertical detail is real.

### Observation model in the two-level solver (23 October 2026)

Step 4, in `solve_two_level.m` and compared by `SB_observation_model.m`. Each part is switchable and updated between the outer rounds from the previous prediction:
- `near_range`: predicted = a(r) F ρ, with a(r) = 1 − exp(−(r/r0)^k).
  - `"fixed"` uses r0 10 km and k 3.8, from the measured lowest-sweep suppression.
  - `"estimate"` leaves the rows within `near_range_fit_max_m` out of round 1, fits r0 and k on their echo values against that held-out prediction, then uses all rows with a(r).
  - The fit gives one curve, or one per elevation group with `near_range_elev_edges`.
- `no_echo="one_sided"`: a no-echo gate means predicted ≤ L(r), applied as an active set (gates predicted below L drop out, the others are pulled to L).
  - L is `"empirical"` (the 1% quantile of the echo values by range: 0.05 / 0.22 / 0.35 / 0.37 birds km⁻³ at 25 / 50 / 100 / 150 km) or `"clear_air"` (−19.5 dBZ at 50 km + 20 log10(r/50 km) with 11 cm²: 0.02 at 50 km, about 10× lower than the weakest echoes).
  - Scan-mode limits are not used: the VCP is not in the case file.
- `error_cv`: multiplicative error, σ = cv·max(pred, 0) + floor, data weights 1/σ² normalized to median 1 on the echo rows.

Every observation-model run uses all rows (no cuts). The two-level settings are those of step 3 (K = 3, radar regional shape, fine 0.2, level 0.03).

| Option | Bias | Misfit | Trusted misfit | Detail | Time |
| --- | ---: | ---: | ---: | ---: | ---: |
| logistic weights (step 3 reference) | 0.061 | 0.401 | 0.299 | 0.31 | 67 s |
| no cuts, no model | 0.144 | 0.338 | 0.299 | 0.36 | 64 s |
| held out < 40 km (diagnostic) | 0.056 | 0.440 | 0.320 | 0.30 | 63 s |
| near fixed (10 km, k 3.8) | 0.104 | 0.368 | 0.300 | 0.34 | 80 s |
| near estimated, one curve (11 km, k 2.75) | 0.093 | 0.377 | 0.299 | 0.34 | 95 s |
| **near estimated by elevation** | **0.073** | 0.374 | 0.303 | 0.33 | 102 s |
| no-echo one-sided | 0.144 | 0.338 | 0.299 | 0.36 | 66 s |
| error cv 0.4 (floor 1) | 0.200 | 0.514 | 0.455 | 0.46 | 88 s |
| near elev + no-echo | 0.073 | 0.374 | 0.303 | 0.33 | 99 s |
| near elev + error cv 0.4 | 0.110 | 0.493 | 0.433 | 0.41 | 120 s |
| near elev + error, floor 15 | 0.078 | 0.415 | 0.341 | 0.32 | 101 s |
| all three | 0.110 | 0.493 | 0.432 | 0.40 | 117 s |
| current solver, logistic, prior 0.3 | 0.056 | 0.324 | 0.226 | 0.34 | 23 s |

The misfit columns are computed on F ρ without a(r), so near rows count as misfit even when a(r) explains them; the trusted rows (30–120 km) are hardly affected.

Held-out near-range ratio, observed/predicted on echo rows; the prediction comes from the run without rows < 40 km:

| Elevation | 0–5 | 5–10 | 10–15 | 15–20 | 20–25 | 25–30 | 30–35 | 35–40 | 40–45 km |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| < 1.2° | 0.13 | 0.48 | 1.12 | 1.33 | 1.23 | 1.18 | 1.15 | 1.09 | 1.04 |
| 1.2–2° | 0.15 | 0.41 | 0.69 | 0.75 | 0.81 | 0.90 | 0.94 | 0.96 | 0.97 |
| 2–3° | 0.08 | 0.13 | 0.27 | 0.45 | 0.64 | 0.76 | 0.83 | 0.90 | 0.93 |
| ≥ 3° | 0.02 | 0.07 | 0.22 | 0.39 | 0.52 | 0.66 | 0.75 | 0.82 | 0.75 |

Fitted curves (r0, k): 7.5 km, 4.5 below 1.2°; 8.5 km, 2.0 at 1.2–2°; 20.5 km, 1.25 at 2–3°; 29.5 km, 1.25 at ≥ 3°. Files: `figures/sandbox/observation_model.png`, `figures/sandbox/observation_model_comparison.csv`, `figures/sandbox/observation_model_near_curve.csv`.

**Reading**
- **The near-range deficit depends on elevation, which is why one curve fails.**
  - The lowest sweeps recover by about 10 km, as in the earlier measurement. Between 10 and 40 km they even read 10–30% *above* the prediction, so the shape puts too few birds in the lowest layer near radars.
  - The sweeps at 2° and above are still suppressed at 20–40 km (0.5–0.8). That explains why a 25–28 km cut worked.
  - A recovery curve per elevation group, estimated from held-out rows, removes about 70% of the imprint without any data cut (bias 0.144 → 0.073; logistic cut 0.061). It leaves the trusted fit unchanged.
  - The cause of the high-elevation suppression is not known (cleaning, MistNet, or near-range sensitivity). It is worth checking in the cleaning before modelling it further.
- **The one-sided no-echo term changes nothing on this night.** The empirical detection limit (about 0.2–0.4 birds km⁻³) is far below the densities, so "below the limit" is almost "zero". It only drops 2–30% of the no-echo gates.
- **Multiplicative reweighting makes things worse, at every floor.** Down-weighting the dense rows costs fit and brings the imprint back: trusted misfit 0.34–0.46, bias 0.08–0.20. An unweighted L2 score is not the right judge of a weighted model, but the bias worsens too. Keep it off until the leave-one-radar-out test with standardized errors.
- **What remains is not near range.** The radar-distance curve of the best run is 0.94 at 0–25 km, 1.07–1.09 at 25–100 km and 0.91–0.95 at 150–250 km. The same 150–250 km deficit is in the held-out run, so it comes from far range, where beams only see the upper profile. The current solver with the radar prior does not have it (1.04–1.07 there).

### Radar cross-validation (24 October 2026)

`SB_validation.m` (step 5).
- **Folds.** The 28 sandbox radars form 7 folds of 4. There are 4 longitude bands of 7 radars, and fold j takes the j-th radar by latitude in each band. Every radar is withheld once, giving leave-four-radars-out with 7 solves per method instead of 28.
- **What a fold removes.** `sandbox_fold.m` removes the withheld radars' observations and drops them from the radar prior and the regional shape. The advection wind still uses all radars.
- **Predictions.** Each model predicts the withheld observations with its own observation operator: the fold's near-range curve for two-level, none for the current solver.
- **Caching.** Predictions are cached in `<data>/prod/sandbox/validation/`.
- **Folds of neighbours.** Some folds withhold neighbours (fold 5: KDMX, KDVN, KIND, KOAX), so gaps can be 300+ km.

Scores on the withheld rows:
- **Relative L2:** no-echo gates count as 0; dominated by dense rows.
- **Log bias:** mean of log((obs + 1)/(pred + 1)) on echo rows, a multiplicative error; > 0 means the model predicts too little. Log RMS is its root mean square.

| Method | Near L2 | Near log bias | Near log RMS | Trusted L2 | Trusted log bias | Trusted log RMS | Far L2 | Far log bias | Far log RMS | s/fold |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current, no cuts, prior 0.3 | 0.801 | −1.05 | 1.68 | 0.630 | −0.17 | 0.858 | 0.664 | −0.24 | 0.798 | 31 |
| current, logistic, prior 0.3 | 0.801 | −1.05 | 1.68 | 0.630 | −0.18 | 0.858 | 0.664 | −0.26 | 0.797 | 31 |
| current, logistic, prior 0.5 | 0.803 | −1.05 | 1.68 | **0.623** | −0.19 | 0.843 | **0.651** | −0.29 | 0.788 | 26 |
| 2L K3, logistic, no obs model | 0.775 | −0.95 | 1.62 | 0.638 | +0.03 | 0.819 | 0.657 | −0.04 | 0.775 | 87 |
| **2L K3, near elev (default)** | **0.755** | −0.50 | 1.27 | 0.642 | +0.06 | 0.814 | 0.657 | −0.04 | 0.773 | 126 |
| 2L K6, near elev | 0.762 | −0.49 | 1.28 | 0.648 | +0.06 | **0.806** | 0.656 | −0.04 | **0.755** | 136 |
| 2L K3, near elev + error 0.4 | 0.763 | **−0.43** | **1.26** | 0.650 | +0.11 | 0.819 | 0.667 | **0.00** | 0.772 | 118 |

Per withheld radar (trusted L2): median 0.623–0.630 for the current solver and 0.593–0.612 for two-level. The two-level K3 default beats the current solver (prior 0.5) on 17 of 28 radars (61%).
- **Largest gains** are at radars on the edge of the network, or whose neighbours were withheld too: KABR 0.94 → 0.56, KVNX 0.99 → 0.70, KDMX 0.97 → 0.77.
- **Losses:** KDVN, KOAX, KSGF, KLSX and KVWX (0.42–0.65 → 0.55–0.75).

Files: `figures/sandbox/validation.png`, `figures/sandbox/validation_scores.csv`, `figures/sandbox/validation_per_radar_misfit.csv`.

**Reading**
- **Out-of-sample errors are large for every model,** about 0.62–0.65 relative L2 against 0.23–0.30 in sample. The in-sample differences in trusted misfit (0.226 against 0.30) do not carry over. The extra vertical detail of the current solver (and of K = 6 over K = 3) is mostly fitted to the radar's own data and does not predict a withheld radar better.
- **The two-level model is unbiased; the current solver over-predicts.** The current solver's log bias is −0.17 to −0.29 (predictions about 20–25% too high in multiplicative terms), likely from the radar-profile prior spread from neighbours. Two-level is within ±0.06 at trusted and far range, with a lower log RMS (0.81 against 0.84–0.86). The pooled L2, dominated by the densest rows, slightly favours the current solver with prior 0.5 (0.623 against 0.642).
- **The near-range model is confirmed out of sample.** On the withheld near rows, log bias goes from −1.05 to −0.50 and log RMS from 1.68 to 1.27. The remaining near deficit is in the lowest kilometres (0–5 km), where a(r) is steep.
- **K.** K = 6 gains little over K = 3 (log RMS 0.806 against 0.814, L2 slightly worse) at higher cost. **K = 3 is enough.**
- **Multiplicative reweighting.** No out-of-sample gain in log RMS. It moves the bias up (+0.11 trusted) and worsens L2. **Keep it off.**
- **Caveat:** one night, one region, 7 folds. The ranking between the current solver and two-level depends on the score: L2 slightly favours the current solver, multiplicative scores favour two-level. The multi-night library (section 8) is needed before a final choice.

### Far-range deficit (24 October 2026)

`SB_far_range.m`. Columns 150–250 km from the nearest radar have no data at all: no observed beam crosses them, so they are V_L plus whatever the deviation and advection carry there. Two causes:
1. **The relative amplitude of the deviation biases V_L low.** The amplitude term weights columns by 1/level². Dense patches that the 600 km level cannot follow get small weights, so V_L sits below them: the mean of x_1 − V_L is +5.5 birds km⁻² (+12%) on observed columns, but +0.3 with the weights. Far columns fall back to V_L and end up low against their regional mean.
2. **Advection of the whole density keeps earlier, lower densities in the data-free areas.** The level rises through these frames (00:30–03:00), so x_1 − V_L there goes from about 0 at the first frame to −5 birds km⁻² at the last.

Fixes, now the `solve_two_level` defaults:
- `fine_amplitude="uniform"`: the amplitude part of the deviation prior is unweighted, while the smoothness part stays relative to the level, so V_L is the plain local mean.
- `advection_on="deviation"`: only the 3D deviation B1·e is advected. The level's own change through the night acts as the source term, as in Nussbaumer 2021.

| Option | Bias | Trusted misfit | Detail | x_1 − V_L observed | x_1 − V_L data-free, last frame | Ratio 150–200 km | 200–300 km |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| before: relative amplitude, advect density | 0.073 | 0.303 | 0.33 | +5.5 | −5.2 | 0.906 | 0.948 |
| uniform weights, advect density | 0.069 | 0.313 | 0.34 | +0.3 | −4.6 | 0.965 | 0.991 |
| uniform amplitude, advect density | 0.062 | 0.300 | 0.33 | +0.7 | −5.2 | 0.950 | 0.983 |
| uniform amplitude, no advection | 0.044 | 0.296 | 0.33 | +0.5 | −0.2 | 0.973 | 1.010 |
| **uniform amplitude, advect deviation (default)** | **0.043** | 0.299 | 0.33 | +0.3 | −0.3 | 0.976 | 1.012 |

**Cross-validation of the fix** (`SB_validation.m`, same 7 folds, new rows "2L K3, far fix" and "… no advection"):

| Method | Near L2 | Near log bias | Trusted L2 | Trusted log bias | Trusted log RMS | Far L2 | Far log bias | Far log RMS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current, logistic, prior 0.5 | 0.803 | −1.05 | 0.623 | −0.19 | 0.843 | 0.651 | −0.29 | 0.788 |
| 2L K3, near elev (before the fix) | 0.755 | −0.50 | 0.642 | +0.06 | 0.814 | 0.657 | −0.04 | 0.773 |
| **2L K3, far fix (default)** | **0.724** | −0.67 | **0.596** | −0.06 | 0.799 | **0.626** | −0.12 | 0.774 |
| 2L K3, far fix, no advection | 0.727 | −0.66 | 0.598 | −0.06 | **0.796** | 0.626 | −0.12 | 0.772 |

Per withheld radar (trusted L2) the median is 0.571 for the fixed two-level model, against 0.623 for the current solver (prior 0.5); it is better on 19 of 28 radars. The largest gains are KABR (0.94 → 0.57), KVNX (0.99 → 0.82), KINX (0.88 → 0.74) and KMKX (0.59 → 0.48). The losses are in fold 5, where four neighbours are withheld together: KDVN (0.58 → 0.96), KDMX (0.97 → 1.05), KTWX (0.96 → 1.22). Files: `figures/sandbox/far_range.png`, `figures/sandbox/far_range_comparison.csv`, `figures/sandbox/validation.png`, `figures/sandbox/validation_scores.csv`.

**Reading**
- **Both causes matter.** The uniform amplitude centres V_L, and advecting the deviation stops the data-free areas from lagging the level. Together they bring the far bins to 0.98–1.01 and the bias from 0.073 to 0.043 at the same in-sample fit. That is now better than the logistic cut with the current solver (0.056).
- **Out of sample the fix is a clear gain.** It is the best model on every L2 score (trusted 0.596 against 0.623 for the current solver and 0.642 before the fix), with a small over-prediction (log bias −0.06 to −0.12, against −0.19 to −0.29 for the current solver).
- **Advection adds nothing measurable on this case,** neither in sample nor out of sample (6 frames, 2.5 h, during takeoff). It costs little to keep. It should be judged on a full night, where the time link matters most between radars and over landing.
- **What remains:** the 0–25 km bin (0.91) and a small 25–75 km excess (1.03–1.04), i.e. the near-range model at the lowest kilometres. Fold 5 shows that 300+ km gaps are still poorly filled.

Note (24 October 2026): the solver was moved to `functions/` here (used by `scripts/P4c_invert_two_level.m`). The solver was slimmed for memory, with identical sandbox scores (bias 0.043, trusted 0.299).

## Full us_east night

### Two-level solver on the full night (24 October 2026)

`P4c_invert_two_level` on `us_east_near0_zero150`, the us_east case rebuilt by P3 without data cuts (all near-range rows, no-echo gates to 150 km): 5.83M rows (production 4.72M), nnz(F) 92M. Sandbox settings: K = 3, radar regional shape, fine 0.2, level 0.03, near range per elevation group, unweighted deviation amplitude, advection of the deviation.

**Cost**
- Solve 2213 s (37 min; production lean 782 s), with 278 PCG iterations over 3 rounds (151 / 71 / 56). Rounds take 956 / 688 / 543 s.
- FISTA does not converge in rounds 2–3 (88 and 100 iterations, the cap): negative mass 1.4%.
- Unknowns: 8.9M (3585 nodes and 339 623 columns × 26 frames), against 238M for production.
- Peak MATLAB memory 12.1 GB (client RSS sampled every 15 s), against 7–8 GB for production. This needed slimming `solve_two_level`: basis built on the fly, no scaled copy of F, uint32/single advection links, single row × time arrays. It gives identical sandbox scores.
- Prior and regional shape take 342 s.

**Comparison with production** (`P5c_compare_two_level`, `figures/us_east/two_level_maps.png`, `two_level_bias.png`, `two_level_comparison.csv`)

| | Production (lean, cuts) | Two-level (no cuts) |
| --- | ---: | ---: |
| Radar-distance bias, mean over 26 frames | **0.036** | 0.098 |
| Bias at 02:00–07:00 | 0.008–0.021 | 0.077–0.091 |
| Bias at 08:00–11:00 (landing) | 0.031–0.054 | 0.107–0.156 |
| Column / regional mean at 0–25 km from a radar | 1.01 | **0.79** |
| Column / regional mean at 200–300 km | 0.89 | 1.06 |
| Trusted misfit (30–120 km), mean | **0.271** | 0.395 |
| Domain total at the peak (04:00) | 331 M | 344 M |

Misfit by elevation and range (all frames):

| Elevation | Run | 0–10 km | 10–25 km | 25–50 km | 50–75 km | 75–100 km | 100–150 km |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| < 1.2° | production | 0.87 | 0.65 | **0.24** | **0.25** | **0.29** | **0.43** |
| < 1.2° | two-level | **0.67** | **0.50** | 0.46 | 0.38 | 0.35 | 0.46 |
| 1.2–2° | production | 0.93 | 0.66 | **0.20** | **0.22** | **0.29** | 1.00 |
| 1.2–2° | two-level | **0.66** | **0.55** | 0.35 | 0.30 | 0.32 | **0.93** |
| 2–3° | production | 1.51 | 0.84 | **0.22** | **0.30** | — | — |
| 2–3° | two-level | **0.81** | **0.57** | 0.37 | 0.52 | — | — |
| ≥ 3° | production | 2.77 | 0.99 | **0.25** | **0.67** | — | — |
| ≥ 3° | two-level | **0.87** | **0.60** | 0.48 | 1.58 | — | — |

Near-range curve fitted on this night: r0 5.5 / 6.5 / 19 / 26.5 km (k 3.5 / 2.25 / 1.25 / 1.25), against 7.5 / 8.5 / 21 / 31 km in the sandbox.

**Reading**
- **The two-level model runs on the full night at acceptable cost** (37 min, 12 GB). It agrees with production on the domain total (within 4% from 04:00 to 11:00; 10–20% higher during takeoff) and on the large-scale pattern.
- **It is worse than production on this night, for one main reason: the near-range rows are under-corrected.**
  - With all near rows included, the fitted recovery curve is weaker than in the sandbox, so near-radar columns are pulled down. That shows as a dark spot at every radar (0–25 km bin 0.79) and a worse fit at 25–75 km, where those columns must also explain the farther beams.
  - The production run avoids this simply by cutting the rows within 25 km.
  - The effect grows at landing (08:00–11:00): birds descend into the lowest, most suppressed layers, and one curve for the whole night does not follow that.
- **The sandbox comparison did not show this,** because it covered 6 frames at takeoff and 28 radars.
- **Suggested next tests,** each a 40 min run:
  1. two-level on the production (cut) case with `near_range="off"`, to separate the near-range model from the rest;
  2. two-level with a near-range curve fitted per frame, or a fixed logistic weight below 25 km kept on top of a(r);
  3. a looser FISTA cap or a better non-negativity scheme, since FISTA stops at its cap.
