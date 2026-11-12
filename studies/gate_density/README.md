# Gate density: bird density from radar gates near the radar

The full write-up, with figures, is in **`report.html`** (also published as a private artifact: https://claude.ai/artifact/66mhAz6ZYBLGufLoz58cHH). The near-range deficit is measured as a range bias b(d): one vertical profile per scan, seen through each beam's vertical weighting, and one factor per distance shared by all scans (`gd_fit_range_bias`). It has two parts:
- Part I, sampling: the gate model, each gate as an interval and its likelihood, the combination of gates in a ±1 km window (three interactive panels), comparison with the radar output and its ±1 km moving mean on simulated and real sweeps, the fixed cross-section spread, and the range bias on 100 selected rain-free scans;
- Part II, the radar's clutter filter: where and how much it removes (CFP), and a zero-isodop test of whether it removes birds (a few percent at the zero isodop on 1.8–5.1° cuts, none detectable at the lowest cuts); no correction retained.

The estimator is in production as `bmm4d.radar.gate.*` and the filter `bmm4d.radar.filter.estimate_density` (`read_rd` option `density`, set to `"estimator"` by `bmm4d.config`); the scripts here call those functions (see `MIGRATION.md`).

Background and the clutter-filter findings: Parts II and III of `report.html`.

## Run

```matlab
project_root = setup_paths();
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_A_illustration.m"))   % 1 min
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_B_sensitivity.m"))    % 1 min
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_E_selected_scans.m"))  % 45 min, three spreads (+ old cleaning of 100 volumes, cached); spread sensitivity, caches the old cleaning used by GD_I
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_F_clutter_filter.m"))  % 2 min
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_G_spread_diagnostic.m"))  % 20 min (8 workers); best spread by range, height, birds per gate, clutter filter
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_H_example_maps.m"))  % 2 min; uses GD_J outputs
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_I_range_bias.m"))  % 45 min; uses GD_J outputs
run(fullfile(project_root, "studies", "gate_density", "scripts", "GD_K_cfp_birds.m"))  % 3 min; uses GD_J outputs
```

| Path | Content |
| --- | --- |
| `toolbox/+bmm4d/+radar/+gate/` | The estimator (production): gate and bird model (`beam`: two-way weight f^4 W^2, azimuth widened by antenna rotation ×1.088, cross-section log-sd fixed at 1.25), distribution of the summed bird echo by expected birds per gate (`bird_sum_table`: Monte Carlo up to 2000 birds, normal beyond), gate intervals (`bounds`), log likelihood (`loglik`), window (`window_sum`) and estimate with its uncertainty (`estimate_sweep`). |
| `functions/gd_simulate_sweep.m`, `gd_random_field.m`, `gd_observe.m` | Simulated sweeps: Poisson birds, overlapping radials, density patterns, background echo, 0.5 dB steps, detection limit. |
| `functions/gd_window_mean.m` | Linear mean of the gates in the estimator's window (moving-mean baseline). |
| `functions/gd_fit_range_bias.m` | Vertical profile per scan and one bias per distance bin, through the beam's vertical weighting (the range-bias test). |
| `functions/gd_estimate_tiles.m` | The likelihood on non-overlapping 2 km tiles (spread diagnostic). |
| `functions/gd_cut_velocity.m`, `gd_vad_fit.m` | Radial velocity on a cut's grid (from the Doppler pass for split cuts) and a VAD fit of one ring. |
| `functions/gd_profile_max.m`, `gd_crossing.m` | Parabola-refined maximum of a profile, and the crossing of a level. |
| `functions/gd_download_level2.m`, `gd_plot_band.m` | Level II download for GD_J; mean line and band of replicate curves. |
| `scripts/GD_A_illustration.m` | Why gates are spiky; radar output, moving mean and estimator on simulated sweeps, recovery by range and density. |
| `scripts/GD_B_sensitivity.m` | One factor at a time: cross-section spread, volume, background, pattern, density, detection, floor, window, N_max. |
| `scripts/GD_E_selected_scans.m` | Output mean, estimator and old filter cleaning on the 100 rain-free scans selected for the BMM-within simulator: deficit by height layer and elevation class, and the sensitivity to the spread (1 and 1.75 against 1.25). |
| `scripts/GD_F_clutter_filter.m` | Share of gates that lost ≥ 3 dB to the radar's clutter filter, by range and elevation, for beams below 4 km (sandbox CFP). |
| `scripts/GD_G_spread_diagnostic.m` | Profile likelihood of the cross-section spread by range, height, expected birds per gate and clutter-filter loss (2 km tiles, 100 selected scans). |
| `scripts/GD_H_example_maps.m` | Example sweeps (output, moving mean, estimator) for 12 random GD_J scans, 4 per density class, one PNG per scan. |
| `scripts/GD_I_range_bias.m` | Range bias b(d) of the radar output, its moving mean and the estimator on the 124 post-2020 scans of GD_J, all scans and by density class, with a bootstrap band. |
| `scripts/GD_J_select_post2020.m` | The 124 reference scans (2021–2025, 49 radars, rain-free, three nights per radar). |
| `scripts/GD_K_cfp_birds.m` | Does the clutter filter remove birds: share filtered, output and output + CFP by the birds' expected radial velocity (VAD per 1 km ring), lowest (split) and batch cuts at 4–40 km, with a 90°-turned control, one example scan (KSGF), and the output of filtered gates by CFP outside the notch. |
| `figures/` | Outputs (`GD_*`); `near_range_cfp.png` and `near_range_raw.png` are copied from `studies/sandbox/figures`. |

Level II moments are read with `bmm4d.radar.read_level2_moments` (REF, VEL and CFP).

## Main numbers (1 November 2026)
- **Range bias b(d), 124 post-2020 scans (GD_I, 49 radars)** (1 = no bias with distance, scaled at 30–60 km), scans grouped by their own mean output at 20–60 km on the lowest cut (<20: 30, 20–100: 59, >100 birds km⁻³: 35): radar output (capped) 0.11, 0.24, 0.49 at 2.5, 4.5, 9.5 km; moving mean (of the capped output) the same as the output at every distance; estimator (production, no cap) 0.34, 0.53 (0.46–0.61), 0.78, 0.89 at 15 km, within 10% of 1 from 17 km, band includes 1 from about 29 km. At 4.5 km the estimator reads 0.30 / 0.60 / 0.74 by class; below 20 birds km⁻³ it rises out to 80 km. Reading the radar's point-target flag (CFP code 1, 12–17% of gates at 3–6 km) as "at least the output" gave 0.59 at 4.5 km (earlier report versions); production has no CFP, so it is not used.
- **Spread sensitivity (100 scans, 2015–2021 selection):** spread 1 and 1.75 scale the estimate by ×0.92 and ×1.27 against 1.25, nearly the same at every range and elevation.
- **Clutter filter (sandbox CFP, beams below 4 km):** above 1°, about 60% of gates lose ≥ 3 dB at 5 km, 45% at 10 km, 20–25% at 20 km, nearly the same at every elevation.
- **Does the filter remove birds (GD_K, 124 scans, 4–40 km)?** Batch cuts (1.8–5.1°): yes, at the zero isodop. Gates where the birds' expected |v_r| < 1 m/s lose ≥ 3 dB 1.5× (4–10 km) to 2× (20–40 km) as often as at |v_r| ≥ 4 m/s, their output is 0.65 / 0.79 at 4–10 / 10–20 km, and output + CFP 0.85 / 1.04; turning the flight direction by 90° moves the pattern. Averaged over azimuth the loss is small: output 0.94 at 4–10 km, 0.98 beyond. Lowest cuts (≤ 1.3°, split, VEL from the Doppler pass): 62% of gates within 10 km lose ≥ 3 dB at every azimuth, with no clear dependence on the birds' direction (output 0.85 at |v_r| < 1 m/s within 10 km, about 1 beyond; azimuth average 0.98). Most filtered gates are not explained by birds in the notch, so the test does not support reading all of them as "birds removed".
- **CFP strength (GD_K, outside the notch, |v_r| ≥ 4 m/s):** grouped by CFP, the output of filtered gates falls from about 1 to 0.03–0.4 of the unfiltered gates of the ring as CFP grows, but this is selection (a gate shows a large CFP when its bird power is small next to the clutter). All filtered gates together read 1.07–1.11 of the unfiltered gates on the 1.8–5.1° cuts and 0.91–1.50 on the lowest cuts: no deficit outside the notch. Conditioning on CFP ≥ 3 dB gives 0.34–0.82, which is why "filtered gates read low" and "leaving filtered gates out lifts b(d)" are not evidence of bird loss.
- **Model:** beam from the radar specifications (0.95°, ×1.088 in azimuth for the rotating antenna), range weight W² of Appendix A, ±1 km window. Cross-section log-sd fixed at 1.25, the best fit where gates hold 10–30 birds (GD_G); fitted to all tiles it rises with range (0.5 at 2–5 km to 1.84 at 60–80 km) and with birds per gate, because it absorbs density variation between gates and clutter, so it is not calibrated.
- **Simulated sweeps:** the estimator recovers 92–103% of the true density at 2–60 km for 20–200 birds km⁻³; the ±1 km moving mean of the output is unbiased but has about twice its per-gate error.

## GD_L: removed power as an upper bound (3 November 2026)
`scripts/GD_L_cfp_bounds.m` (8 min, 8 workers; dev, needs `studies/sandbox/dev/cfp_error_toolbox.patch`): the GD_I pipeline with REF and CFP, four estimators (production; filtered gates as [output, output + CFP]; plus point-clutter gates as "at least the output"; the same on batch cuts only), b(d) pooled and per group of cuts (split ≤1°, split 1.3–2.5°, batch 1.8–2.5°, batch 3–4.5°; `gd_fit_range_bias` options `group`, `normalize_groups`). The bound lifts the near-range estimate but lets ground clutter in: ×1.31–1.49 on the lowest cuts at 30–60 km on every cut, and batch cuts 1.3–2.5 times the split cuts near the radar when used on batch cuts only. Not retained; see Part III of `report.html`.

## GD_M: a level per cut type and distance, bounded by the clutter filter (4 November 2026, dev)
`scripts/GD_M_cut_levels.m` (about 3 min; `functions/gd_cfp_power.m`). Each scan gets one vertical profile, 100 m layers smoothed over three layers. It is fitted to the production estimator per sweep and range bin, through the beam's vertical weighting, as density = c(type, band) × β(λ) × profile.
- **c** is a level per cut type (split or batch) and distance band (2–10, 10–20, 20–30, 30–45, 45–80 km). It is kept within [P_out/P_all, 1.1], where P_out/P_all is the share of the echo the clutter filter left. Split cuts at 45–80 km are the anchor.
- **β(λ)** is the sampling deficit, one curve for all scans, by expected birds per gate λ.

Three fits: *before* (profile only), *levels* (bounded c) and *refined* (bounded c pulled toward 1 by one sweep's worth of pseudo-observations, plus β). Disagreement between elevations is the gate-weighted RMS of log((observed + 1)/(fitted + 1)) within 45 km. Figures: `GD_M_examples_1–4.png` (22 scans ordered by filter loss: distance–height before and refined, and the profile of every sweep within 45 km for the three fits), `GD_M_summary.png`; table `GD_M_cut_levels.csv`.
- **Sampling curve.** β = 0.33 / 0.44 / 0.67 / 0.93 / 1.04 at λ = 0.03–0.1 / 0.1–0.3 / 0.3–1 / 1–3 / 3–10 birds per gate, and 1 above. Below 0.03 it hits its cap of 1.2: few data there.
- **Disagreement**, median over scans: 0.37 before, 0.25 with levels only, 0.225 refined (sparse scans 0.49 / 0.28 / 0.25, dense 0.29 / 0.24 / 0.22). It is lower in 123 of 124 scans and higher in none. Extra parameters always fit better, so this alone does not prove the corrections right.
- **Batch cuts**, median mean log ratio at 2–10 / 10–20 / 20–30 km: −0.72 / −0.27 / −0.12 before, −0.17 / −0.08 / −0.05 refined.
- **KDVN 7 May 2022** (80% loss, 22 birds km⁻³): levels only, the batch c were 0.02 / 0.09, a ×10–50 correction. Refined they are 0.25 / 0.28, with β taking the sparse-night part.
- **KAMA 3 Apr 2021** (87% loss, cuts agree): refined, the batch c go to the 1.1 cap, because β over-corrects this scan and c pushes back. Disagreement is 0.34 before and 0.36 refined.
- **Where the filter took little** (KLNX 27 Mar 2022, 4% loss), the bound keeps the batch c near 0.95–1. The batch cuts still read low at 10–45 km: something other than the filter.
- **Near-ground bins with few birds per gate** get large factors (÷ c β down to ÷ 0.1), so single bins spike after correction (KTLX 27 Mar 2022, KSGF 9 Nov 2025).
