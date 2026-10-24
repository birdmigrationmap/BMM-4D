# Gate density: bird density from radar gates near the radar

The full write-up, with figures, is in **`report.html`** (also published as a private artifact: https://claude.ai/artifact/CVqucUraSLFZXzuiAi6bwn). It covers:
- the near-range echo deficit;
- why our current cleaning loses birds near radars and at low densities;
- the homogeneous-density estimator (step 1: sampling; step 2: the radar's clutter-filter record, CFP);
- sensitivity, calibration, and results on the sandbox sweeps.

Background: `prod/NEAR_RANGE_DEFICIT.md`.

## Run

```matlab
project_root = setup_paths();
run(fullfile(project_root, "prod", "gate_density", "scripts", "GD_A_illustration.m"))   % 30 s
run(fullfile(project_root, "prod", "gate_density", "scripts", "GD_B_sensitivity.m"))    % 40 s
run(fullfile(project_root, "prod", "gate_density", "scripts", "GD_C_calibration.m"))    % 17 min; writes GD_C_calibration.mat
run(fullfile(project_root, "prod", "gate_density", "scripts", "GD_D_sandbox.m"))        % 10 min; uses the calibration
```

| Path | Content |
| --- | --- |
| `functions/gd_beam.m` | Gate and bird model (beam width, gate length, cross-section mean and spread, volumes). |
| `functions/gd_bird_sum_table.m` | Distribution of the summed contribution of N birds in a gate (Monte Carlo). |
| `functions/gd_simulate_sweep.m`, `gd_random_field.m`, `gd_observe.m` | Simulated sweeps: Poisson birds, overlapping radials, density patterns, background echo, detection limit. |
| `functions/gd_gate_bounds.m` | Each gate as an interval of its bird contribution (exact, below detection, CFP bounds). |
| `functions/gd_estimate_sweep.m` | The estimator: density around every gate from all gates within a radius (step 1, or step 2 with CFP bounds). |
| `functions/gd_window_sum.m`, `gd_window_mean.m` | The window (range x arc box, as the read_rd smoothing) and the linear window mean. |
| `functions/gd_estimate_density.m`, `gd_estimate_rings.m` | Same likelihood on independent range rings (used for calibration). |
| `functions/gd_read_rd_filters.m` | Copy of our current read_rd cleaning steps, for simulation tests. |
| `functions/gd_radar_filters.m` | Stand-ins for the radar's point and clutter filters (assumptions). |
| `functions/gd_plot_band.m` | Mean line and band of replicate curves. |
| `scripts/GD_A_illustration.m` | Why gates are spiky, what our cleaning does, recovery by range and density. |
| `scripts/GD_B_sensitivity.m` | One factor at a time: cross-section spread, volume, background, pattern, density, detection, floor, window, N_max. |
| `scripts/GD_C_calibration.m` | Profile likelihood of cross-section spread and volume factor: recovery test, then the real sweeps. |
| `scripts/GD_D_sandbox.m` | All methods on the sandbox's real sweeps: deficit per layer and radar, example maps. |
| `figures/` | Outputs (`GD_*`); `near_range_cfp.png` and `near_range_raw.png` are copied from `prod/sandbox/figures`. |

Level II moments are read with `prod/functions/read_level2_moments.m` (REF and CFP).

## Main numbers (24 October 2026)
- **Calibrated model:** cross-section log-sd 1.75, beam volume ×1.68 (rings at 15–45 km). With it, the estimate at 30–60 km matches the radar output mean: 52 against 51 birds km⁻³.
- **Near-range deficit at 300–500 m** (density / 30–60 km density, median of 28 radars):

  | Range | Output mean | Our cleaned volumes | Step 1 |
  | --- | ---: | ---: | ---: |
  | 4.5 km | 0.06 | 0.01 | 0.20 |
  | 13 km | 0.43 | 0.18 | 0.87 |
  | 21 km | 0.89 | 0.68 | 1.04 |

- **Step 2 (CFP):** it helps above about 800 m (0.78 against 0.55 at 11 km). Low down it inflates the reference, because CFP includes ground clutter at 0.5–0.9°.
