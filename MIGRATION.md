# Migration log

What moved, changed or was removed since tag `paper-v1`, step by step (see `RESTRUCTURE.md`). Use it to port a fix to the paper repository, or to find where a paper function lives now.

## Step 3: duplicate implementations removed

Each removed function had a kept counterpart that computes the same method. Equality was checked on the Key West regression case (three KBYX volumes, paper settings) before removal; the regression test then passes with the paper path at 1e-11 of `paper-v1`.

| Removed | Use instead | Check |
| --- | --- | --- |
| `toolbox/prior/build_radar_density_prior.m` | `build_density_prior_fast` (default `interpolation="idw"`) | Land profiles: relative difference 7e-16. Separate land and water profiles: 7e-16. |
| `toolbox/solve/compute_scans_advection.m` + `prepare_scans_advection.m` | `solve_advection_lean` with `true(g.sz)` active and `store_class="double"` | Same PCG iterates: 3e-12 after 50 iterations, 4e-6 after 99. |
| `toolbox/solve/solve_advection_fast.m` | as above | It was a copy of the PCG path of `compute_scans_advection` (difference 0). |
| `toolbox/solve/compute_scans_advection_active.m` + `prepare_scans_advection_active.m` | `solve_advection_lean` (active cells, single precision) | Used only by `studies/legacy_us_night`, which was already broken. |

**Moved into the kept functions:**
- `build_density_prior_fast` gained the paper builder's methods:
  - `split_land_water` (separate land and water profiles, the paper's Cape May case);
  - `plot_diag`;
  - finding each row's radar from `row_ids` when `row_radar_id` is missing.
- `solve_advection_lean` gained the paper solver's check on the number of advection steps (`solve_advection_lean:InvalidStepCount`).

**One difference left:** when a solve stops before converging, the paper solver (MATLAB `pcg`) returned the iterate with the lowest residual, and `solve_advection_lean` returns the last one. At 1000 iterations without convergence on Key West the two returned fields differ by 1e-3.

**Fix:** `solve_advection_lean` returned NaN when its starting field already solves the system (prior only, no data), because the first step divided 0 by 0. It now stops at once, as `pcg` did.

**Callers updated:**
- `pipeline/scripts/P4_invert.m`, `solver_path = "paper"`, now uses the lean solver with every cell active, in double precision, without advection relaxation or non-negativity.
- The unit tests and `tests/benchmarks/benchmark_joint_mosaic_vp_inversion.m` now call the kept functions.
- The regression's paper chain is split in two:
  - `regression_case_paper_v1` keeps the paper-v1 calls, used only to rebuild the reference on a checkout of `paper-v1`;
  - `regression_case_paper` runs the same method with today's functions.

## Step 2: new layout (no code change)

Files moved with `git mv`, so `git log --follow <new path>` reaches their paper history. Code is unchanged except for folder paths in `setup_paths`, scripts, tests and comments. The regression test against `paper-v1` passes with zero difference on both paths.

**Kept on purpose:**
- Generated products still go to `local_data_dir()/prod/...`, so existing caches and cases are reused.
- Log lines keep their `[prod/<case>]` tag.
- `studies/legacy_us_night` still points to its old `us_night/functions` folder; it did not run before the move either.

**Paths:**
- `setup_paths` now adds `toolbox/` (with subfolders) instead of `functions/`.
- Scripts add `pipeline/functions` (case settings, caches, logging) and their study folder themselves.
- Synthetic tools are not on the default path; tests that need them add `synthetic/`.

**Removed:** paper-figure helpers that only the paper workflows used. The paper repository keeps its copies.

- `functions/compute/build_velocity_diagnostic_slice.m`
- `functions/compute/compute_scan_calibration_grid.m`
- `functions/plot/add_figure_bottom_space.m`
- `functions/plot/axis_label.m`
- `functions/plot/get_case_clim.m`
- `functions/plot/layout_paper_pair.m`
- `functions/util/get_case_name.m`

### Old path → new path

| Old | New |
| --- | --- |
| `functions/radar/list_nexrad_day.m` | `toolbox/radar/list_nexrad_day.m` |
| `functions/radar/load_rd.m` | `toolbox/radar/load_rd.m` |
| `functions/radar/nexrad_raw_dir.m` | `toolbox/radar/nexrad_raw_dir.m` |
| `functions/radar/read_rd.m` | `toolbox/radar/read_rd.m` |
| `functions/radar/select_radars_in_bbox.m` | `toolbox/radar/select_radars_in_bbox.m` |
| `prod/functions/read_level2_moments.m` | `toolbox/radar/read_level2_moments.m` |
| `prod/functions/clean_volume.m` | `toolbox/radar/clean_volume.m` |
| `prod/functions/static_clutter_mask.m` | `toolbox/radar/static_clutter_mask.m` |
| `prod/functions/local_contrast.m` | `toolbox/radar/local_contrast.m` |
| `prod/functions/download_volumes.m` | `toolbox/radar/download_volumes.m` |
| `functions/wrslib/rsl2mat_wrap.m` | `toolbox/radar/rsl2mat_wrap.m` |
| `functions/wrslib/mistnet_polar_radar.m` | `toolbox/radar/mistnet_polar_radar.m` |
| `functions/wrslib/matconvnet_fix/vl_nnpool.m` | `toolbox/radar/matconvnet_fix/vl_nnpool.m` |
| `functions/construct_F.m` | `toolbox/forward/construct_F.m` |
| `functions/util/format_data.m` | `toolbox/forward/format_data.m` |
| `functions/F_elev/compute_F_elev.m` | `toolbox/forward/compute_F_elev.m` |
| `functions/F_elev/load_F_elev.m` | `toolbox/forward/load_F_elev.m` |
| `functions/F_elev/generate_F_elev.m` | `toolbox/forward/generate_F_elev.m` |
| `functions/compute/prepare_radar_observations.m` | `toolbox/forward/prepare_radar_observations.m` |
| `prod/functions/attach_kernels.m` | `toolbox/forward/attach_kernels.m` |
| `functions/util/create_grid.m` | `toolbox/grid/create_grid.m` |
| `prod/functions/coarsen_grid.m` | `toolbox/grid/coarsen_grid.m` |
| `functions/util/get_water_mask.m` | `toolbox/grid/get_water_mask.m` |
| `functions/util/get_conus_mask.m` | `toolbox/grid/get_conus_mask.m` |
| `prod/functions/build_domain_mask.m` | `toolbox/grid/build_domain_mask.m` |
| `prod/functions/build_active_mask.m` | `toolbox/grid/build_active_mask.m` |
| `functions/util/find_land_water_transitions.m` | `toolbox/grid/find_land_water_transitions.m` |
| `functions/util/LaplacianMatrix.m` | `toolbox/grid/LaplacianMatrix.m` |
| `functions/util/LaplacianMatrix1D.m` | `toolbox/grid/LaplacianMatrix1D.m` |
| `prod/functions/active_laplacian.m` | `toolbox/grid/active_laplacian.m` |
| `functions/compute/compute_vp.m` | `toolbox/prior/compute_vp.m` |
| `functions/compute/build_radar_density_prior.m` | `toolbox/prior/build_radar_density_prior.m` |
| `prod/functions/build_density_prior_fast.m` | `toolbox/prior/build_density_prior_fast.m` |
| `prod/functions/build_gaussian_weights.m` | `toolbox/prior/build_gaussian_weights.m` |
| `functions/util/build_horizontal_idw_weights.m` | `toolbox/prior/build_horizontal_idw_weights.m` |
| `functions/compute/build_radar_velocity_prior.m` | `toolbox/prior/build_radar_velocity_prior.m` |
| `functions/compute/radar_uv_on_grid.m` | `toolbox/prior/radar_uv_on_grid.m` |
| `functions/compute/build_profile_advection_steps.m` | `toolbox/advect/build_profile_advection_steps.m` |
| `functions/compute/build_radar_advection_steps.m` | `toolbox/advect/build_radar_advection_steps.m` |
| `functions/compute_continuity_balance.m` | `toolbox/advect/compute_continuity_balance.m` |
| `functions/compute/compute_scan.m` | `toolbox/solve/compute_scan.m` |
| `functions/compute/prepare_scans_advection.m` | `toolbox/solve/prepare_scans_advection.m` |
| `functions/compute/prepare_scans_advection_active.m` | `toolbox/solve/prepare_scans_advection_active.m` |
| `functions/compute/compute_scans_advection.m` | `toolbox/solve/compute_scans_advection.m` |
| `functions/compute/compute_scans_advection_active.m` | `toolbox/solve/compute_scans_advection_active.m` |
| `prod/functions/solve_advection_fast.m` | `toolbox/solve/solve_advection_fast.m` |
| `prod/functions/solve_advection_lean.m` | `toolbox/solve/solve_advection_lean.m` |
| `functions/plot/plot_xy.m` | `toolbox/plot/plot_xy.m` |
| `functions/plot/plot_xz.m` | `toolbox/plot/plot_xz.m` |
| `functions/plot/plot_xyt.m` | `toolbox/plot/plot_xyt.m` |
| `functions/plot/plot_volume.m` | `toolbox/plot/plot_volume.m` |
| `functions/plot/plot_sweep.m` | `toolbox/plot/plot_sweep.m` |
| `functions/plot/plot_obs.m` | `toolbox/plot/plot_obs.m` |
| `prod/functions/plot_night_map.m` | `toolbox/plot/plot_night_map.m` |
| `functions/plot/plot_satellite_density.m` | `toolbox/plot/plot_satellite_density.m` |
| `functions/plot/plot_map_scale.m` | `toolbox/plot/plot_map_scale.m` |
| `functions/plot/add_velocity_arrow.m` | `toolbox/plot/add_velocity_arrow.m` |
| `functions/plot/get_plot_radar_xy.m` | `toolbox/plot/get_plot_radar_xy.m` |
| `functions/plot/set_axis_ticks_minmax_zero.m` | `toolbox/plot/set_axis_ticks_minmax_zero.m` |
| `functions/util/colormap_wct.m` | `toolbox/plot/colormap_wct.m` |
| `functions/plot/create_export_figure.m` | `toolbox/plot/create_export_figure.m` |
| `functions/plot/export_fixed_figure.m` | `toolbox/plot/export_fixed_figure.m` |
| `functions/plot/create_diagnostic_figure.m` | `toolbox/plot/create_diagnostic_figure.m` |
| `functions/plot/style_paper_figure.m` | `toolbox/plot/style_paper_figure.m` |
| `functions/plot/plot_radar_interpolation_diag.m` | `toolbox/plot/plot_radar_interpolation_diag.m` |
| `prod/functions/radar_distance_bias.m` | `toolbox/plot/radar_distance_bias.m` |
| `prod/functions/sun_elevation.m` | `toolbox/util/sun_elevation.m` |
| `functions/util/getSunriseSunset.m` | `toolbox/util/getSunriseSunset.m` |
| `functions/util/format_hour_labels.m` | `toolbox/util/format_hour_labels.m` |
| `prod/functions/local_data_dir.m` | `toolbox/util/local_data_dir.m` |
| `functions/regtools/l_corner.m` | `toolbox/util/regtools/l_corner.m` |
| `functions/regtools/l_curve_wo_fig.m` | `toolbox/util/regtools/l_curve_wo_fig.m` |
| `functions/README.md` | `toolbox/README.md` |
| `functions/construct_F_compact.m` | `synthetic/construct_F_compact.m` |
| `prod/functions/prod_case.m` | `pipeline/functions/prod_case.m` |
| `prod/functions/build_radar_cache.m` | `pipeline/functions/build_radar_cache.m` |
| `prod/functions/build_radar_rows.m` | `pipeline/functions/build_radar_rows.m` |
| `prod/functions/volume_cache_dir.m` | `pipeline/functions/volume_cache_dir.m` |
| `prod/functions/stage_log.m` | `pipeline/functions/stage_log.m` |
| `prod/README.md` | `pipeline/README.md` |
| `prod/.gitignore` | `pipeline/.gitignore` |
| `prod/HANDOFF.md` | `studies/notes/HANDOFF.md` |
| `prod/MODEL_REVIEW.md` | `studies/notes/MODEL_REVIEW.md` |
| `prod/NEAR_RANGE_DEFICIT.md` | `studies/notes/NEAR_RANGE_DEFICIT.md` |
| `preprocessing/beam_blockage` | `studies/beam_blockage` |
| `prod/gate_density` | `studies/gate_density` |
| `prod/cleaning_study` | `studies/cleaning_study` |
| `prod/sandbox` | `studies/sandbox` |
| `prod/experiments/two_level` | `studies/two_level` |
| `prod/legacy_us_night` | `studies/legacy_us_night` |
| `tests/wrslib` | `tests/radar` |
| `prod/scripts` | `pipeline/scripts` |
| `functions/synthetic/advect_density_sequence.m` | `synthetic/advect_density_sequence.m` |
| `functions/synthetic/advect_sample_coordinates.m` | `synthetic/advect_sample_coordinates.m` |
| `functions/synthetic/aggregate_static_recovery_metrics.m` | `synthetic/aggregate_static_recovery_metrics.m` |
| `functions/synthetic/average_velocity_time_profile.m` | `synthetic/average_velocity_time_profile.m` |
| `functions/synthetic/build_advection_generation_grid.m` | `synthetic/build_advection_generation_grid.m` |
| `functions/synthetic/build_single_radar_scan_sequence.m` | `synthetic/build_single_radar_scan_sequence.m` |
| `functions/synthetic/build_static_range_height_bins.m` | `synthetic/build_static_range_height_bins.m` |
| `functions/synthetic/compute_static_observation_support.m` | `synthetic/compute_static_observation_support.m` |
| `functions/synthetic/evaluate_synthetic_velocity.m` | `synthetic/evaluate_synthetic_velocity.m` |
| `functions/synthetic/fit_altitude_transform.m` | `synthetic/fit_altitude_transform.m` |
| `functions/synthetic/fit_kbessel_covariance.m` | `synthetic/fit_kbessel_covariance.m` |
| `functions/synthetic/forward_observation_sequence.m` | `synthetic/forward_observation_sequence.m` |
| `functions/synthetic/make_advection_synthetic.m` | `synthetic/make_advection_synthetic.m` |
| `functions/synthetic/make_mosaic_synthetic.m` | `synthetic/make_mosaic_synthetic.m` |
| `functions/synthetic/sample_synthetic_velocity_profile.m` | `synthetic/sample_synthetic_velocity_profile.m` |
| `functions/synthetic/simulate_gaussian_truth.m` | `synthetic/simulate_gaussian_truth.m` |
| `functions/synthetic/summarize_sequence_recovery_metrics.m` | `synthetic/summarize_sequence_recovery_metrics.m` |
| `functions/synthetic/summarize_static_recovery_metrics.m` | `synthetic/summarize_static_recovery_metrics.m` |
| `functions/synthetic/summarize_static_truth_descriptors.m` | `synthetic/summarize_static_truth_descriptors.m` |
| `functions/regtools/l_corner.m` (untracked) | `toolbox/util/regtools/l_corner.m` |
| `functions/regtools/regu/` (untracked) | `toolbox/util/regtools/regu/` |
