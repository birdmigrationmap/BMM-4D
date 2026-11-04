# Migration log

What moved, changed or was removed since tag `paper-v1`, step by step. Use it to port a fix to the paper repository, or to find where a paper function lives now.

## No density cap or patch filter after the estimator

**Method change (production only).** `bmm4d.radar.clean_volume` no longer caps densities at 3000 birds km^-3 or removes echo patches smaller than 25 km^2; it only thins the gates. `bmm4d.config` keeps `c.filter.data_res` only (`max_density` and `min_patch_km2` removed), so cleaned volumes go to a new cache folder. Both filters were made for the paper's chain, where no-echo gates were 0. The gate estimator gives every gate with information a positive density, so the patch filter (echo = data > 0) could only drop valid islands between missing gates; it removed nothing on 4 mid-night volumes of 12 May 2021. The cap removed 0.13% of the birds on one of the 4 (KDMX) and nothing on the others (`studies/sandbox/README.md`, "Post-estimator filters of clean_volume"). The paper's chain (`read_rd` filter `max_density`, `removeOutliers`) is unchanged.

## Gate density estimator in production

**Method change (production only).** Near the radar a gate holds too few birds for gate-by-gate filters: the outlier filters and the moving median read the empty gates and lose most of the birds within about 20 km (`studies/gate_density/report.html`). Production now turns gates into density with the gate estimator: the density that makes all gates within ±1 km most probable together.

**New:**
- `bmm4d.radar.gate`: `beam`, `bird_sum_table`, `bounds`, `loglik`, `window_sum`, `estimate_sweep` (from `studies/gate_density/functions/gd_*`).
- `bmm4d.radar.filter.estimate_density`: the filter. It replaces `max_density`, `outlier_mad`, `outlier_continuity`, `no_echo`, `to_density` and `smooth_median`. Missing gates (rain, blocked beam) stay missing.
- `read_rd` option `density`: `"median"` (default, the paper chain) or `"estimator"`. `bmm4d.config` sets `"estimator"` and no longer sets `removeOutliers`; `smoothing` (2000 m) is the estimator's window.
- `tests/radar/gateEstimatorTest.m`.

**Changed against the study code:**
- The bird-sum table is tabulated on a log-lambda grid (lambda = density x counting volume) and interpolated, which makes a volume about 4x faster (11 s against 45 s for a 150 km volume up to 7°).
- Beyond the Monte Carlo table (N > 2000 birds per gate) the sum is taken as normal instead of flagging the density as out of range. On a real sweep to 80 km the estimate is the same as the study's to 0.03% (99% of gates). On a synthetic 1500 birds km^-3 sweep at 60-80 km it reads 1493 against 1460 for the study table of 10 000 birds.
- The table uses its own random stream (seed 1) instead of `rng(1)`, so it no longer resets the caller's random numbers.
- The estimate is computed only on the rows `range_limits` keeps (plus the window), so the shape correction is fitted on those gates; on one KSGF volume this rescales the lowest sweep by about 1%.

**Unchanged:** the paper path (`read_rd` default) and the frozen production settings of the regression reference, which do not set `density`, so `paperRegressionTest` still compares the median chain. Cleaned volumes with the new settings go to a new volume cache folder (its name hashes `read_options`).

**Study cleanup (`studies/gate_density`):** the scripts call `bmm4d.radar.gate.*`; removed `gd_beam`, `gd_bird_sum_table`, `gd_gate_bounds`, `gd_gate_loglik`, `gd_estimate_sweep`, `gd_window_sum`, the unused `gd_estimate_density`, `gd_estimate_rings`, `gd_gate_posterior`, `gd_read_rd_filters`, and with the clutter mixture `gd_radar_filters`, `gd_high_persistence`, `GD_D_sandbox.m`; `GD_I_clutter_mixture.m` became `GD_I_range_bias.m` (no mixture).

## The gate estimator returns its sd

**New output (densities unchanged).** `bmm4d.radar.gate.estimate_sweep` returns `sd`, the standard deviation of each window's density (birds km^-3) from its likelihood normalized over the density grid (uniform prior in density), and at least density × log_sd where the likelihood is narrower than the grid step. `estimate_density` stores it in each gate's `sd` field, `read_rd` copies it to `dens.sweeps.sd` (missing where a later filter removed the gate). `clean_volume` does not save it yet. Re-cleaning the 168 sandbox volumes reproduces the cached densities exactly (`studies/sandbox/SB_estimator_sd.m`).

## Repository tidying

Removed from the repository (still in the history, e.g. at `v1.0.1`): `RESTRUCTURE.md` (the restructure plan, done; Steps 2-8 below), the internal notes `studies/notes/HANDOFF.md`, `MODEL_REVIEW.md` and `NEAR_RANGE_DEFICIT.md` (the near-range findings are in the gate-density report), and the studies set aside: `studies/two_level` (two-level separable model), `studies/legacy_us_night` (the US night before the pipeline) and `studies/pipeline_review`. Paths to them in the tables below are historical.

## Version 1.0.1

Bug fixes; the paper method is unchanged.

**`bmm4d.util.getSunriseSunset`.** It wrapped the UTC hour with mod 24, so an event after 00 UTC (in the Midwest, every sunset from about late March to late September) was placed on the previous UTC day: the night of 12 May 2021 got its sunset at 01:17 UTC on 12 May instead of 13 May. The local solar time is now wrapped instead, and the UTC time may cross into the next or previous day. Callers: `examples/case_settings.m` (night cases), `studies/sandbox/SB_night_maps.m` (the volumes of the five nights with a sunset after 00 UTC, 2021-05-12, 2023-05-17, 2023-05-21, 2021-09-08, 2025-09-25, came from the evening before, which was not screened for rain; the overlap and propagation tests used them) and `studies/gate_density/scripts/GD_J_select_post2020.m` (spring scans taken 3 h after the previous evening's sunset, i.e. from the night before the Dark Ecology night; each scan was still checked for rain and birds itself). Production nights (`bmm4d.config`, sun elevation) were not affected.

**Column sums count the air part of the lowest cell.** A cell is active when its top is above the ground, and the solved density of the lowest cell is about that of the cell above (median ratio 1.01). Column sums counted whole cells, so a column jumped by a layer of birds where the ground crosses a 100 m layer boundary, and the maps showed blocky contours following the terrain: between neighbouring 4 km columns a one-layer ground step changed the column by 8% (21 October 2026) and 14% (12 May 2021, dense low layer). New `bmm4d.grid.air_fraction(g, active_id)`: the share of each cell above the ground (max(g.dem, 0)). `bmm4d.export_maps` and `bmm4d.plot.radar_distance_bias` now weight the lowest cell by it; the step drops to 1.8% and 4% (`studies/sandbox/SB_ground.m`). The inversion is unchanged: the same share in the forward operator, and radar profiles in height above the ground for the prior, were tested and changed nothing measurable, so they were not kept. Radar-distance bias values computed before this change counted whole cells (sandbox: 0.043 then, 0.046 now).

**The radar-profile prior is clipped at zero** (`bmm4d.prior.build_density_prior_fast`): a weakly smoothed profile can dip below zero, which is not a density and breaks the advection relaxation (ratios of the prior between linked cells). No change with the paper's smoothing (`lambda_profile` 10).

**Correlated prior keeps the white prior's variance** (study option). `solve_advection_lean` option `prior_correlation_m` (and `prior_correlation_vertical_m`) used the precision (I + (l/dx)^2 L_h + (l_v/dz)^2 L_v) without normalization, so a longer correlation length also made the prior much stiffer: at l = 30 / 50 / 100 km on a 4 km grid the prior standard deviation of a cell was 9.7 / 15 / 28 times smaller than with the white prior of the same weight. The precision is now multiplied by the marginal variance on an unbounded lattice (mean of 1 / (1 + a_h lambda_h + a_v lambda_v) over the frequencies), so `prior` sets how far a cell may deviate and l only how far a deviation is correlated. Production does not set these options, so production and the paper path are unchanged. Sandbox results before 2 November 2026 with these options (`studies/sandbox/README.md`, 23 October) used the old, stiffer form.

**`bmm4d.radar.read_level2_moments`** reads Level II archives stored without bzip2 records (some 2015-2016 files).

**`setup_paths`** puts the package folder on the path on parallel workers (batch jobs), which cannot install packages.

## Version 1.0.0

The code behind the submitted manuscript (tag `paper-v1`) as the MATLAB package `bmm4d`: Steps 2-8 below, with CI and a regression test against `paper-v1`. The manifest said 0.1.0 during the restructure.

## Step 8: package and build files

**Package:**
- `toolbox/resources/mpackage.json` makes `toolbox/` the MATLAB package `bmm4d` (version 0.1.0, R2024b or later).
- `setup_paths` now installs the package in place for the session (`mpminstall(..., InPlace=true, Temporary=true)`) instead of adding `toolbox/` to the path. With the manifest present, MATLAB warns when a package folder is on the path without being installed.
- The unit tests install it the same way in their class setup.
- The paper-v1 regression script no longer removes the package from the path: today's code is all in the namespace, so the paper-v1 functions can sit beside it.

**Build:** `buildfile.m` with four `buildtool` tasks:
- `check`: Code Analyzer on toolbox, examples, tests and studies; fails on errors.
- `test`: unit tests, without the Regression tag.
- `regression`: the paper-v1 test.
- `package`: `release/BMM-4D_<version>.mltbx` (not tracked).
- `buildtool` alone runs `check` and `test`.

**Checks:**
- `check`: 0 errors.
- `test`: 62 tests pass.
- `regression`: both paths pass.
- `package`: the toolbox file holds the 92 package files and the manifest.
- A temporary in-place install from an empty path makes `bmm4d.config` available, and uninstalling removes it.

**Docs:**
- README: setup, package and build.
- AGENTS.md: a repository-layout section (package, namespaces, filters, MIGRATION.md and regression rules, `buildtool`); paths and wording no longer refer to the paper's numbered workflows.

## Step 7: settings and stage functions

The pipeline scripts became package functions. `pipeline/` is gone.

| Old | New |
| --- | --- |
| `pipeline/functions/prod_case.m` (defaults and case table) | `bmm4d.config(name, ...)` (defaults, overrides, night frames, folders) and `examples/case_settings.m` (the test cases and their suffixes) |
| `P1_get_volumes.m` | `bmm4d.get_volumes(c)` |
| `P2_clean_volumes.m` | `bmm4d.clean_volumes(c)` |
| `P2b_static_clutter.m` | `bmm4d.static_clutter(c)` |
| `P3_build_case.m` | `bmm4d.build_observations(c)` |
| `P4_invert.m` (`solver_path` variable) | `bmm4d.invert(c, solver_path)` |
| `P5_plot.m` | `bmm4d.export_maps(c)` |
| `P4a_plot_prior.m`, `P4b_weight_tests.m` | `studies/pipeline_review/plot_prior.m`, `weight_tests.m` |
| `pipeline/scripts/generate_F_elev.m` | `studies/kernels/generate_F_elev.m` |
| `pipeline/functions/stage_log.m`, `volume_cache_dir.m`, `build_radar_rows.m` | `bmm4d.util.stage_log`, `bmm4d.radar.volume_cache_dir`, `bmm4d.forward.build_radar_rows` |
| `pipeline/functions/build_radar_cache.m` | `studies/legacy_us_night/functions/` (used only there) |
| `pipeline/README.md` | `examples/README.md` |

**`bmm4d.config`:**
- It takes any setting as a name-value pair. A struct setting (e.g. `inversion`) replaces only the fields it contains; a `read_options` cell replaces only the options it names.
- `radars`, `times` and `grid_bounds` are required.
- `examples/run_case.m` runs the stages for one case. A new case is one `bmm4d.config` call.
- New `tests/configTest.m`.

**Unchanged:**
- File products and locations: `local_data_dir()/prod/<case>/`, and the volume cache keyed by the cleaning settings.
- Stage log lines now say `[bmm4d/<case>]` with the new stage names.
- Figures moved from `pipeline/figures/<case>` to `figures/<case>`: setting `figure_dir`, which replaces `pipeline_figure_dir`.

**Checks:**
- `case_settings` returns the same settings as `prod_case` for nine cases (keywest, keywest_rho, midwest, midwest_coarse, midwest_none, us_east, us_east_near15_zero75, sandbox, sandbox_dr). The only difference is `figure_dir`.
- On Key West, the old scripts P1–P5 and the new stages gave the same outputs:
  - manifest and static clutter masks identical;
  - a volume cleaned again by `clean_volumes`, identical;
  - `case.mat` identical (the grid differs only in its function handles, which give the same values);
  - same active cells; the solve took 300 iterations to the same residual in both runs.
- The solved field differs by 0.07%, within the run-to-run variation of the non-negative FISTA step (0.1–0.5%, see step 1).
- All 64 tests pass.

## MatConvNet pooling fallback removed

`toolbox/matconvnet_fix/vl_nnpool.m` (a MATLAB max/average pooling that replaced MatConvNet's on Apple Silicon) is removed, with its test.

**The cause** was in wsrlib's build, not in MatConvNet.
- wsrlib's maca64 build compiled MatConvNet with `-ffast-math`, which implies `-ffinite-math-only`.
- Clang then drops the -Inf starting value of max pooling, so windows with only negative values returned 0.
- The fix is to build with `-ffast-math -fno-finite-math-only`: commit "Keep infinities in maca64 MatConvNet builds" (482e6a8) on `master` of Rafnuss/wsrlib, part of darkecology/wsrlib#3.

**Check of the rebuilt `vl_nnpool`:**
- It matches a reference max pooling exactly.
- MistNet labels and probabilities on a KBYX volume are identical to those obtained with the fallback.

**The fallback had its own error:** average pooling with padding counted the padded zeros in the mean, which MatConvNet does not. MistNet does not use that case.

**Now in BMM-4D:**
- `mistnet_polar_radar` checks once per session that `vl_nnpool` keeps negative maxima, and stops with `mistnet_polar_radar:BrokenPooling` otherwise. It no longer changes the path.
- New test `tests/radar/matconvnetPoolingTest`.

## Step 6: cleaning split into filters

`bmm4d.radar.read_rd` (790 → 270 lines) now reads the volume, then applies cleaning filters in order. Each filter is a function in `bmm4d.radar.filter` (`toolbox/+bmm4d/+radar/+filter/`), with the signature `scan = filter(scan, settings)`.
- The filters work on `scan.gates(i)`, one per sweep, which holds:
  - `data`, in dBZ until `to_density` and in birds km⁻³ after;
  - `range`, `az`, `elev` and `vert_half_bw`;
  - `no_echo`: gates without echo in the raw sweep;
  - `rain_mask`.
- A new filter is one new file plus its name in the list.

| Filter | Was in `read_rd` | Default when |
| --- | --- | --- |
| `rain_rho`, `rain_dr`, `rain_mistnet` | rain mask by `removeRain` (`mistnet_rho` and `mistnetScreen` handled in `rain_mistnet`) | `removeRain` is not `none` |
| `rain_dilation` | growing the rain mask | `rainDilation` > 0 |
| `beam_blockage` | terrain-blocked gates | always |
| `max_density`, `outlier_mad`, `outlier_continuity` | the three outlier steps | `removeOutliers` > 0 |
| `no_echo` | no-echo gates kept as zero within `noEchoRange` | always |
| `to_density` | dBZ to birds km⁻³ | always |
| `shape` | anisotropy correction (fit in `scan.shape_fit`, returned as `shape_fit`) | `shape` true |
| `smooth_median` | moving median in range and azimuth | `smoothing` > 0 |
| `range_limits` | `minRange`, `maxRange`, `maxAlt` | always |
| `thin` | `data_res` thinning | always |

**Options:**
- The new `read_rd` option `filters` gives the order. Without it, the list is built from the existing options as above, so existing calls are unchanged.
- Each volume gets `filter_log`, a table of the net number of gates each filter made missing. Filters that fill gates (`no_echo`, `smooth_median`) count negative.
- Helpers used by the filters: `remove_rain` (adds to the rain mask and removes those gates) and `finish_plot` (exports or shows a diagnostic figure).

**Checks:**
- The regression test passes, both paths unchanged.
- Old and new `read_rd` give identical density sweeps on one KBYX volume for five option sets:
  - DR rain;
  - MistNet+ρHV with dilation;
  - no-echo gates, no outliers, no shape, full resolution;
  - no rain removal, no smoothing;
  - screened MistNet with dilation.
- New `tests/radar/radarFilterOrderTest`:
  - the explicit default list equals the implicit one;
  - leaving out `smooth_median` equals `smoothing=0`.

**Changed besides moving code:**
- The `rain` diagnostic plot is drawn after all filters, from the final rain mask.
- The `shape` plot's ρHV map uses the ρHV sweep's own ranges. Before, it took the ranges of the last reflectivity sweep.
- Unused volume fields (`sw`, `dr`, `ph`, …) are removed at the end instead of after the rain mask, so `rain_dr` can run at any position.

## Step 5: bmm4d namespace

The toolbox modules became namespaces: `toolbox/<module>/` moved to `toolbox/+bmm4d/+<module>/`, and every call was rewritten to its full name. 391 calls in 85 files changed.
- The rewrite was mechanical: code only, not comments, strings, field names or function definition lines.
- Calls built inside `evalc` strings (sandbox, two-level study, a forward-operator test) and the `which` check of the regression script were updated by hand.
- Function names are unchanged, and no code changed otherwise. All 59 tests pass, including the paper-v1 regression test.

**Not in the namespace:**
- `toolbox/matconvnet_fix/vl_nnpool.m`: it must override MatConvNet's function of the same name, so it stays a plain folder on the path. `mistnet_polar_radar` looks for it there.
- `generate_F_elev.m`: it is a script, which a namespace cannot hold, so it moved to `pipeline/scripts/`.
- `tests/regression/regression_case_paper_v1.m`: it keeps the paper-v1 names, since it runs only on a checkout of `paper-v1`.

**Short names in a script:** add `import bmm4d.radar.*` (or another module) at the top.

### Old name → new name

| Old | New |
| --- | --- |
| `build_profile_advection_steps` | `bmm4d.advect.build_profile_advection_steps` |
| `build_radar_advection_steps` | `bmm4d.advect.build_radar_advection_steps` |
| `compute_continuity_balance` | `bmm4d.advect.compute_continuity_balance` |
| `attach_kernels` | `bmm4d.forward.attach_kernels` |
| `compute_F_elev` | `bmm4d.forward.compute_F_elev` |
| `construct_F` | `bmm4d.forward.construct_F` |
| `format_data` | `bmm4d.forward.format_data` |
| `load_F_elev` | `bmm4d.forward.load_F_elev` |
| `prepare_radar_observations` | `bmm4d.forward.prepare_radar_observations` |
| `LaplacianMatrix` | `bmm4d.grid.LaplacianMatrix` |
| `LaplacianMatrix1D` | `bmm4d.grid.LaplacianMatrix1D` |
| `active_laplacian` | `bmm4d.grid.active_laplacian` |
| `build_active_mask` | `bmm4d.grid.build_active_mask` |
| `build_domain_mask` | `bmm4d.grid.build_domain_mask` |
| `coarsen_grid` | `bmm4d.grid.coarsen_grid` |
| `create_grid` | `bmm4d.grid.create_grid` |
| `find_land_water_transitions` | `bmm4d.grid.find_land_water_transitions` |
| `get_conus_mask` | `bmm4d.grid.get_conus_mask` |
| `get_water_mask` | `bmm4d.grid.get_water_mask` |
| `add_velocity_arrow` | `bmm4d.plot.add_velocity_arrow` |
| `colormap_wct` | `bmm4d.plot.colormap_wct` |
| `create_diagnostic_figure` | `bmm4d.plot.create_diagnostic_figure` |
| `create_export_figure` | `bmm4d.plot.create_export_figure` |
| `export_fixed_figure` | `bmm4d.plot.export_fixed_figure` |
| `get_plot_radar_xy` | `bmm4d.plot.get_plot_radar_xy` |
| `plot_map_scale` | `bmm4d.plot.plot_map_scale` |
| `plot_night_map` | `bmm4d.plot.plot_night_map` |
| `plot_obs` | `bmm4d.plot.plot_obs` |
| `plot_radar_interpolation_diag` | `bmm4d.plot.plot_radar_interpolation_diag` |
| `plot_satellite_density` | `bmm4d.plot.plot_satellite_density` |
| `plot_sweep` | `bmm4d.plot.plot_sweep` |
| `plot_volume` | `bmm4d.plot.plot_volume` |
| `plot_xy` | `bmm4d.plot.plot_xy` |
| `plot_xyt` | `bmm4d.plot.plot_xyt` |
| `plot_xz` | `bmm4d.plot.plot_xz` |
| `radar_distance_bias` | `bmm4d.plot.radar_distance_bias` |
| `set_axis_ticks_minmax_zero` | `bmm4d.plot.set_axis_ticks_minmax_zero` |
| `style_paper_figure` | `bmm4d.plot.style_paper_figure` |
| `build_density_prior_fast` | `bmm4d.prior.build_density_prior_fast` |
| `build_gaussian_weights` | `bmm4d.prior.build_gaussian_weights` |
| `build_horizontal_idw_weights` | `bmm4d.prior.build_horizontal_idw_weights` |
| `build_radar_velocity_prior` | `bmm4d.prior.build_radar_velocity_prior` |
| `compute_vp` | `bmm4d.prior.compute_vp` |
| `l_curve_wo_fig` | `bmm4d.prior.l_curve_wo_fig` |
| `radar_uv_on_grid` | `bmm4d.prior.radar_uv_on_grid` |
| `clean_volume` | `bmm4d.radar.clean_volume` |
| `download_volumes` | `bmm4d.radar.download_volumes` |
| `list_nexrad_day` | `bmm4d.radar.list_nexrad_day` |
| `load_rd` | `bmm4d.radar.load_rd` |
| `local_contrast` | `bmm4d.radar.local_contrast` |
| `mistnet_polar_radar` | `bmm4d.radar.mistnet_polar_radar` |
| `nexrad_raw_dir` | `bmm4d.radar.nexrad_raw_dir` |
| `read_level2_moments` | `bmm4d.radar.read_level2_moments` |
| `read_rd` | `bmm4d.radar.read_rd` |
| `rsl2mat_wrap` | `bmm4d.radar.rsl2mat_wrap` |
| `select_radars_in_bbox` | `bmm4d.radar.select_radars_in_bbox` |
| `static_clutter_mask` | `bmm4d.radar.static_clutter_mask` |
| `compute_scan` | `bmm4d.solve.compute_scan` |
| `solve_advection_lean` | `bmm4d.solve.solve_advection_lean` |
| `data_path` | `bmm4d.util.data_path` |
| `difference_operator` | `bmm4d.util.difference_operator` |
| `format_hour_labels` | `bmm4d.util.format_hour_labels` |
| `getSunriseSunset` | `bmm4d.util.getSunriseSunset` |
| `local_data_dir` | `bmm4d.util.local_data_dir` |
| `local_folder` | `bmm4d.util.local_folder` |
| `sun_elevation` | `bmm4d.util.sun_elevation` |

## Step 4: data paths and external code

**Data folder:**
- New `data_path(...)` (`toolbox/util`) gives a path inside the data folder: `data_dir` in `local_settings.json`, or `data/` in the repository root.
- Every read of `data/...` relative to the current folder now goes through it, so the code works from any current folder. The changed files:
  - `create_grid`, `get_water_mask`, `read_rd`, `load_rd`, `load_F_elev`, `generate_F_elev`;
  - `prod_case`, P3;
  - the beam-blockage script, cleaning study scripts A and D, the legacy US-night scripts;
  - the regression test and its reference script.
- The generated files are unchanged.

**Settings:**
- New `local_folder(name, default)` reads one folder from `local_settings.json` in the repository root, whatever the current folder; relative values are taken from the repository root.
- `nexrad_raw_dir()` now uses it. Before, it read `local_settings.json` from the current folder.
- `local_data_dir()` now uses it and takes no argument. Before, it was `local_data_dir(project_root)`; callers are updated.
- `prod_case` puts figures in `pipeline_figure_dir` (default `pipeline/figures` in the repository root). Before, it used the current folder.
- `local_settings.example.json` lists every setting.

**External code**, now from `local_settings.json` instead of fixed paths in `setup_paths`:
- `wsrlib_dir` (required; clear error when missing);
- `regtools_dir` (Hansen's Regularization Tools);
- `fast_gaussian_simulation_dir`.

All three default to folders under `external/`, which is not tracked.
- The untracked Regularization Tools files (`regu/`, `l_corner.m`) moved from `toolbox/util/regtools/` to `external/regtools/`.
- Our wrapper `l_curve_wo_fig.m` moved to `toolbox/prior/`.
- `setup_paths` adds the package last, so the MatConvNet pooling fix still overrides wsrlib's copy.

**`get_l` replaced:**
- `difference_operator(n, order)` (`toolbox/util`) replaces Regularization Tools' `get_l` in `build_density_prior_fast` and `compute_vp`.
- Identical matrices, checked for n = 1 to 30 and orders 0 to 3.
- Regularization Tools is still needed by `compute_vp` (`cgsvd`, `tikhonov`, L-curve corner), so it stays an external dependency.

**Tests:**
- `load_rd` gained a `search_dirs` option. `radarFileSelectionTest` uses it to search its temporary folder instead of relying on the current folder.
- `laplacianAndGridTest` adds the Regularization Tools folder for `compute_vp` and skips when it is missing.
- `paperRegressionTest` adds the repository root itself.

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
