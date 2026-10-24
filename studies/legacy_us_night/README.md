# US Night Reconstruction

Purpose: reconstruct one full night of bird migration over the United States east of the Rocky Mountains from real NEXRAD data, from download to maps. This folder is separate from the manuscript evidence blocks (`01_`–`11_`). It reuses the shared code in `functions/` and keeps its own helpers, data, and figures here.

## Run

Set `night_date` (the local evening date) at the top of each script, then run from the project root in MATLAB:

```matlab
project_root = setup_paths();
run(fullfile(project_root, "us_night", "scripts", "UN_A_download_scans.m"))
run(fullfile(project_root, "us_night", "scripts", "UN_B_build_observations.m"))
run(fullfile(project_root, "us_night", "scripts", "UN_C_run_inversion.m"))
run(fullfile(project_root, "us_night", "scripts", "UN_D_plot_night.m"))
```

| Script | Does | Output in `data/<yyyyMMdd>/` |
| --- | --- | --- |
| `UN_A_download_scans.m` | Selects radars in `[-104 -66 24 50]` (107 stations), lists files over HTTPS, picks the volume within 10 min of each requested time, and downloads them in parallel with `curl` into the local raw cache | `manifest.mat` (raw files in `nexrad_raw_dir()`) |
| `UN_B_build_observations.m` | Builds the shared 2 km × 100 m grid, then for each radar: cleans volumes with `read_rd`, estimates velocity profiles with `load_F_elev`, builds that radar's forward-operator rows with `construct_F`, and caches the compact rows | `radar_cache/<station>.mat`, `case.mat` |
| `UN_C_run_inversion.m` | Fits the radar density prior, builds advection mappings, and solves the joint inversion on the active domain | `inversion.mat` |
| `UN_D_plot_night.m` | Maps of vertically integrated density: a nine-time overview PNG and an MP4 of every time step | `figures/<yyyyMMdd>/` |

## Default case

- Night of 21–22 October 2026, every 30 min from 23:30 to 12:00 UTC (26 times).
- 107 radars east of −104° longitude. The western limit leaves out the Front Range and mountain stations. KFDX and KUDX are the westernmost stations kept.
- Grid up to 3500 m above sea level (35 layers of 100 m). It spans all radars plus 160 km of padding, about 1750 × 1490 × 35 cells. Only cells between the ground and 3 km above it are solved.

## Terrain, radar altitude and beam blockage

Checked on 22 October 2026 with KDMX (antenna 334 m), KFDX (1432 m), KAMA and KLBB.

- **Heights are above sea level.** Kernel z is height above the antenna (cell 1 = 0–100 m). `construct_F` shifts each kernel by the radar's grid layer, so `g.z` is above sea level and `g.dem` (GTOPO30, about 1 km) gives the ground.
- **Grid top (shared fix).** `create_grid` used to take `g.z` from the kernel axis, capping the grid at 2950 m ASL whatever `heightmax` was. It now extends `g.z` with the same spacing when `heightmax` is higher. Workflows with `heightmax` ≤ 2950 m are unchanged.
- **Antenna height (shared fix).** The Level II header height is the site ground elevation, but beam heights are relative to the antenna, which is 15–112 m higher. `read_rd` now sets `scan.height` to the antenna height from `radar.csv` (ground plus tower), the same height the beam-blockage maps use. `construct_F` shifts kernels by `round(height / dz)` layers, so the placement error is now within ±50 m. Before, it was 0 to +100 m from the rounding, plus the missing tower height.
- **Velocity profiles (shared fix).** `epvvp` bins velocity by height above the antenna, but the grid layers are above sea level. `radar_uv_on_grid` shifts each profile by the antenna height; `build_radar_advection_steps`, `build_radar_velocity_prior` and `build_velocity_diagnostic_slice` all use it. For synthetic scans (height 0) it is the identity.
- **Active domain.** For each column within 160 km of a radar, the solved cells run from the ground (cell top above the DEM) to 3 km above ground, capped by the 3500 m grid top. Ground within coverage has a median of 272 m, a 99.9th percentile of 1773 m and a maximum of 2156 m (Black Hills); 0.5% of the area is above 1500 m.
- **Partial beams (shared fix).** Kernel rows sum to 1 when the whole beam is inside the grid. A row cut by the grid top, the grid edges or the active domain pushes its full value into the cells that remain; this inflated the top layers about twofold in a high-Plains test. `prepare_radar_observations` (grid) and `prepare_scans_advection_active` (active domain) now keep only rows with at least 90% of their beam weight inside (`min_row_mass = 0.9`). On the high Plains this keeps about 40% of rows.
- **Beam blockage.** `read_rd` removes gates whose lower half-power edge is not clear of the terrain horizon (see `studies/beam_blockage`). Maps exist for all 107 radars. In the test, no forward-operator weight fell below ground. Partial blockage is not corrected: the gallery shows radial streaks of lower density, for example at KFFC.
- **Density prior.** Radar profiles are fitted on sea-level layers and interpolated horizontally (IDW), so they do not follow the terrain. This is acceptable east of the Plains but not ideal across the 0–1400 m rise to the west.

## Design choices and assumptions

- **No AWS CLI.** The installed `/usr/local/bin/aws` is an Intel binary that no longer runs on this Mac. Listing goes through the public `unidata-nexrad-level2` S3 endpoint over HTTPS (shared `functions/radar/list_nexrad_day.m`, also used by `load_rd`). Downloads use `curl`.
- **Raw cache outside OneDrive.** Raw volumes go to the folder returned by `nexrad_raw_dir()`. That folder is set in the gitignored `local_settings.json` at the project root, for example `{"nexrad_raw_dir": "/Users/rafnuss/Data/nexrad_raw"}`. Without the file, the project folder `data/raw_scans/nexrad` is used. The layout is `<RADAR>_<yyyymmdd>/<file>`, shared with `load_rd`, which also searches the old project folder.
- **Radar by radar.** Full volumes take about 25 MB each in memory. To stay within memory, `UN_B` holds only one radar's volumes at a time and keeps just its compact observation rows. The rows are stacked across radars, and `observations.row_radar_id` records which radar each row comes from. The density prior uses this field in place of the global row layout of `construct_F`. Each radar is cached under `radar_cache/`, so an interrupted run resumes where it stopped. The cache is invalidated when the grid bounds or `heightmax` change.
- **Missing volumes.** A radar is kept only if at least half of the requested times are readable. A missing time leaves that radar's observation column empty. Its velocity profile is interpolated linearly in time, and the density-prior profile is filled through the temporal smoothness term (`lambda_time_density_prior = 1`). Radars whose scan pattern has no precomputed kernel in `data/derived/F_elev/150km` are skipped and reported. In that case `load_F_elev` also writes `data/reference/scan_pattern_new.mat`.
- **Water.** Water columns are excluded by default (`include_water = false`), because the density prior is zero over water.
- **Prior in time blocks.** The full-grid prior returned by `build_radar_density_prior` is too large for the whole night on 35 layers, so it is built in blocks of 4 times. Temporal smoothing of the radar profiles acts within each block only.
- **Memory chunks.** `UN_C` estimates the memory per time slice and splits the night into consecutive time chunks within `memory_limit_gb` (default 22 GB of solver workspace on this 36 GB machine; advection mappings are built per chunk). Each chunk is solved jointly with advection. Advection does not link the last time of one chunk to the first time of the next, so there can be a visible step at a chunk boundary. When the whole night fits, it is solved as a single chunk.
- **Solver.** The solver is `compute_scans_advection_active`, as in `08_us_synthetic`. It has data, prior, and advection terms but no spatial-smoothness term, unlike the Midwest solver (`smoothness = 0.3`). Negative values are clipped to zero.

## Expected cost and memory (full default night)

- **Where files go.** Raw volumes go to `nexrad_raw_dir()`. Night products (radar caches, `case.mat`, `inversion.mat`) go to `local_data_dir()/us_night/data/<night>`. Both locations are set in the gitignored `local_settings.json` and sit outside OneDrive, which otherwise spends CPU uploading several GB. Figures stay in `us_night/figures/`.
- **Download.** About 2,800 files, roughly 25 GB. Set `delete_raw = true` in `UN_B` to remove raw files once all radars are cached.
- **Processing.** Whole radars run in parallel on `n_workers = 6`, each serially through `build_radar_cache`. A worker holds one radar's raw volumes, kernels and rows (about 2–3 GB). The earlier version spread one radar's volumes over 12 workers and kept only about 3 busy; it took 36–45 s per radar (about 1.2 h for 101 radars with 15 times).
- **Density prior.** The full-grid prior returned by `build_radar_density_prior` is built in blocks of 4 times, peaking at about 3 GB (double) plus a temporary copy.
- **Inversion.** Advection mappings are built per chunk. `memory_limit_gb = 16` covers the solver workspace; the active prior (about 5 GB for the night), the observation system (2–3 GB), MATLAB and the OS come on top.

## Test run (22 October 2026)

All four scripts ran end to end on a reduced case: KILX, KIND, and KLOT on the night of 20–21 October 2026, at 02:00, 03:00, and 04:00 UTC. Download took 6 s. Processing took 14–51 s per radar. PCG converged in 157 iterations (flag 0), and the maps exported correctly. The full default night has not been run yet.

High-terrain test (22 October 2026): KFDX, KAMA and KLBB, 03:00–05:00 UTC on 21 October, 35-layer grid. No forward-operator weight below ground; 85,100 of 223,000 rows kept by the beam-weight filter. PCG converged in 169 iterations (flag 0). The reconstructed profile (2.0 → 3.9 → 2.6 birds km⁻³, peaking at 1.75–2 km above ground) matches the shape of the radar profiles.
