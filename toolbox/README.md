# Toolbox overview

This folder is the BMM-4D package: the reusable code of the reconstruction, grouped by module (`radar`, `forward`, `grid`, `prior`, `advect`, `solve`, `plot`, `util`). Synthetic-case tools are in `../synthetic/`.

## Canonical grid convention

Use one Cartesian grid convention across the whole project.

- Arrays are ordered as `[x y z]`.
- `g.sz = [numel(g.x) numel(g.y) numel(g.z)]`.
- `reshape(vector, g.sz)` is the standard vector-to-grid conversion.
- `sub2ind(g.sz, x_id, y_id, z_id)` is the matching linear indexing rule.
- Use transpose or `permute(...)` only in local plotting code, not in shared computational helpers.

## First-party project functions

### Forward model and inversion

- `forward/construct_F.m`: build the sparse forward operator that maps grid cells to radar observations.
- `solve/compute_scan.m`: reconstruct a single 3D density field from one observation vector with optional regularization; residual blocks are normalized by row count before weighting.
- `forward/prepare_radar_observations.m`: compact the forward operator and observations to active rows for storage and inversion, keeping only rows with at least 90% of their beam weight inside the grid.
- `prior/build_density_prior_fast.m`: fit one vertical density profile per radar and time (land only, all cells, or separate land and water profiles) and interpolate them horizontally (IDW or Gaussian weights) into the prior field, on the full grid or on active cells.
- `advect/build_radar_advection_steps.m`: interpolate radar velocity one interval at a time and retain only compact advection mappings.
- `prior/radar_uv_on_grid.m`: map a radar velocity profile (height above the antenna) to the sea-level grid layers.
- `solve/solve_advection_lean.m`: reconstruct the 4D field (data, prior, 3D smoothness, advection) with batched PCG on active cells, single or double precision; optional projected FISTA for non-negative densities and study options (row weights, correlated prior, shape term).
- `advect/build_profile_advection_steps.m`: build one advection step per time interval from vertical-profile or grid-wise speed input.
- `advect/compute_continuity_balance.m`: compute within-grid mass balance, signed source-sink residuals, and positive take-off, landing, entering, and leaving count/rate products.

### Synthetic-only helpers

- `../synthetic/`: synthetic truth generation, synthetic mosaic/radar layouts, and generator-specific transforms and diagnostics.
- `../synthetic/make_mosaic_synthetic.m`: generate truth, radar-sampled density, and VVP-footprint velocity profiles for simple or realistic mosaic designs.
- `../synthetic/evaluate_synthetic_velocity.m`: evaluate the compact uniform or smooth synthetic truth-velocity model for generation and diagnostics.
- synthetic mosaic generators should call `create_grid(...)` and `construct_F(...)` explicitly after setting their flat-surface grid fields.
- `solve/`, `forward/` and `prior/` hold the shared inversion and forward-model code used by the real-data workflows.

### Radar I/O and cleaning

- `radar/load_rd.m`: locate or download (HTTPS, Unidata S3 bucket) raw radar files for requested radar-time combinations.
- `radar/list_nexrad_day.m`: list Level II files for one station and UTC day from the public bucket.
- `radar/nexrad_raw_dir.m`: local raw-volume cache folder, set in the gitignored `local_settings.json`.
- `radar/read_rd.m`: read radar files, remove precipitation and outliers, convert reflectivity to bird density, and optionally correct azimuthal shape bias. `scan.height` is the antenna height from `radar.csv`.
- `radar/select_radars_in_bbox.m`: read radar metadata and keep sites inside a lon/lat bounding box.

### Precomputed F-elevation support

- `forward/compute_F_elev.m`: precompute one elevation-specific forward kernel on a Cartesian grid.
- `forward/load_F_elev.m`: load the precomputed kernels needed by the current scans and attach matching metadata, including `radar_scans{i}.kernel_idx`.
- `forward/generate_F_elev.m`: script used to generate missing `F_elev` kernels from saved scan-pattern metadata.

### Vertical-profile methods

- `prior/compute_vp.m`: estimate a vertical profile with Tikhonov regularization and automatic or fixed lambda.

### Grid and matrix utilities

- `grid/create_grid.m`: build the Cartesian analysis grid and attach shared masks and DEM.
- `prior/build_horizontal_idw_weights.m`: build reusable sparse radar-to-grid horizontal IDW interpolation weights with duplicate-site handling.
- `grid/LaplacianMatrix.m`: build the 3D smoothness operator used in inversion.
- `grid/LaplacianMatrix1D.m`: build the 1D Laplacian used for vertical smoothing.
- `forward/format_data.m`: extract `obs_values`, `nonzero_obs_idx`, and optional geometry metadata from sparse observation vectors.
- `plot/colormap_wct.m`: project colormap helper.
- `plot/plot_obs.m`: processed-observation plotter for 3D observation samples or elevation-specific PPI panels.
- `plot/plot_sweep.m`: low-level raw sweep plotter for azimuth-range and local Cartesian x/y diagnostics before data are packed into observation vectors.
- `plot/plot_radar_interpolation_diag.m`: shared low-level diagnostic plot for interpolated horizontal, vertical, and time slices with raw radar values overlaid.

### Wrappers and external helpers

- `radar/mistnet_polar_radar.m` and `radar/rsl2mat_wrap.m`: wrappers or adapted helpers around upstream radar tooling.
- `util/regtools/l_corner.m` and `util/regtools/l_curve_wo_fig.m`: local regularization helpers.
- `regtools/regu/`: only the vendored Regularization Tools routines required by `compute_vp.m`.

## Current Shared Naming

These are the names that first-party code now uses across folder boundaries:

- `radar_scans`: cleaned scans and scan metadata
- `elevation_kernels`: precomputed elevation-specific kernels
- `g`: Cartesian analysis grid
- `obs_matrix`: sparse observation matrix returned by `construct_F(...)`
- `kernel_idx`: scan-to-kernel index stored on each radar scan
- `nonzero_obs_idx`: linear indices of nonzero observations
- `obs_values`: extracted nonzero observation values
- `vertical_aggregation`: descriptive local name for vertical summation operators when expanded in scripts

## Shared Saved Products

- `data/derived/simulator/scans_data.mat` stores `g`, `F`, `obs_matrix`, `elevation_kernels`, `info`, and `scans`.
- `data/derived/simulator/radar_scans.mat` stores `radar_scans`.
- `data/derived/simulator/shape_fit.mat` stores scan-cleaning shape fits.
- `data/derived/simulator/cov.mat` stores covariance fits, vertical profiles, transforms, time windows, and target scan indices.
- `info.obs_grid_size` is the canonical reshape contract for one observation vector.
- `info.obs_block_offset` stores the start offset of each observation block in the stacked system.
