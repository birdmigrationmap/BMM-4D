# Functions Overview

This folder contains the reusable code used by the reconstruction workflows.

## Canonical grid convention

Use one Cartesian grid convention across the whole project.

- Arrays are ordered as `[x y z]`.
- `grid.sz = [numel(grid.x) numel(grid.y) numel(grid.z)]`.
- `reshape(vector, grid.sz)` is the standard vector-to-grid conversion.
- `sub2ind(grid.sz, x_id, y_id, z_id)` is the matching linear indexing rule.
- Use transpose or `permute(...)` only in local plotting code, not in shared computational helpers.

## First-party project functions

### Forward model and inversion

- `construct_F.m`: build the sparse forward operator that maps grid cells to radar observations.
- `compute/compute_scan.m`: reconstruct a single 3D density field from one observation vector with optional regularization; residual blocks are normalized by row count before weighting.
- `compute/prepare_radar_observations.m`: compact the forward operator and observations to active rows for storage and inversion.
- `compute/build_radar_density_prior.m`: estimate radar vertical density profiles from compact observations and interpolate them into a 4D prior field.
- `compute/build_radar_advection_steps.m`: interpolate radar velocity one interval at a time and retain only compact advection mappings.
- `compute/prepare_scans_advection.m`: prepare reusable matrix-free data, smoothing, and advection components.
- `compute/compute_scans_advection.m`: reconstruct a 4D field with matrix-free PCG applied to the original least-squares normal equations.
- `compute/build_profile_advection_steps.m`: build one advection step per time interval from vertical-profile or grid-wise speed input.
- `sinksource.m`: compute within-grid mass balance, take-off, landing, entering, and leaving flux terms.

### Synthetic-only helpers

- `synthetic/`: synthetic truth generation, synthetic mosaic/radar layouts, and generator-specific transforms and diagnostics.
- `synthetic/make_mosaic_synthetic.m`: generate truth, radar-sampled density, and VVP-footprint velocity profiles for simple or realistic mosaic designs.
- `synthetic/evaluate_synthetic_velocity.m`: evaluate the compact uniform or smooth synthetic truth-velocity model for generation and diagnostics.
- synthetic mosaic generators should call `create_grid(...)` and `construct_F(...)` explicitly after setting their flat-surface grid fields.
- `compute/` should stay limited to shared inversion and forward-model code used by the real-data workflows.

### Radar I/O and cleaning

- `radar/load_rd.m`: locate or download raw radar files for requested radar-time combinations.
- `radar/read_rd.m`: read radar files, remove precipitation and outliers, convert reflectivity to bird density, and optionally correct azimuthal shape bias.
- `radar/select_radars_in_bbox.m`: read radar metadata and keep sites inside a lon/lat bounding box.

### Precomputed F-elevation support

- `F_elev/compute_F_elev.m`: precompute one elevation-specific forward kernel on a Cartesian grid.
- `F_elev/load_F_elev.m`: load the precomputed kernels needed by the current scans and attach matching metadata, including `radar_scans{i}.kernel_idx`.
- `F_elev/generate_F_elev.m`: script used to generate missing `F_elev` kernels from saved scan-pattern metadata.

### Vertical-profile methods

- `compute/compute_vp.m`: estimate a vertical profile with Tikhonov regularization and automatic or fixed lambda.
- `compute/compute_vp_thiko.m`: estimate a vertical profile with standard-form Tikhonov regularization.
- `vp/compute_vp_flex.m`: estimate a vertical profile with a direct LSQR system and simple penalties.
- `vp/compute_vp_mean.m`: compute a mean profile by altitude binning observations.
- `vp/compute_vp_ratio.m`: estimate a normalized profile from matched inter-elevation observation ratios.
- `vp/compute_vp_buler.m`: implement the Buler-style profile reconstruction based on beam-height weighting.
- `vp/compute_vp_hu.m`: estimate a vertical profile with standard-form Tikhonov regularization.

### Grid and matrix utilities

- `util/create_grid.m`: build the Cartesian analysis grid and attach shared masks and DEM.
- `util/build_horizontal_idw_weights.m`: build reusable sparse radar-to-grid horizontal IDW interpolation weights with duplicate-site handling.
- `util/LaplacianMatrix.m`: build the 3D smoothness operator used in inversion.
- `util/LaplacianMatrix1D.m`: build the 1D Laplacian used for vertical smoothing.
- `util/format_data.m`: extract `obs_values`, `nonzero_obs_idx`, and optional geometry metadata from sparse observation vectors.
- `util/lldistkm.m`: geographic distance helper.
- `util/colormap_wct.m`: project colormap helper.
- `plot/plot_obs.m`: processed-observation plotter for 3D observation samples or elevation-specific PPI panels.
- `plot/plot_sweep.m`: low-level raw sweep plotter for azimuth-range and local Cartesian x/y diagnostics before data are packed into observation vectors.
- `plot/plot_radar_interpolation_diag.m`: shared low-level diagnostic plot for interpolated horizontal, vertical, and time slices with raw radar values overlaid.

### Wrappers and external helpers

- `wrslib/aws_list_datetime.m`, `wrslib/mistnet_polar_radar.m`, `wrslib/radar2mat_dens.m`, `wrslib/rsl2mat_wrap.m`: wrappers or adapted helpers around upstream radar tooling.
- `regtools/get_ln.m`, `regtools/l_corner.m`, `regtools/l_curve_wo_fig.m`: local regularization helpers.
- `regtools/regu/`: vendored Regularization Tools library. Treat as external code.

## Current Shared Naming

These are the names that first-party code now uses across folder boundaries:

- `radar_scans`: cleaned scans and scan metadata
- `elevation_kernels`: precomputed elevation-specific kernels
- `grid`: Cartesian analysis grid
- `obs_matrix`: sparse observation matrix returned by `construct_F(...)`
- `kernel_idx`: scan-to-kernel index stored on each radar scan
- `nonzero_obs_idx`: linear indices of nonzero observations
- `obs_values`: extracted nonzero observation values
- `vertical_aggregation`: descriptive local name for vertical summation operators when expanded in scripts

## Shared Saved Products

- `data/derived/simulator/scans_data.mat` stores `grid`, `F`, `obs_matrix`, `elevation_kernels`, `info`, and `scans`.
- `data/derived/simulator/radar_scans.mat` stores `radar_scans`.
- `data/derived/simulator/shape_fit.mat` stores scan-cleaning shape fits.
- `data/derived/simulator/cov.mat` stores covariance fits, vertical profiles, transforms, time windows, and target scan indices.
- `info.obs_grid_size` is the canonical reshape contract for one observation vector.
- `info.obs_block_offset` stores the start offset of each observation block in the stacked system.
