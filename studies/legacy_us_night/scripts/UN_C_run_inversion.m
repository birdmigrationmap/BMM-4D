%% US night reconstruction: 4D inversion
% Reconstruct the night on the cells covered by at least one radar. The
% joint advection solve is split into consecutive time chunks sized to the
% memory budget; with enough memory the whole night is one chunk.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "us_night", "functions"))
script_tic = tic;

night_date = datetime("2026-10-21");
night_tag = string(night_date, "yyyyMMdd");
night_dir = fullfile(bmm4d.util.local_data_dir(), "us_night", "data", night_tag);

%% Controls ---------------------------------------------------------------
coverage_range_m = 160000;
flight_layer_agl_m = 3000;
min_row_mass = 0.75;  % beam weight required inside the active domain
include_water = false;
% Solver workspace only. The night's active prior (~5 GB), the observation
% system (~2-3 GB), MATLAB and the OS come on top; 16 GB keeps a 36 GB
% laptop out of swap.
memory_limit_gb = 16;
lambda_density_prior = 10;
lambda_time_density_prior = 1;
prior_block_length = 4;
density_prior_idw_neighbors = 3;
density_prior_idw_power = 2.25;
velocity_idw_neighbors = 10;
velocity_idw_power = 1.6;
prior = 0.1;
advection = 1;
pcg_tol = 1e-5;
pcg_maxit = 200;

load(fullfile(night_dir, "case.mat"), "observations", "g", "radar_scans", "case_controls")
n_time = size(radar_scans, 2);

%% Define the active domain -----------------------------------------------
% Keep columns within radar range and, in each column, the cells between
% the ground and flight_layer_agl_m above it (g.z is above sea level, so
% cells below the terrain are never solved). Water columns are excluded by
% default because the density prior is zero over water.
radar_xy = [g.f_X(cellfun(@(s) s.lon, radar_scans(:, 1))) g.f_Y(cellfun(@(s) s.lat, radar_scans(:, 1)))];
[X, Y] = ndgrid(g.x, g.y);
distance_to_radar = inf(g.sz(1:2));
for i_r = 1:size(radar_xy, 1)
    distance_to_radar = min(distance_to_radar, hypot(X - radar_xy(i_r, 1), Y - radar_xy(i_r, 2)));
end
clear X Y
active_xy = distance_to_radar <= coverage_range_m;
if ~include_water
    active_xy = active_xy & ~g.mask_water;
end
cell_top = g.z + g.dz / 2;
cell_bottom = g.z - g.dz / 2;
active_mask_3d = active_xy & cell_top > g.dem & cell_bottom < g.dem + flight_layer_agl_m;
active_id = find(active_mask_3d(:));
n_active = numel(active_id);
fprintf("[UN/C] active cells %.3g of %.3g (%.0f%% of columns), %d times\n", ...
    n_active, prod(g.sz), 100 * mean(active_xy(:)), n_time)

%% Build the density prior ------------------------------------------------
% build_radar_density_prior returns a full-grid double field, which does
% not fit in memory for the whole night, so it is built in short time
% blocks and reduced to active cells in single precision. Temporal
% smoothing of the radar profiles therefore acts within each block only.
prior_active = zeros(n_active, n_time, "single");
radar_profiles = cell(1, ceil(n_time / prior_block_length));
for i_block = 1:ceil(n_time / prior_block_length)
    time_id = (i_block - 1) * prior_block_length + 1:min(i_block * prior_block_length, n_time);
    fprintf("[UN/C] fitting radar density profiles, times %d-%d\n", time_id(1), time_id(end))
    obs_block = observations;
    obs_block.values = observations.values(:, time_id);
    [prior_field, density_prior_info] = build_radar_density_prior( ...
        g, radar_scans(:, time_id), obs_block, ...
        lambda_profile=lambda_density_prior, ...
        lambda_time_profile=lambda_time_density_prior, ...
        idw_neighbors=density_prior_idw_neighbors, ...
        idw_power=density_prior_idw_power, ...
        maxit=pcg_maxit);
    prior_field = reshape(prior_field, [], numel(time_id));
    prior_active(:, time_id) = single(prior_field(active_id, :));
    radar_profiles{i_block} = density_prior_info.radar_profiles;
    clear prior_field obs_block
end

%% Split the night into chunks that fit in memory -------------------------
% Per time slice: prior, data diagonal and rhs, and about seven PCG work
% vectors in single precision, plus two advection index vectors held twice
% (full-grid and active-domain ids). Advection mappings are built per chunk
% so that only one chunk's are in memory.
gb_per_time = (10 * 4 + 4 * 4) * n_active / 1024 ^ 3;
chunk_length = max(1, min(n_time, floor(memory_limit_gb / gb_per_time)));
chunk_start = 1:chunk_length:n_time;
fprintf("[UN/C] %.2f GB per time slice -> %d chunk(s) of up to %d times\n", ...
    gb_per_time, numel(chunk_start), chunk_length)

%% Reconstruct ------------------------------------------------------------
dt_seconds = seconds(diff(case_controls.night_controls.scan_time_requested));
active_field = zeros(n_active, n_time, "single");
inversion_info = cell(numel(chunk_start), 1);
for i_chunk = 1:numel(chunk_start)
    chunk_tic = tic;
    time_id = chunk_start(i_chunk):min(chunk_start(i_chunk) + chunk_length - 1, n_time);
    fprintf("[UN/C] chunk %d/%d | building advection mappings\n", i_chunk, numel(chunk_start))
    steps = bmm4d.advect.build_radar_advection_steps( ...
        g, radar_scans(:, time_id), dt_seconds(time_id(1:end - 1)), ...
        idw_neighbors=velocity_idw_neighbors, ...
        idw_power=velocity_idw_power, ...
        boundary_mode="open", ...
        shoreline_crossing_coeff=1, ...
        active_mask_3d=active_mask_3d);
    obs_chunk = observations;
    obs_chunk.values = observations.values(:, time_id);
    system = prepare_scans_advection_active( ...
        obs_chunk, g, active_mask_3d, steps, store_class="single", min_row_mass=min_row_mass);
    [active_field(:, time_id), inversion_info{i_chunk}] = compute_scans_advection_active( ...
        system, prior_active(:, time_id), ...
        prior=prior, advection=advection, tol=pcg_tol, maxit=pcg_maxit, nonnegative=true);
    inversion_info{i_chunk}.time_id = time_id;
    clear system obs_chunk steps
    fprintf("[UN/C] chunk %d/%d (times %d-%d) | %d iterations, flag %d | %.1f min\n", ...
        i_chunk, numel(chunk_start), time_id(1), time_id(end), ...
        inversion_info{i_chunk}.solver.iterations, inversion_info{i_chunk}.solver.flag, toc(chunk_tic) / 60)
end

%% Save -------------------------------------------------------------------
inversion_controls = struct( ...
    "coverage_range_m", coverage_range_m, ...
    "flight_layer_agl_m", flight_layer_agl_m, ...
    "min_row_mass", min_row_mass, ...
    "prior_block_length", prior_block_length, ...
    "include_water", include_water, ...
    "memory_limit_gb", memory_limit_gb, ...
    "chunk_length", chunk_length, ...
    "lambda_density_prior", lambda_density_prior, ...
    "lambda_time_density_prior", lambda_time_density_prior, ...
    "density_prior_idw_neighbors", density_prior_idw_neighbors, ...
    "density_prior_idw_power", density_prior_idw_power, ...
    "velocity_idw_neighbors", velocity_idw_neighbors, ...
    "velocity_idw_power", velocity_idw_power, ...
    "prior", prior, ...
    "advection", advection, ...
    "pcg_tol", pcg_tol, ...
    "pcg_maxit", pcg_maxit);
density_prior_profiles = single(cat(3, radar_profiles{:}));

inversion_file = fullfile(night_dir, "inversion.mat");
save(inversion_file, "active_field", "active_id", "inversion_info", "inversion_controls", ...
    "density_prior_profiles", "-v7.3")
fprintf("[UN/C] saved %s in %.1f min\n", inversion_file, toc(script_tic) / 60)
