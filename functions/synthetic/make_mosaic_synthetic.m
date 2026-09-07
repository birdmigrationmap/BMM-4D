function [truth_sequences, obs_sequences, mosaic] = make_mosaic_synthetic(elevation_kernels, template_scan, radar_table, cov_parm, vp_case, trans_finv_case, target_id, n_sim, varargin)
%MAKE_MOSAIC_SYNTHETIC Generate synthetic mosaic truth and radar measurements.

p = inputParser;
addParameter(p, "time_datetimes", datetime("2020-10-01 22:00:00"):minutes(30):datetime("2020-10-02 05:00:00"));
addParameter(p, "grid_bounds", []);
addParameter(p, "heightmax", []);
addParameter(p, "seed", "shuffle");
addParameter(p, "large_scale_range_km", 250);
addParameter(p, "large_scale_log_sd", 0.6);
addParameter(p, "velocity_model", "smooth");
addParameter(p, "base_speed_mps", 12);
addParameter(p, "base_direction_deg", 35);
addParameter(p, "speed_variation_mps", 2);
addParameter(p, "direction_variation_deg", 8);
addParameter(p, "velocity_range_km", 250);
addParameter(p, "velocity_sample_min_range_km", 5);
addParameter(p, "velocity_sample_max_range_km", 37.5);
addParameter(p, "velocity_sample_spacing_km", 5);
addParameter(p, "speed_padding_safety_km", 20);
addParameter(p, "temporal_floor", 0);
addParameter(p, "temporal_steepness", 1);
addParameter(p, "observation_constructor", "global");
addParameter(p, "use_water_mask", false);
addParameter(p, "use_conus_mask", false);
addParameter(p, "min_range_m", 5000);
addParameter(p, "max_range_m", Inf);
addParameter(p, "max_alt_m", Inf);
addParameter(p, "data_res_m", []);
parse(p, varargin{:});

seed = p.Results.seed;
if isstring(seed) || ischar(seed)
    rng(seed);
    seed = randi(2^31 - 1);
end

time_datetimes = p.Results.time_datetimes(:)';
i_center = ceil(numel(time_datetimes) / 2);
dt_hours = hours(time_datetimes(2) - time_datetimes(1));
time_offsets = hours(time_datetimes - time_datetimes(i_center)) ./ dt_hours;

% One full template volume defines the shared scan geometry and valid gates.
radar_scan_geometry = build_mosaic_scan_geometry(template_scan, radar_table, time_datetimes(1), p.Results);

if isempty(p.Results.grid_bounds)
    if isempty(p.Results.heightmax)
        g = create_grid(elevation_kernels, radar_scan_geometry);
    else
        g = create_grid(elevation_kernels, radar_scan_geometry);
        g.z = g.z(g.z <= p.Results.heightmax);
        g.sz = [numel(g.x) numel(g.y) numel(g.z)];
        g.dz = median(diff(g.z));
    end
elseif isempty(p.Results.heightmax)
    g = create_grid( ...
        elevation_kernels, radar_scan_geometry, ...
        p.Results.grid_bounds(1), p.Results.grid_bounds(2), ...
        p.Results.grid_bounds(3), p.Results.grid_bounds(4));
else
    g = create_grid( ...
        elevation_kernels, radar_scan_geometry, ...
        p.Results.grid_bounds(1), p.Results.grid_bounds(2), ...
        p.Results.grid_bounds(3), p.Results.grid_bounds(4), ...
        p.Results.heightmax);
end
g.time = time_datetimes;
if ~p.Results.use_water_mask
    g.mask_water = false(g.sz(1:2));
end
g.mask_water_raw = g.mask_water;
if p.Results.use_conus_mask
    g.mask_conus = get_conus_mask(g);
    g.mask_water = g.mask_water | ~g.mask_conus;
else
    g.mask_conus = true(g.sz(1:2));
end
g.dem = zeros(g.sz(1:2));
if string(p.Results.observation_constructor) == "compact"
    [F_mosaic, ~, info_mosaic] = construct_F_compact(radar_scan_geometry, elevation_kernels, g);
else
    [F_mosaic, ~, info_mosaic] = construct_F(radar_scan_geometry, elevation_kernels, g);
end

% Saved scans contain only the velocity measurements needed by inversion.
radar_scans_mosaic = build_radar_sample_table(radar_table, time_datetimes, g.sz(3));

max_duration_seconds = max(abs(seconds(time_datetimes - time_datetimes(i_center))));
max_speed_mps = p.Results.base_speed_mps + 3 * p.Results.speed_variation_mps;
pad_m = max_speed_mps * max_duration_seconds + p.Results.speed_padding_safety_km * 1000;
generation_grid = pad_generation_grid(g, pad_m);

large_scale_factor_padded = make_large_scale_factor(generation_grid, p.Results.large_scale_range_km * 1000, p.Results.large_scale_log_sd, seed);
velocity_config = make_velocity_config(time_datetimes, p.Results, generation_grid, seed);
temporal_factor = make_temporal_factor(time_datetimes, p.Results);

truth_fields = simulate_gaussian_truth( ...
    generation_grid, cov_parm, trans_finv_case{target_id}, vp_case(1:generation_grid.sz(3), target_id), ...
    n_sim, seed=seed);

truth_sequences = cell(n_sim, 1);
obs_sequences = cell(n_sim, 1);
for i_sim = 1:n_sim
    % Transport the spatial density pattern, then apply night-wide scaling.
    reference_field = truth_fields{i_sim} .* large_scale_factor_padded;
    truth_sequence_padded = advect_density_sequence_velocity_field( ...
        reference_field, generation_grid, velocity_config, time_datetimes);

    for i_t = 1:numel(time_datetimes)
        truth_sequence_padded(:, :, :, i_t) = truth_sequence_padded(:, :, :, i_t) .* temporal_factor(i_t);
    end

    truth_sequence = truth_sequence_padded(generation_grid.crop.x_idx, generation_grid.crop.y_idx, :, :);
    inactive_mask_3d = repmat(g.mask_water, 1, 1, g.sz(3));
    for i_t = 1:numel(time_datetimes)
        truth_t = truth_sequence(:, :, :, i_t);
        truth_t(inactive_mask_3d) = 0;
        truth_sequence(:, :, :, i_t) = truth_t;
    end
    truth_sequences{i_sim} = truth_sequence;
    obs_sequences{i_sim} = forward_project_sequence(F_mosaic, truth_sequence);
end

% Sample the true velocity over each radar's VVP footprint and scan times.
uv_profiles_radar = extract_radar_velocity_profiles( ...
    generation_grid, radar_scans_mosaic, velocity_config, ...
    p.Results.velocity_sample_min_range_km * 1000, ...
    p.Results.velocity_sample_max_range_km * 1000, ...
    p.Results.velocity_sample_spacing_km * 1000);
for i_r = 1:size(radar_scans_mosaic, 1)
    for i_t = 1:size(radar_scans_mosaic, 2)
        radar_scans_mosaic{i_r, i_t}.uv = uv_profiles_radar(:, :, i_r, i_t);
    end
end

mosaic = struct( ...
    "g", g, ...
    "F", F_mosaic, ...
    "info", info_mosaic, ...
    "radar_scans", {radar_scans_mosaic}, ...
    "time_offsets", time_offsets, ...
    "time_datetimes", time_datetimes, ...
    "dt", median(diff(time_datetimes)), ...
    "seed", seed, ...
    "temporal_factor", temporal_factor, ...
    "truth_velocity_model", velocity_config, ...
    "velocity_sampling", struct( ...
        "min_range_km", p.Results.velocity_sample_min_range_km, ...
        "max_range_km", p.Results.velocity_sample_max_range_km, ...
        "spacing_km", p.Results.velocity_sample_spacing_km), ...
    "large_scale_range_km", p.Results.large_scale_range_km, ...
    "large_scale_log_sd", p.Results.large_scale_log_sd, ...
    "radar_table", radar_table, ...
    "grid_bounds", p.Results.grid_bounds);
end

function radar_scan_geometry = build_mosaic_scan_geometry(template_scan, radar_table, scan_datetime, controls)

n_radar = height(radar_table);
radar_scan_geometry = cell(n_radar, 1);
template_scan = thin_synthetic_scan_geometry(template_scan, controls);
for i_r = 1:n_radar
    radar_scan_geometry{i_r} = template_scan;
    radar_scan_geometry{i_r}.station = radar_table.name(i_r);
    radar_scan_geometry{i_r}.lat = radar_table.lat(i_r);
    radar_scan_geometry{i_r}.lon = radar_table.lon(i_r);
    radar_scan_geometry{i_r}.height = 0;
    radar_scan_geometry{i_r}.datetime = scan_datetime;
end
end

function radar_scan = thin_synthetic_scan_geometry(radar_scan, controls)

data_res_m = controls.data_res_m;
if isempty(data_res_m) && isinf(controls.max_range_m) && isinf(controls.max_alt_m)
    return
end
if isempty(data_res_m)
    data_res_m = [0 0];
end

for i_elev = 1:numel(radar_scan.dens.sweeps)
    sweep = radar_scan.dens.sweeps(i_elev);
    [~, sweep_range, sweep_az] = sweep2mat(sweep);
    [~, sweep_height] = slant2ground(sweep_range, sweep.elev);

    keep = true(numel(sweep_range), numel(sweep_az));
    keep(sweep_range < controls.min_range_m | sweep_range > controls.max_range_m, :) = false;
    keep(sweep_height > controls.max_alt_m, :) = false;

    if any(data_res_m > 0)
        res_range = max(1, floor(data_res_m(1) ./ median(diff(sweep_range))));
        res_az = max(1, floor(data_res_m(2) ./ (deg2rad(median(diff(sweep_az))) .* sweep_range)));
        keep_res = false(size(keep));
        for i_range = 1:res_range:size(keep_res, 1)
            keep_res(i_range, 1:res_az(i_range):end) = true;
        end
        keep = keep & keep_res;
    end

    sweep_data = ones(size(sweep.data));
    sweep_data(~keep) = nan;
    radar_scan.dens.sweeps(i_elev).data = sweep_data;
end
end

function radar_scans = build_radar_sample_table(radar_table, time_datetimes, nz)

n_radar = height(radar_table);
n_time = numel(time_datetimes);
radar_scans = cell(n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        radar_scans{i_r, i_t} = struct( ...
            "station", radar_table.name(i_r), ...
            "lat", radar_table.lat(i_r), ...
            "lon", radar_table.lon(i_r), ...
            "height", 0, ...
            "datetime", time_datetimes(i_t), ...
            "uv", zeros(nz, 2));
    end
end
end

function generation_grid = pad_generation_grid(g, pad_m)

pad_x = ceil(pad_m / g.dx);
pad_y = ceil(pad_m / g.dy);

generation_grid = g;
generation_grid.x = g.x(1) + ((-pad_x):(g.sz(1) + pad_x - 1)) .* g.dx;
generation_grid.y = g.y(1) + ((-pad_y):(g.sz(2) + pad_y - 1))' .* g.dy;
generation_grid.sz = [numel(generation_grid.x) numel(generation_grid.y) numel(generation_grid.z)];
generation_grid.lat = generation_grid.f_lat(generation_grid.y);
generation_grid.lon = generation_grid.f_lon(generation_grid.x);
if isfield(g, "mask_water") && any(g.mask_water(:))
    generation_grid.mask_water = get_water_mask(generation_grid);
else
    generation_grid.mask_water = false(generation_grid.sz(1:2));
end
generation_grid.dem = zeros(generation_grid.sz(1:2));
generation_grid.crop = struct( ...
    "x_idx", (1:g.sz(1)) + pad_x, ...
    "y_idx", (1:g.sz(2)) + pad_y, ...
    "pad_x", pad_x, ...
    "pad_y", pad_y);
end

function factor = make_large_scale_factor(g, range_m, log_sd, seed)

if log_sd == 0
    factor = ones(g.sz(1:2));
    return
end

rng(seed + 1000, "twister");
coarse_x = g.x(1):range_m:g.x(end);
coarse_y = g.y(1):range_m:g.y(end);
coarse_x = unique([coarse_x g.x(end)]);
coarse_y = unique([coarse_y g.y(end)]);
coarse_noise = randn(numel(coarse_x), numel(coarse_y));
F = griddedInterpolant({coarse_x, coarse_y}, coarse_noise, "linear", "nearest");
smooth_noise = F({g.x(:), g.y(:)});
smooth_noise = (smooth_noise - mean(smooth_noise(:))) ./ max(std(smooth_noise(:)), eps);
factor = exp(log_sd .* smooth_noise);
factor = factor ./ mean(factor(:));
end

function velocity_config = make_velocity_config(time_datetimes, controls, g, seed)

rng(seed + 2000, "twister");

velocity_config = struct( ...
    "model", string(controls.velocity_model), ...
    "n_interval", numel(time_datetimes) - 1, ...
    "base_direction_rad", deg2rad(controls.base_direction_deg), ...
    "direction_variation_rad", deg2rad(controls.direction_variation_deg), ...
    "range_m", controls.velocity_range_km * 1000, ...
    "base_speed_mps", controls.base_speed_mps, ...
    "speed_variation_mps", controls.speed_variation_mps);

if velocity_config.model == "uniform"
    return
end

field_scale_m = velocity_config.range_m;
node_spacing_xy_m = max([g.dx * 4, g.dy * 4, field_scale_m / 2]);
node_spacing_z_m = max(400, min(1000, field_scale_m / 250));

coarse_x = build_coarse_axis(g.x(1), g.x(end), node_spacing_xy_m);
coarse_y = build_coarse_axis(g.y(1), g.y(end), node_spacing_xy_m);
coarse_z = build_coarse_axis(g.z(1), g.z(end), node_spacing_z_m);

speed_nodes = randn(numel(coarse_x), numel(coarse_y), numel(coarse_z), velocity_config.n_interval);
direction_nodes = randn(numel(coarse_x), numel(coarse_y), numel(coarse_z), velocity_config.n_interval);
speed_nodes = normalize_gaussian_nodes(smooth_gaussian_nodes(speed_nodes));
direction_nodes = normalize_gaussian_nodes(smooth_gaussian_nodes(direction_nodes));

mean_direction_rad = velocity_config.base_direction_rad + 0.5 * velocity_config.direction_variation_rad * randn;
speed_fields = max(1, velocity_config.base_speed_mps + velocity_config.speed_variation_mps .* speed_nodes);
direction_fields = mean_direction_rad + velocity_config.direction_variation_rad .* direction_nodes;

velocity_config.coarse_x = coarse_x;
velocity_config.coarse_y = coarse_y;
velocity_config.coarse_z = coarse_z;
velocity_config.mean_direction_rad = mean_direction_rad;
velocity_config.node_spacing_xy_m = node_spacing_xy_m;
velocity_config.node_spacing_z_m = node_spacing_z_m;
velocity_config.speed_fields = speed_fields;
velocity_config.direction_fields = direction_fields;
end

function axis_values = build_coarse_axis(axis_min, axis_max, spacing)

if axis_min == axis_max
    axis_values = axis_min;
    return
end

axis_values = axis_min:spacing:axis_max;
axis_values = unique([axis_values axis_max]);
end

function nodes = smooth_gaussian_nodes(nodes)

kernel = [1 2 1] ./ 4;
nodes = smooth_nodes_along_dim(nodes, 1, kernel);
nodes = smooth_nodes_along_dim(nodes, 2, kernel);
nodes = smooth_nodes_along_dim(nodes, 3, kernel);
nodes = smooth_nodes_along_dim(nodes, 4, kernel);
nodes = smooth_nodes_along_dim(nodes, 1, kernel);
nodes = smooth_nodes_along_dim(nodes, 2, kernel);
end

function nodes = smooth_nodes_along_dim(nodes, dim, kernel)

if size(nodes, dim) == 1
    return
end

kernel_shape = ones(1, ndims(nodes));
kernel_shape(dim) = numel(kernel);
nodes = convn(nodes, reshape(kernel, kernel_shape), "same");
end

function nodes = normalize_gaussian_nodes(nodes)

nodes = nodes - mean(nodes, "all");
nodes = nodes ./ max(std(nodes, 0, "all"), eps);
end

function temporal_factor = make_temporal_factor(time_datetimes, controls)

if isscalar(time_datetimes)
    temporal_factor = 1;
else
    time_fraction = linspace(0, 1, numel(time_datetimes));
    temporal_factor = sin(pi .* time_fraction) .^ controls.temporal_steepness;
end
temporal_factor = controls.temporal_floor + (1 - controls.temporal_floor) .* temporal_factor;
end

function fields = advect_density_sequence_velocity_field(reference_field, g, velocity_config, time_datetimes)

n_time = numel(time_datetimes);
n_interval = n_time - 1;
fields = zeros([size(reference_field) n_time], "like", reference_field);
center_id = ceil(n_time / 2);
[X, Y, Z] = ndgrid(g.x(:), g.y(:), g.z(:));

fields(:, :, :, center_id) = reference_field;

for i_t = center_id + 1:n_time
    uv_interval = make_velocity_interval(X, Y, Z, g.sz, velocity_config, i_t - 1);
    fields(:, :, :, i_t) = advect_one_step(fields(:, :, :, i_t - 1), g, X, Y, Z, uv_interval, seconds(time_datetimes(i_t) - time_datetimes(i_t - 1)), -1);
end

for i_t = center_id - 1:-1:1
    uv_interval = make_velocity_interval(X, Y, Z, g.sz, velocity_config, i_t);
    fields(:, :, :, i_t) = advect_one_step(fields(:, :, :, i_t + 1), g, X, Y, Z, uv_interval, seconds(time_datetimes(i_t + 1) - time_datetimes(i_t)), 1);
end

assert(velocity_config.n_interval == n_interval)
end

function uv_interval = make_velocity_interval(X, Y, Z, grid_sz, velocity_config, i_interval)

[u, v] = evaluate_synthetic_velocity(X, Y, Z, velocity_config, i_interval);
uv_interval = zeros([grid_sz 2]);
uv_interval(:, :, :, 1) = u;
uv_interval(:, :, :, 2) = v;
end

function field_out = advect_one_step(field_in, g, X, Y, Z, uv_interval, dt_seconds, direction)

F = griddedInterpolant({g.x(:), g.y(:), g.z(:)}, field_in, "linear", "none");
query_x = X + direction .* uv_interval(:, :, :, 1) .* dt_seconds;
query_y = Y + direction .* uv_interval(:, :, :, 2) .* dt_seconds;
field_out = F(query_x, query_y, Z);
field_out(~isfinite(field_out)) = 0;
end

function obs_sequence = forward_project_sequence(F, truth_sequence)

n_time = size(truth_sequence, 4);
obs_rows = cell(n_time, 1);
obs_cols = cell(n_time, 1);
obs_vals = cell(n_time, 1);
for i_t = 1:n_time
    field_t = truth_sequence(:, :, :, i_t);
    obs_t = F * field_t(:);
    obs_rows{i_t} = find(obs_t);
    obs_cols{i_t} = repmat(i_t, numel(obs_rows{i_t}), 1);
    obs_vals{i_t} = obs_t(obs_rows{i_t});
end
obs_sequence = sparse(vertcat(obs_rows{:}), vertcat(obs_cols{:}), vertcat(obs_vals{:}), size(F, 1), n_time);
end

function uv_profiles_radar = extract_radar_velocity_profiles(g, radar_scans, velocity_config, min_range_m, max_range_m, spacing_m)

n_radar = size(radar_scans, 1);
n_time = size(radar_scans, 2);
nz = g.sz(3);
n_interval = velocity_config.n_interval;
uv_profiles_interval = nan(nz, 2, n_radar, n_interval);
radar_xy = get_radar_xy(g, radar_scans(:, 1));
sample_axis = -max_range_m:spacing_m:max_range_m;
[sample_dx, sample_dy] = ndgrid(sample_axis, sample_axis);
sample_range = hypot(sample_dx, sample_dy);
sample_id = sample_range >= min_range_m & sample_range <= max_range_m;
sample_dx = sample_dx(sample_id);
sample_dy = sample_dy(sample_id);

X = reshape(radar_xy(:, 1), 1, n_radar, 1) + reshape(sample_dx, 1, 1, []);
Y = reshape(radar_xy(:, 2), 1, n_radar, 1) + reshape(sample_dy, 1, 1, []);
Z = g.z(:);

for i_interval = 1:n_interval
    [u, v] = evaluate_synthetic_velocity(X, Y, Z, velocity_config, i_interval);
    uv_profiles_interval(:, 1, :, i_interval) = reshape(mean(u, 3), nz, 1, n_radar);
    uv_profiles_interval(:, 2, :, i_interval) = reshape(mean(v, 3), nz, 1, n_radar);
end

uv_profiles_radar = nan(nz, 2, n_radar, n_time);
uv_profiles_radar(:, :, :, 1) = uv_profiles_interval(:, :, :, 1);
uv_profiles_radar(:, :, :, n_time) = uv_profiles_interval(:, :, :, n_interval);
for i_t = 2:(n_time - 1)
    uv_profiles_radar(:, :, :, i_t) = 0.5 .* (uv_profiles_interval(:, :, :, i_t - 1) + uv_profiles_interval(:, :, :, i_t));
end
end

function radar_xy = get_radar_xy(g, radar_scans)

n_radar = numel(radar_scans);
radar_xy = nan(n_radar, 2);
for i_r = 1:n_radar
    radar_xy(i_r, :) = [g.f_X(radar_scans{i_r}.lon) g.f_Y(radar_scans{i_r}.lat)];
end
end
