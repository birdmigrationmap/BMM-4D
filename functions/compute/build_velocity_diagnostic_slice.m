function velocity_diagnostic = build_velocity_diagnostic_slice(g, radar_scans, truth_velocity_model, controls, varargin)
%BUILD_VELOCITY_DIAGNOSTIC_SLICE Build a compact velocity comparison slice.

p = inputParser;
addParameter(p, "z_slice_m", 350);
addParameter(p, "interval_id", []);
addParameter(p, "target_nx", 45);
parse(p, varargin{:});

n_time = size(radar_scans, 2);
if isempty(p.Results.interval_id)
    interval_id = min(ceil(n_time / 2), n_time - 1);
else
    interval_id = p.Results.interval_id;
end
[~, z_id] = min(abs(g.z - p.Results.z_slice_m));
velocity_stride = max(1, round(min(g.sz(1:2)) / p.Results.target_nx));
x_id = 1:velocity_stride:g.sz(1);
y_id = 1:velocity_stride:g.sz(2);
[X, Y] = ndgrid(g.x(x_id), g.y(y_id));
Z = g.z(z_id) + zeros(size(X));

[u_truth, v_truth] = evaluate_synthetic_velocity(X, Y, Z, truth_velocity_model, interval_id);
[u_used, v_used] = interpolate_used_velocity_slice(g, radar_scans, x_id, y_id, z_id, interval_id, controls);

if isfield(g, "mask_water")
    active_mask = ~logical(g.mask_water(x_id, y_id));
else
    active_mask = true(numel(x_id), numel(y_id));
end

velocity_diagnostic = struct( ...
    "x_id", x_id, ...
    "y_id", y_id, ...
    "lon", g.f_lon(g.x(x_id)), ...
    "lat", g.f_lat(g.y(y_id)), ...
    "z_id", z_id, ...
    "z_m", g.z(z_id), ...
    "interval_id", interval_id, ...
    "interval_start", radar_scans{1, interval_id}.datetime, ...
    "interval_end", radar_scans{1, interval_id + 1}.datetime, ...
    "active_mask", active_mask, ...
    "u_truth", single(u_truth), ...
    "v_truth", single(v_truth), ...
    "u_used", single(u_used), ...
    "v_used", single(v_used), ...
    "speed_truth_kmh", single(hypot(u_truth, v_truth) .* 3.6), ...
    "speed_used_kmh", single(hypot(u_used, v_used) .* 3.6), ...
    "speed_error_kmh", single((hypot(u_used, v_used) - hypot(u_truth, v_truth)) .* 3.6), ...
    "idw_neighbors", controls.velocity_idw_neighbors, ...
    "idw_power", controls.velocity_idw_power);
end

function [u_used, v_used] = interpolate_used_velocity_slice(g, radar_scans, x_id, y_id, z_id, interval_id, controls)

n_radar = size(radar_scans, 1);
radar_xy = nan(n_radar, 2);
uv_radar = nan(n_radar, 2);
for i_r = 1:n_radar
    radar_xy(i_r, :) = [g.f_X(radar_scans{i_r, 1}.lon) g.f_Y(radar_scans{i_r, 1}.lat)];
    uv_radar(i_r, :) = 0.5 .* (radar_scans{i_r, interval_id}.uv(z_id, :) + radar_scans{i_r, interval_id + 1}.uv(z_id, :));
end

[radar_xy_unique, ~, site_id] = unique(radar_xy, "rows", "stable");
n_site = size(radar_xy_unique, 1);
uv_site = nan(n_site, 2);
for i_site = 1:n_site
    uv_site(i_site, :) = mean(uv_radar(site_id == i_site, :), 1, "omitnan");
end

[X, Y] = ndgrid(g.x(x_id), g.y(y_id));
grid_xy = [X(:) Y(:)];
n_neighbors = min(controls.velocity_idw_neighbors, n_site);
u_used = nan(size(grid_xy, 1), 1);
v_used = nan(size(grid_xy, 1), 1);
chunk_size = 10000;
for row_start = 1:chunk_size:size(grid_xy, 1)
    row_end = min(row_start + chunk_size - 1, size(grid_xy, 1));
    chunk_ids = row_start:row_end;
    dx = grid_xy(chunk_ids, 1) - radar_xy_unique(:, 1)';
    dy = grid_xy(chunk_ids, 2) - radar_xy_unique(:, 2)';
    [distance_squared, site_ids] = mink(dx .^ 2 + dy .^ 2, n_neighbors, 2);
    weights = distance_squared .^ (-controls.velocity_idw_power / 2);
    zero_distance = distance_squared == 0;
    rows_with_zero = any(zero_distance, 2);
    weights(rows_with_zero, :) = zero_distance(rows_with_zero, :);
    weights = weights ./ sum(weights, 2);
    u_site = reshape(uv_site(site_ids(:), 1), size(site_ids));
    v_site = reshape(uv_site(site_ids(:), 2), size(site_ids));
    u_used(chunk_ids) = sum(weights .* u_site, 2);
    v_used(chunk_ids) = sum(weights .* v_site, 2);
end
u_used = reshape(u_used, numel(x_id), numel(y_id));
v_used = reshape(v_used, numel(x_id), numel(y_id));
end
