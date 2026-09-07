function balance = compute_continuity_balance(g, density_field, velocity_field_interval)
%COMPUTE_CONTINUITY_BALANCE Estimate source-sink terms from density and velocity.
%   BALANCE = COMPUTE_CONTINUITY_BALANCE(G, DENSITY_FIELD,
%   VELOCITY_FIELD_INTERVAL) computes a continuity residual from
%   reconstructed density and horizontal velocity.
%   DENSITY_FIELD is the full [x y z time] volumetric density field in
%   birds/km^3, not a vertically summed column-density field.
%   VELOCITY_FIELD_INTERVAL is [x y z 2 interval] in m/s, with component 1
%   along x and component 2 along y. Output source and boundary terms are
%   vertically integrated bird counts per horizontal grid cell per scan
%   interval. The *_rate_birds_per_hour fields are the corresponding rates.
%   Signed fields retain the continuity-equation sign convention; take-off,
%   landing, entering, and leaving fields are positive magnitudes.

time_start = reshape(g.time(1:end-1), 1, []);
time_end = reshape(g.time(2:end), 1, []);
dt_seconds_vector = seconds(time_end - time_start);
dt_hours = reshape(dt_seconds_vector ./ 3600, 1, 1, []);
dt_seconds = reshape(dt_seconds_vector, 1, 1, 1, []);

cell_area_km2 = g.dx .* g.dy ./ 1e6;
cell_volume_km3 = cell_area_km2 .* g.dz ./ 1000;
face_area_x_m2 = g.dy .* g.dz;
face_area_y_m2 = g.dx .* g.dz;

density_start = density_field(:, :, :, 1:end-1);
density_end = density_field(:, :, :, 2:end);
density_mid = 0.5 .* (density_start + density_end);
density_mid_m3 = density_mid ./ 1e9;

nx = g.sz(1);
ny = g.sz(2);
nz = g.sz(3);
n_interval = numel(g.time) - 1;
u = reshape(velocity_field_interval(:, :, :, 1, :), [nx ny nz n_interval]);
v = reshape(velocity_field_interval(:, :, :, 2, :), [nx ny nz n_interval]);

transport_x_count_birds = zeros(nx + 1, ny, nz, n_interval);
transport_y_count_birds = zeros(nx, ny + 1, nz, n_interval);

transport_x_count_birds(1, :, :, :) = density_mid_m3(1, :, :, :) .* u(1, :, :, :) .* face_area_x_m2 .* dt_seconds;
transport_x_count_birds(2:nx, :, :, :) = ...
    0.5 .* (density_mid_m3(1:nx-1, :, :, :) + density_mid_m3(2:nx, :, :, :)) .* ...
    0.5 .* (u(1:nx-1, :, :, :) + u(2:nx, :, :, :)) .* face_area_x_m2 .* dt_seconds;
transport_x_count_birds(nx + 1, :, :, :) = density_mid_m3(nx, :, :, :) .* u(nx, :, :, :) .* face_area_x_m2 .* dt_seconds;

transport_y_count_birds(:, 1, :, :) = density_mid_m3(:, 1, :, :) .* v(:, 1, :, :) .* face_area_y_m2 .* dt_seconds;
transport_y_count_birds(:, 2:ny, :, :) = ...
    0.5 .* (density_mid_m3(:, 1:ny-1, :, :) + density_mid_m3(:, 2:ny, :, :)) .* ...
    0.5 .* (v(:, 1:ny-1, :, :) + v(:, 2:ny, :, :)) .* face_area_y_m2 .* dt_seconds;
transport_y_count_birds(:, ny + 1, :, :) = density_mid_m3(:, ny, :, :) .* v(:, ny, :, :) .* face_area_y_m2 .* dt_seconds;

net_outflow_xyz = ...
    transport_x_count_birds(2:end, :, :, :) - transport_x_count_birds(1:end-1, :, :, :) + ...
    transport_y_count_birds(:, 2:end, :, :) - transport_y_count_birds(:, 1:end-1, :, :);
net_outflow_count_birds = reshape(sum(net_outflow_xyz, 3), nx, ny, n_interval);

mass_birds = reshape(sum(density_field .* cell_volume_km3, 3), nx, ny, numel(g.time));
mass_change_count_birds = diff(mass_birds, 1, 3);
signed_source_count_birds = mass_change_count_birds + net_outflow_count_birds;

boundary_outflow_xyz = zeros(nx, ny, nz, n_interval);
boundary_outflow_xyz(1, :, :, :) = boundary_outflow_xyz(1, :, :, :) - transport_x_count_birds(1, :, :, :);
boundary_outflow_xyz(nx, :, :, :) = boundary_outflow_xyz(nx, :, :, :) + transport_x_count_birds(nx + 1, :, :, :);
boundary_outflow_xyz(:, 1, :, :) = boundary_outflow_xyz(:, 1, :, :) - transport_y_count_birds(:, 1, :, :);
boundary_outflow_xyz(:, ny, :, :) = boundary_outflow_xyz(:, ny, :, :) + transport_y_count_birds(:, ny + 1, :, :);
signed_boundary_count_birds = reshape(sum(boundary_outflow_xyz, 3), nx, ny, n_interval);

balance = struct();
balance.signed_source_count_birds = signed_source_count_birds;
balance.takeoff_count_birds = max(signed_source_count_birds, 0);
balance.landing_count_birds = max(-signed_source_count_birds, 0);
balance.signed_boundary_count_birds = signed_boundary_count_birds;
balance.leaving_count_birds = max(signed_boundary_count_birds, 0);
balance.entering_count_birds = max(-signed_boundary_count_birds, 0);

balance.signed_source_rate_birds_per_hour = signed_source_count_birds ./ dt_hours;
balance.takeoff_rate_birds_per_hour = balance.takeoff_count_birds ./ dt_hours;
balance.landing_rate_birds_per_hour = balance.landing_count_birds ./ dt_hours;
balance.signed_boundary_rate_birds_per_hour = signed_boundary_count_birds ./ dt_hours;
balance.leaving_rate_birds_per_hour = balance.leaving_count_birds ./ dt_hours;
balance.entering_rate_birds_per_hour = balance.entering_count_birds ./ dt_hours;

balance.column_density_birds_per_km2 = reshape( ...
    sum(density_field, 3) .* g.dz ./ 1000, nx, ny, numel(g.time));
balance.mass_birds = mass_birds;
balance.mass_change_count_birds = mass_change_count_birds;
balance.net_outflow_count_birds = net_outflow_count_birds;
balance.net_outflow_rate_birds_per_hour = net_outflow_count_birds ./ dt_hours;
balance.mask_water = g.mask_water;
balance.time_start = time_start;
balance.time_end = time_end;
balance.time_mid = time_start + seconds(dt_seconds_vector ./ 2);
balance.units = struct( ...
    "density_field", "birds/km^3", ...
    "velocity_field_interval", "m/s", ...
    "signed_source_count_birds", "birds per horizontal grid cell per scan interval", ...
    "takeoff_count_birds", "birds per horizontal grid cell per scan interval", ...
    "landing_count_birds", "birds per horizontal grid cell per scan interval", ...
    "signed_boundary_count_birds", "birds per horizontal grid cell per scan interval", ...
    "signed_source_rate_birds_per_hour", "birds per horizontal grid cell per hour", ...
    "takeoff_rate_birds_per_hour", "birds per horizontal grid cell per hour", ...
    "landing_rate_birds_per_hour", "birds per horizontal grid cell per hour", ...
    "signed_boundary_rate_birds_per_hour", "birds per horizontal grid cell per hour", ...
    "column_density_birds_per_km2", "birds/km^2", ...
    "mass_birds", "birds per horizontal grid cell", ...
    "mass_change_count_birds", "birds per horizontal grid cell per scan interval", ...
    "net_outflow_count_birds", "birds per horizontal grid cell per scan interval", ...
    "net_outflow_rate_birds_per_hour", "birds per horizontal grid cell per hour", ...
    "positive_signed_source", "candidate take-off/source", ...
    "negative_signed_source", "candidate landing/sink", ...
    "positive_signed_boundary", "leaving domain", ...
    "negative_signed_boundary", "entering domain");
end
