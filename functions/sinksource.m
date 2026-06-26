function flux = sinksource(grid, density_field, velocity_field_interval)
%SINKSOURCE Estimate local source/sink terms from density and velocity.
%   FLUX = SINKSOURCE(GRID, DENSITY_FIELD, VELOCITY_FIELD_INTERVAL) computes
%   a continuity residual from reconstructed density and horizontal velocity.
%   DENSITY_FIELD is [x y z time] in birds/km^3. VELOCITY_FIELD_INTERVAL is
%   [x y z 2 interval] in m/s, with component 1 along x and component 2
%   along y. Positive source residuals are candidate take-off, and negative
%   residuals are candidate landing.

time_start = reshape(grid.time(1:end-1), 1, []);
time_end = reshape(grid.time(2:end), 1, []);
dt_seconds_vector = seconds(time_end - time_start);
dt_hours = reshape(dt_seconds_vector ./ 3600, 1, 1, []);
dt_seconds = reshape(dt_seconds_vector, 1, 1, 1, []);

cell_area_km2 = grid.dx .* grid.dy ./ 1e6;
cell_volume_km3 = cell_area_km2 .* grid.dz ./ 1000;
face_area_x_m2 = grid.dy .* grid.dz;
face_area_y_m2 = grid.dx .* grid.dz;

density_start = density_field(:, :, :, 1:end-1);
density_end = density_field(:, :, :, 2:end);
density_mid = 0.5 .* (density_start + density_end);
density_mid_m3 = density_mid ./ 1e9;

nx = grid.sz(1);
ny = grid.sz(2);
nz = grid.sz(3);
n_interval = numel(grid.time) - 1;
u = reshape(velocity_field_interval(:, :, :, 1, :), [nx ny nz n_interval]);
v = reshape(velocity_field_interval(:, :, :, 2, :), [nx ny nz n_interval]);

flux_x_birds = zeros(nx + 1, ny, nz, n_interval);
flux_y_birds = zeros(nx, ny + 1, nz, n_interval);

flux_x_birds(1, :, :, :) = density_mid_m3(1, :, :, :) .* u(1, :, :, :) .* face_area_x_m2 .* dt_seconds;
flux_x_birds(2:nx, :, :, :) = ...
    0.5 .* (density_mid_m3(1:nx-1, :, :, :) + density_mid_m3(2:nx, :, :, :)) .* ...
    0.5 .* (u(1:nx-1, :, :, :) + u(2:nx, :, :, :)) .* face_area_x_m2 .* dt_seconds;
flux_x_birds(nx + 1, :, :, :) = density_mid_m3(nx, :, :, :) .* u(nx, :, :, :) .* face_area_x_m2 .* dt_seconds;

flux_y_birds(:, 1, :, :) = density_mid_m3(:, 1, :, :) .* v(:, 1, :, :) .* face_area_y_m2 .* dt_seconds;
flux_y_birds(:, 2:ny, :, :) = ...
    0.5 .* (density_mid_m3(:, 1:ny-1, :, :) + density_mid_m3(:, 2:ny, :, :)) .* ...
    0.5 .* (v(:, 1:ny-1, :, :) + v(:, 2:ny, :, :)) .* face_area_y_m2 .* dt_seconds;
flux_y_birds(:, ny + 1, :, :) = density_mid_m3(:, ny, :, :) .* v(:, ny, :, :) .* face_area_y_m2 .* dt_seconds;

net_outflow_xyz = ...
    flux_x_birds(2:end, :, :, :) - flux_x_birds(1:end-1, :, :, :) + ...
    flux_y_birds(:, 2:end, :, :) - flux_y_birds(:, 1:end-1, :, :);
net_outflow_birds = reshape(sum(net_outflow_xyz, 3), nx, ny, n_interval);

mass_birds = reshape(sum(density_field .* cell_volume_km3, 3), nx, ny, numel(grid.time));
mass_change_birds = diff(mass_birds, 1, 3);
source_birds = mass_change_birds + net_outflow_birds;

boundary_outflow_xyz = zeros(nx, ny, nz, n_interval);
boundary_outflow_xyz(1, :, :, :) = boundary_outflow_xyz(1, :, :, :) - flux_x_birds(1, :, :, :);
boundary_outflow_xyz(nx, :, :, :) = boundary_outflow_xyz(nx, :, :, :) + flux_x_birds(nx + 1, :, :, :);
boundary_outflow_xyz(:, 1, :, :) = boundary_outflow_xyz(:, 1, :, :) - flux_y_birds(:, 1, :, :);
boundary_outflow_xyz(:, ny, :, :) = boundary_outflow_xyz(:, ny, :, :) + flux_y_birds(:, ny + 1, :, :);
boundary_outflow_birds = reshape(sum(boundary_outflow_xyz, 3), nx, ny, n_interval);

flux = struct();
flux.source_birds = source_birds;
flux.takingoff_birds = max(source_birds, 0);
flux.landing_birds = min(source_birds, 0);
flux.boundary_flux_birds = boundary_outflow_birds;
flux.leaving_birds = max(boundary_outflow_birds, 0);
flux.entering_birds = min(boundary_outflow_birds, 0);

flux.source_birds_per_hour = source_birds ./ dt_hours;
flux.takingoff_birds_per_hour = flux.takingoff_birds ./ dt_hours;
flux.landing_birds_per_hour = flux.landing_birds ./ dt_hours;
flux.boundary_flux_birds_per_hour = boundary_outflow_birds ./ dt_hours;
flux.leaving_birds_per_hour = flux.leaving_birds ./ dt_hours;
flux.entering_birds_per_hour = flux.entering_birds ./ dt_hours;

flux.column_density_birds_per_km2 = reshape( ...
    sum(density_field, 3) .* grid.dz ./ 1000, nx, ny, numel(grid.time));
flux.mass_birds = mass_birds;
flux.mass_change_birds = mass_change_birds;
flux.net_outflow_birds = net_outflow_birds;
flux.net_outflow_birds_per_hour = net_outflow_birds ./ dt_hours;
flux.mask_water = grid.mask_water;
flux.time_start = time_start;
flux.time_end = time_end;
flux.time_mid = time_start + seconds(dt_seconds_vector ./ 2);
flux.units = struct( ...
    "density_field", "birds/km^3", ...
    "velocity_field_interval", "m/s", ...
    "source_birds", "birds per scan interval", ...
    "source_birds_per_hour", "birds/hour", ...
    "column_density_birds_per_km2", "birds/km^2", ...
    "mass_birds", "birds per horizontal grid cell", ...
    "positive_source", "candidate take-off/source", ...
    "negative_source", "candidate landing/sink", ...
    "positive_boundary_flux", "leaving domain", ...
    "negative_boundary_flux", "entering domain");
end
