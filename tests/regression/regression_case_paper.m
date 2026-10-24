function out = regression_case_paper(raw_files, kernel_dir, s)
%REGRESSION_CASE_PAPER Paper path on a small Key West case, today's code.
%   Same method and settings S as regression_case_paper_v1 (paper cleaning,
%   forward operator, IDW land-profile prior, advection steps, full-grid
%   solve in double precision without relaxation or non-negativity),
%   computed with today's functions: build_density_prior_fast and
%   solve_advection_lean with every cell active. OUT holds the compared
%   outputs.

%% Clean the volumes --------------------------------------------------------
radar_scans = bmm4d.radar.read_rd(raw_files, s.read_options{:});
radar_scans = reshape(radar_scans, 1, []);
% Sweeps differ in size: one flat list of all sweeps of all volumes
out.sweep_data = {};
for i_scan = 1:numel(radar_scans)
    out.sweep_data = [out.sweep_data; {radar_scans{i_scan}.dens.sweeps.data}'];
end

%% Build the observation system --------------------------------------------
% load_F_elev also adds each volume's velocity profile (uv)
[elevation_kernels, radar_scans] = bmm4d.forward.load_F_elev(kernel_dir, radar_scans);
out.uv = cellfun(@(scan) scan.uv, radar_scans(:), UniformOutput=false);
g = bmm4d.grid.create_grid(elevation_kernels, radar_scans, s.grid_bounds(1), s.grid_bounds(2), ...
    s.grid_bounds(3), s.grid_bounds(4), s.heightmax);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
g.mask_water(inpolygon(lon_grid, lat_grid, s.water_override(1, :), s.water_override(2, :))) = true;
[F, obs_matrix, info] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);
observations = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info, min_row_mass=s.min_row_mass);
out.grid_size = g.sz;
out.mask_water = g.mask_water;
out.F = observations.F;
out.values = observations.values;

%% Prior, advection, solve ---------------------------------------------------
inv = s.inversion;
out.prior = single(bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, ...
    lambda_profile=inv.lambda_profile, lambda_time_profile=inv.lambda_time_profile, ...
    idw_neighbors=inv.density_prior_idw_neighbors, idw_power=inv.density_prior_idw_power, maxit=inv.maxit));
dt_seconds = seconds(diff(cellfun(@(scan) scan.datetime, radar_scans(1, :))));
steps = bmm4d.advect.build_radar_advection_steps(g, radar_scans, dt_seconds, ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, ...
    boundary_mode="open", shoreline_crossing_coeff=inv.shoreline_crossing_coeff);
field = bmm4d.solve.solve_advection_lean(observations, g, true(g.sz), steps, reshape(double(out.prior), [], size(out.prior, 4)), ...
    prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, coeff_horver=inv.coeff_horver, ...
    maxit=inv.maxit, tol=1e-6, store_class="double");
out.field = single(reshape(field, size(out.prior)));
end
