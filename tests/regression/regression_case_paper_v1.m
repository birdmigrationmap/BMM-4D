function out = regression_case_paper_v1(raw_files, kernel_dir, s)
%REGRESSION_CASE_PAPER_V1 Paper path on a small Key West case, paper-v1 code.
%   Runs read_rd with the paper cleaning (rhoHV rain removal), the forward
%   operator, the IDW radar-profile prior, the advection steps and the
%   full-grid solver, with the settings S (see make_regression_reference).
%   Calls the paper-v1 functions; only make_regression_reference runs it,
%   on a checkout of paper-v1. regression_case_paper is the same method on
%   today's code. OUT holds the compared outputs.

%% Clean the volumes --------------------------------------------------------
radar_scans = read_rd(raw_files, s.read_options{:});
radar_scans = reshape(radar_scans, 1, []);
% Sweeps differ in size: one flat list of all sweeps of all volumes
out.sweep_data = {};
for i_scan = 1:numel(radar_scans)
    out.sweep_data = [out.sweep_data; {radar_scans{i_scan}.dens.sweeps.data}'];
end

%% Build the observation system --------------------------------------------
% load_F_elev also adds each volume's velocity profile (uv)
[elevation_kernels, radar_scans] = load_F_elev(kernel_dir, radar_scans);
out.uv = cellfun(@(scan) scan.uv, radar_scans(:), UniformOutput=false);
g = create_grid(elevation_kernels, radar_scans, s.grid_bounds(1), s.grid_bounds(2), ...
    s.grid_bounds(3), s.grid_bounds(4), s.heightmax);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
g.mask_water(inpolygon(lon_grid, lat_grid, s.water_override(1, :), s.water_override(2, :))) = true;
[F, obs_matrix, info] = construct_F(radar_scans, elevation_kernels, g);
observations = prepare_radar_observations(F, obs_matrix, info, min_row_mass=s.min_row_mass);
out.grid_size = g.sz;
out.mask_water = g.mask_water;
out.F = observations.F;
out.values = observations.values;

%% Prior, advection, solve ---------------------------------------------------
inv = s.inversion;
out.prior = single(build_radar_density_prior(g, radar_scans, observations, ...
    lambda_profile=inv.lambda_profile, lambda_time_profile=inv.lambda_time_profile, ...
    idw_neighbors=inv.density_prior_idw_neighbors, idw_power=inv.density_prior_idw_power, maxit=inv.maxit));
dt_seconds = seconds(diff(cellfun(@(scan) scan.datetime, radar_scans(1, :))));
steps = build_radar_advection_steps(g, radar_scans, dt_seconds, ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, ...
    boundary_mode="open", shoreline_crossing_coeff=inv.shoreline_crossing_coeff);
system = prepare_scans_advection(observations, g, steps, water=true, coeff_horver=inv.coeff_horver);
field = compute_scans_advection(system, double(out.prior), ...
    prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, maxit=inv.maxit);
out.field = single(field);
end
