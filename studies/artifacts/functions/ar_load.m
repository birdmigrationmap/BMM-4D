function S = ar_load(c)
%AR_LOAD Load a built case for the artefact tests (dev).
%   S = AR_LOAD(C) holds the case settings, observations, grid, radar
%   scans, the active mask of the lean solver (geographic domain, ground to
%   the flight layer), the production advection steps (height above sea
%   level), the production density prior with its radar profiles, and the
%   production observation weights (estimator sd at the prior's prediction,
%   floor c.inversion.observation_error_floor), as bmm4d.invert builds them.

inv = c.inversion;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
column_mask = [];
if ~isempty(c.domain)
    column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
end
active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));
steps = bmm4d.advect.build_radar_advection_steps(g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=active_mask_3d);
[prior_active, prior_info] = bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id);
row_weights = [];
if inv.observation_weights && isfield(observations, "sd")
    row_weights = bmm4d.solve.observation_weights(observations, ...
        observations.F(:, active_id) * double(prior_active), inv.observation_error_floor);
end
S = struct("c", c, "observations", observations, "g", g, "radar_scans", {radar_scans}, ...
    "radar_selected", radar_selected, "active_mask_3d", active_mask_3d, "active_id", active_id, ...
    "steps", {steps}, "prior", prior_active, "radar_profiles", prior_info.radar_profiles, ...
    "row_weights", row_weights, "cell_air", bmm4d.grid.air_fraction(g, active_id));
end
