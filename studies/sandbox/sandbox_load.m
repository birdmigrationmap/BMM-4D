function S = sandbox_load(case_name)
%SANDBOX_LOAD Load a prod case once for fast option tests (sandbox_run).
%   S holds the case settings, observations (with row range, elevation and
%   a no-echo flag), grid, active mask, advection steps, radar table and an
%   in-memory cache of density priors keyed by their settings.

c = case_settings(case_name);
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
column_mask = [];
if ~isempty(c.domain)
    column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
end
active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
inv = c.inversion;
steps = bmm4d.advect.build_radar_advection_steps(g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=active_mask_3d);
% No-echo rows: every observed value of the row is a zero (stored as eps)
[r, ~, v] = find(observations.values);
row_max = accumarray(r, double(v), [size(observations.values, 1) 1], @max, 0);
S = struct("c", c, "observations", observations, "g", g, "radar_scans", {radar_scans}, ...
    "radar_selected", radar_selected, "active_mask_3d", active_mask_3d, ...
    "active_id", find(active_mask_3d(:)), "steps", {steps}, ...
    "row_is_zero", row_max < 1e-6, "prior_cache", containers.Map());
end
