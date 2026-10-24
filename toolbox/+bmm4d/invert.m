function [active_field, active_id, inversion_info] = invert(c, solver_path)
%INVERT Reconstruct the density sequence of a case.
%   [ACTIVE_FIELD, ACTIVE_ID, INVERSION_INFO] = BMM4D.INVERT(C) builds the
%   density prior, the advection mappings and solves the joint advection
%   inversion (bmm4d.solve.solve_advection_lean) on the case of
%   build_observations, with the inversion settings C.inversion.
%   ACTIVE_FIELD [n_active n_time] holds the solved cells ACTIVE_ID (linear
%   indices of g); all are saved to <C.data_dir>/inversion.mat.
%
%   SOLVER_PATH:
%   - "lean" (default): active cells only (radar coverage or the case's
%     geographic domain, ground to the flight layer), single precision,
%     prior built directly on the active cells;
%   - "paper": every cell active, double precision, without advection
%     relaxation or non-negativity, as for the paper figures.

arguments
    c (1,1) struct
    solver_path (1,1) string {mustBeMember(solver_path, ["lean" "paper"])} = "lean"
end

inv = c.inversion;

stage_tic = tic;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans")
n_time = size(radar_scans, 2);
if solver_path == "lean"
    % Geographic domain when the case defines one, else radar ranges
    column_mask = [];
    if ~isempty(c.domain)
        column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
    end
    active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
else
    active_mask_3d = true(g.sz);
end
active_id = find(active_mask_3d(:));
bmm4d.util.stage_log(c, "invert_load_case", stage_tic, sprintf("%s solver, %.3g of %.3g cells active", solver_path, numel(active_id), prod(g.sz)))

%% Density prior
stage_tic = tic;
prior_args = {"lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, "include_water", inv.prior_include_water};
if solver_path == "lean"
    prior_active = bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, prior_args{:}, active_id=active_id);
else
    prior_field = bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, prior_args{:});
end
bmm4d.util.stage_log(c, "invert_density_prior", stage_tic)

%% Advection mappings
stage_tic = tic;
steps = bmm4d.advect.build_radar_advection_steps( ...
    g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, ...
    idw_power=inv.velocity_idw_power, ...
    boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, ...
    active_mask_3d=active_mask_3d);
bmm4d.util.stage_log(c, "invert_advection_steps", stage_tic)

%% Solve
stage_tic = tic;
% The paper path keeps the paper objective: no advection relaxation, no
% non-negativity.
relax_exponent = inv.advection_relax_exponent;
nonnegative = inv.nonnegative;
store_class = "single";
if solver_path == "paper"
    prior_active = reshape(prior_field, [], n_time);
    clear prior_field
    relax_exponent = 0;
    nonnegative = false;
    store_class = "double";
end
[active_field, inversion_info] = bmm4d.solve.solve_advection_lean( ...
    observations, g, active_mask_3d, steps, prior_active, ...
    prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, ...
    coeff_horver=inv.coeff_horver, tol=inv.tol, maxit=inv.maxit, store_class=store_class, ...
    advection_relax_exponent=relax_exponent, ...
    nonnegative=nonnegative, projected_maxit=inv.projected_maxit);
active_field = single(active_field);
bmm4d.util.stage_log(c, "invert_solve", stage_tic, sprintf("%d iterations, flag %d, relres %.2g, %.3g unknowns", ...
    inversion_info.solver.iterations, inversion_info.solver.flag, inversion_info.solver.relres, numel(active_field)))

%% Save
% Only active cells are stored; export_maps expands them to the grid.
stage_tic = tic;
inversion_controls = inv;
save(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id", "solver_path", ...
    "inversion_info", "inversion_controls", "-v7.3")
bmm4d.util.stage_log(c, "invert_save", stage_tic)
end
