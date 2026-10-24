%% Production test 4: reconstruct the density sequence
% Density prior, advection mappings and the joint advection inversion, with
% each stage timed. Two solver paths:
% - "lean" (default): active cells only (radar coverage, ground to the flight
%   layer), single precision, all times batched, prior built directly on the
%   active cells; same objective and term weights as the paper.
% - "paper": full grid in double precision (prepare_scans_advection +
%   batched compute_scans_advection), as for the paper figures.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"))
% Case: "keywest", "midwest" or "us_east"; suffix _rho or _none for other rain
% removal. A case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
if ~exist("solver_path", "var")
    solver_path = "lean";   % "lean" or "paper"
end
c = prod_case(case_name);
inv = c.inversion;

stage_tic = tic;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans")
n_time = size(radar_scans, 2);
if solver_path == "lean"
    % Geographic domain when the case defines one, else radar ranges
    column_mask = [];
    if ~isempty(c.domain)
        column_mask = build_domain_mask(g, c.domain);
    end
    active_mask_3d = build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
else
    active_mask_3d = true(g.sz);
end
active_id = find(active_mask_3d(:));
stage_log(c, "P4_load_case", stage_tic, sprintf("%s solver, %.3g of %.3g cells active", solver_path, numel(active_id), prod(g.sz)))

%% Density prior ----------------------------------------------------------------
stage_tic = tic;
prior_args = {"lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, "include_water", inv.prior_include_water};
if solver_path == "lean"
    prior_active = build_density_prior_fast(g, radar_scans, observations, prior_args{:}, active_id=active_id);
else
    prior_field = build_density_prior_fast(g, radar_scans, observations, prior_args{:});
end
stage_log(c, "P4_density_prior", stage_tic)

%% Advection mappings ----------------------------------------------------------
stage_tic = tic;
steps = build_radar_advection_steps( ...
    g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, ...
    idw_power=inv.velocity_idw_power, ...
    boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, ...
    active_mask_3d=active_mask_3d);
stage_log(c, "P4_advection_steps", stage_tic)

%% Solve ------------------------------------------------------------------------------
stage_tic = tic;
if solver_path == "lean"
    [active_field, inversion_info] = solve_advection_lean( ...
        observations, g, active_mask_3d, steps, prior_active, ...
        prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, ...
        coeff_horver=inv.coeff_horver, tol=inv.tol, maxit=inv.maxit, store_class="single", ...
        advection_relax_exponent=inv.advection_relax_exponent, ...
        nonnegative=inv.nonnegative, projected_maxit=inv.projected_maxit);
else
    system = prepare_scans_advection(observations, g, steps, water=true, coeff_horver=inv.coeff_horver);
    [active_field, inversion_info] = solve_advection_fast( ...
        system, prior_field, ...
        prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, tol=inv.tol, maxit=inv.maxit);
    active_field = single(reshape(active_field, [], n_time));
    clear system prior_field
end
stage_log(c, "P4_solve", stage_tic, sprintf("%d iterations, flag %d, relres %.2g, %.3g unknowns", ...
    inversion_info.solver.iterations, inversion_info.solver.flag, inversion_info.solver.relres, numel(active_field)))

%% Save --------------------------------------------------------------------------
% Only active cells are stored; P5 expands them to the grid.
stage_tic = tic;
inversion_controls = inv;
save(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id", "solver_path", ...
    "inversion_info", "inversion_controls", "-v7.3")
stage_log(c, "P4_save", stage_tic)
