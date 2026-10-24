%% Production test 4c: two-level separable inversion (prototype)
% The sandbox model (MODEL_REVIEW section 7b) on a full case: column amount
% = coarse level V_L + fine deviation, 3 vertical shape patterns above the
% smoothed ground with coarse shape coefficients started from the radar
% profiles, near-range recovery per elevation group estimated from
% held-out rows, advection of the deviation. Meant for the uncut case
% ("us_east_near0_zero150": all near-range rows, no-echo gates to 150 km),
% since the observation model replaces the data cuts. Saves
% inversion_two_level.mat next to the production inversion.mat.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "two_level", "functions"))
if ~exist("case_name", "var")
    case_name = "us_east_near0_zero150";
end
c = prod_case(case_name);
inv = c.inversion;

%% Load the case ------------------------------------------------------------
stage_tic = tic;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans")
column_mask = [];
if ~isempty(c.domain)
    column_mask = build_domain_mask(g, c.domain);
end
active_mask_3d = build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));
stage_log(c, "P4c_load_case", stage_tic, sprintf("%.3g active cells, %.3g rows", numel(active_id), size(observations.F, 1)))

%% Radar prior and regional shape -------------------------------------------------
% The radar-profile prior only sets the advection relaxation (takeoff and
% landing); the radar profiles with weak vertical smoothing start the
% regional shape and centre its prior. Patterns: the us_east night EOFs
% (studies/sandbox/SB_double_layer.m).
stage_tic = tic;
prior_args = {"lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, "include_water", inv.prior_include_water};
prior_active = build_density_prior_fast(g, radar_scans, observations, prior_args{:}, active_id=active_id);
[~, profile_info] = build_density_prior_fast(g, radar_scans, observations, "lambda_profile", 1, ...
    "lambda_time_profile", inv.lambda_time_profile, "maxit", inv.maxit, "active_id", 1, "include_water", true);
patterns = load(fullfile(prod_case("us_east").data_dir, "shape_patterns.mat"));
n_pattern = 3;
shape_initial = regional_shape(g, radar_scans, profile_info.radar_profiles, patterns.shape_patterns, patterns.z_agl, ...
    n_pattern, 150000, 30000);
stage_log(c, "P4c_prior_and_shape", stage_tic)

%% Advection mappings ----------------------------------------------------------
stage_tic = tic;
steps = build_radar_advection_steps( ...
    g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, ...
    idw_power=inv.velocity_idw_power, ...
    boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, ...
    active_mask_3d=active_mask_3d);
stage_log(c, "P4c_advection_steps", stage_tic)

%% Solve ------------------------------------------------------------------------------
% Sandbox settings: fine deviation 0.2 over 30 km, level 0.03 over 600 km,
% shape prior 0.03 over 150 km, defaults of solve_two_level otherwise
% (near range per elevation, unweighted deviation amplitude, advection of
% the deviation, 3 rounds).
stage_tic = tic;
[active_field, inversion_info] = solve_two_level(observations, g, active_mask_3d, steps, prior_active, ...
    patterns.shape_patterns, patterns.z_agl, n_pattern=n_pattern, shape_sd=patterns.shape_sd, shape_initial=shape_initial, ...
    fine_prior=0.2, shape_prior=0.03, shape_length_m=150000, level_smoothness=0.03, ...
    advection=inv.advection, advection_relax_exponent=inv.advection_relax_exponent, ...
    tol=inv.tol, maxit=inv.maxit, projected_maxit=inv.projected_maxit);
stage_log(c, "P4c_solve", stage_tic, sprintf("%d PCG iterations over %d rounds, %.3g unknowns", ...
    inversion_info.solver.iterations, height(inversion_info.solver.history), ...
    numel(inversion_info.level_nodes) + numel(inversion_info.column_amount)))
disp(inversion_info.solver.history)
disp(inversion_info.near_curve)

%% Save --------------------------------------------------------------------------
stage_tic = tic;
clear observations prior_active steps
save(fullfile(c.data_dir, "inversion_two_level.mat"), "active_field", "active_id", "inversion_info", "-v7.3")
stage_log(c, "P4c_save", stage_tic)
