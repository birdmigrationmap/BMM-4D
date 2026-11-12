function solve_info = ar_invert_variant(night_name, variant, prior_options, holdout)
%AR_INVERT_VARIANT Solve one night of the night set with one prior variant (dev, cluster job).
%   SOLVE_INFO = AR_INVERT_VARIANT(NIGHT_NAME, VARIANT, PRIOR_OPTIONS,
%   HOLDOUT) solves the us_east night NIGHT_NAME (ar_night_config; its
%   observations must exist, ar_night_stage1) as bmm4d.invert (lean path:
%   profiles and links by height above the smoothed ground, observation
%   weights, advection relaxation, non-negativity), with the prior options
%   PRIOR_OPTIONS of ar_solve_advection_lean (name-value cell: smoothness,
%   prior_correlation_m, prior_correlation_vertical_m, prior_matern_m).
%   maxit is raised to 5000 so that every variant converges to the
%   production tolerance.
%   HOLDOUT true leaves out the radars of ar_holdout_radars: their rows get
%   weight 0, the profile field is rebuilt without them, and their rows at
%   30-120 km are then predicted (relative L2 error and median absolute
%   log ratio, both floored at 0.1 birds km^-3, per radar).
%   Saves inversion_<VARIANT>.mat or inversion_<VARIANT>_holdout.mat
%   (active_field, active_id, solve_info) in the night's data folder.

project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
c = ar_night_config(night_name);
inv = c.inversion;
inv.maxit = 5000;

%% Case, active cells, held-out radars
stage_tic = tic;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));
keep_radar = true(size(radar_scans, 1), 1);
if holdout
    keep_radar = ~ismember(string(radar_selected.name), ar_holdout_radars());
end
held_row = ~keep_radar(observations.row_radar_id);

%% Prior mean from the kept radars only (as bmm4d.invert, lean)
% Held-out rows are dropped and the radar index of the others renumbered
new_index = cumsum(keep_radar);
observations_prior = observations;
observations_prior.F = observations.F(~held_row, :);
observations_prior.values = observations.values(~held_row, :);
observations_prior.row_radar_id = new_index(observations.row_radar_id(~held_row));
prior_active = bmm4d.prior.build_density_prior_fast(g, radar_scans(keep_radar, :), observations_prior, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id, "ground_smoothing_m", inv.ground_smoothing_m);
clear observations_prior
% The winds of the advection still come from all radars
steps = bmm4d.advect.build_radar_advection_steps(g, radar_scans, seconds(diff(c.times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=active_mask_3d, ...
    ground_smoothing_m=inv.ground_smoothing_m);
row_weights = double(observations.values ~= 0);
if inv.observation_weights && isfield(observations, "sd")
    row_weights = bmm4d.solve.observation_weights(observations, ...
        observations.F(:, active_id) * double(prior_active), inv.observation_error_floor);
end
row_weights(held_row, :) = 0;
setup_s = toc(stage_tic);

%% Solve with the variant's prior
stage_tic = tic;
[active_field, inversion_info] = ar_solve_advection_lean(observations, g, active_mask_3d, steps, prior_active, ...
    "prior", inv.prior, "advection", inv.advection, "coeff_horver", inv.coeff_horver, "tol", inv.tol, ...
    "maxit", inv.maxit, "store_class", "single", "advection_relax_exponent", inv.advection_relax_exponent, ...
    "nonnegative", inv.nonnegative, "projected_maxit", inv.projected_maxit, "row_weights", row_weights, ...
    prior_options{:});
solve_s = toc(stage_tic);
active_field = single(active_field);

%% Scores
predicted = observations.F(:, active_id) * double(active_field);
observed = full(observations.values);
used = observed ~= 0 & ~held_row;
trusted = used & observations.row_range_m >= 30000 & observations.row_range_m <= 120000;
holdout_scores = table();
if holdout
    for i_r = find(~keep_radar)'
        test = observed ~= 0 & observations.row_radar_id == i_r & observations.row_range_m >= 30000 & ...
            observations.row_range_m <= 120000;
        p = max(predicted(test), 0.1);
        o = max(observed(test), 0.1);
        holdout_scores = [holdout_scores; table(string(radar_selected.name(i_r)), nnz(test), ...
            norm(p - o) / norm(o), median(abs(log(p ./ o))), ...
            VariableNames=["radar" "n_rows" "relative_error" "median_abs_log_ratio"])]; %#ok<AGROW>
    end
    disp(holdout_scores)
end
solve_info = struct("night", night_name, "variant", variant, "holdout", holdout, "setup_s", setup_s, ...
    "solve_s", solve_s, "iterations", inversion_info.solver.iterations, "flag", inversion_info.solver.flag, ...
    "relres", inversion_info.solver.relres, "n_unknowns", numel(active_field), ...
    "misfit", norm(predicted(used) - observed(used)) / norm(observed(used)), ...
    "misfit_30_120km", norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted)), ...
    "holdout_scores", holdout_scores, "prior_options", {prior_options}, "solver", inversion_info.solver);
suffix = "";
if holdout
    suffix = "_holdout";
end
save(fullfile(c.data_dir, "inversion_" + variant + suffix + ".mat"), "active_field", "active_id", "solve_info", "-v7.3")
fprintf("[variant] %s %s%s: %d iterations (flag %d), solve %.0f s, setup %.0f s, misfit %.3f (30-120 km %.3f)\n", ...
    night_name, variant, suffix, solve_info.iterations, solve_info.flag, solve_s, setup_s, solve_info.misfit, ...
    solve_info.misfit_30_120km)
end
