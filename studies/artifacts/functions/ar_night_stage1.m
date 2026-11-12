function fit = ar_night_stage1(night_name)
%AR_NIGHT_STAGE1 Observations and prior-range fit of one night (dev, cluster job).
%   FIT = AR_NIGHT_STAGE1(NIGHT_NAME) builds the observations of the
%   us_east night NIGHT_NAME (ar_night_config) when case.mat is missing,
%   builds the production profile field (as bmm4d.invert, lean) and fits
%   the horizontal and vertical ranges of the anomaly from it for the
%   correlated (GP) and the one-term Matern prior (ar_fit_ranges). Saves
%   range_fit.mat in the night's data folder.

project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
c = ar_night_config(night_name);
inv = c.inversion;
if ~isfile(fullfile(c.data_dir, "case.mat"))
    stage_tic = tic;
    bmm4d.build_observations(c);
    fprintf("[stage1] %s observations built in %.0f s\n", night_name, toc(stage_tic))
end
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans")
column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));
prior_active = bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id, "ground_smoothing_m", inv.ground_smoothing_m);
fit_tic = tic;
fit = ar_fit_ranges(observations, g, prior_active, active_id);
fit.seconds = toc(fit_tic);
fit.night = night_name;
save(fullfile(c.data_dir, "range_fit.mat"), "fit")
end
