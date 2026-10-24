%% Production test 4b: solver weight tests on a short window
% Solves a few consecutive frames of the case with several solver settings
% (term weights, non-negativity) and compares maps, a vertical section through one radar, negative
% densities and data misfit. Uses case.mat (P3); the density prior of the
% whole night is computed once and cached.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"))
if ~exist("case_name", "var")
    case_name = "us_east";
end
c = prod_case(case_name);
inv = c.inversion;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
column_mask = [];
if ~isempty(c.domain)
    column_mask = build_domain_mask(g, c.domain);
end
active_mask_3d = build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));

% Test window and section radar
frame_ids = 3:8;             % 00:30-03:00 UTC for us_east (radar discs strongest early)
section_radar = "KDGX";      % strongest beam imprint in the first us_east run
test_name = "radar_bias_tests";   % output figure and table name

%% Density prior (cached) -------------------------------------------------
prior_file = fullfile(c.data_dir, "prior_active.mat");
if isfile(prior_file)
    load(prior_file, "prior_active")
else
    prior_active = build_density_prior_fast(g, radar_scans, observations, ...
        "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
        "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, ...
        "maxit", inv.maxit, "active_id", active_id, "interpolation", inv.prior_interpolation, ...
        "length_scale_m", inv.prior_length_scale_m, "include_water", inv.prior_include_water);
    save(prior_file, "prior_active", "-v7.3")
end

%% Window observations and advection steps --------------------------------
window_obs = observations;
window_obs.values = observations.values(:, frame_ids);
window_times = c.times(frame_ids);
steps = build_radar_advection_steps(g, radar_scans(:, frame_ids), seconds(diff(window_times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=active_mask_3d);

%% Weight configurations --------------------------------------------------
configs = table( ...
    ["current"; "smoothness 1"; "smoothness 3"; "horizontal smoothing"; "advection 3"; "advection 3 + smoothness 1"], ...
    repmat(inv.prior, 6, 1), ...
    [inv.smoothness; 1; 3; 1; inv.smoothness; 1], ...
    [inv.advection; inv.advection; inv.advection; inv.advection; 3; 3], ...
    [inv.coeff_horver; inv.coeff_horver; inv.coeff_horver; 0.3; inv.coeff_horver; inv.coeff_horver], ...
    true(6, 1), ...
    repmat(inv.projected_maxit, 6, 1), ...
    VariableNames=["name" "prior" "smoothness" "advection" "coeff_horver" "nonnegative" "projected_maxit"]);

% A config_names list in the workspace runs only those configurations
if exist("config_names", "var")
    configs = configs(ismember(configs.name, config_names), :);
end

%% Solve each configuration ----------------------------------------------
F_active = observations.F(:, active_id);
results = cell(height(configs), 1);
metrics = table(strings(height(configs), 1), zeros(height(configs), 1), zeros(height(configs), 1), ...
    zeros(height(configs), 1), zeros(height(configs), 1), zeros(height(configs), 1), zeros(height(configs), 1), ...
    VariableNames=["name" "seconds" "iterations" "negative_share" "negative_mass" "data_misfit" "radar_bias"]);
bias_curves = cell(height(configs), 1);
for i_c = 1:height(configs)
    stage_tic = tic;
    [field, info] = solve_advection_lean(window_obs, g, active_mask_3d, steps, prior_active(:, frame_ids), ...
        prior=configs.prior(i_c), smoothness=configs.smoothness(i_c), advection=configs.advection(i_c), ...
        coeff_horver=configs.coeff_horver(i_c), tol=inv.tol, maxit=inv.maxit, store_class="single", ...
        advection_relax_exponent=inv.advection_relax_exponent, ...
        nonnegative=configs.nonnegative(i_c), projected_maxit=max(configs.projected_maxit(i_c), 1));
    % Relative data misfit on observed rows (zeros included)
    predicted = F_active * double(field);
    observed = full(window_obs.values);
    mask = observed ~= 0;
    misfit = norm(predicted(mask) - observed(mask)) / norm(observed(mask));
    % Column density against distance to the nearest radar (1 = no imprint)
    [bias_curves{i_c}, bias_edges_km, bias_index] = radar_distance_bias(field, active_id, g, radar_selected);
    metrics(i_c, :) = {configs.name(i_c), toc(stage_tic), info.solver.iterations, mean(field(:) < 0), ...
        -sum(field(field < 0)) / sum(field(field > 0)), misfit, mean(bias_index, "omitnan")};
    results{i_c} = field;
    stage_log(c, "P4b_" + strrep(configs.name(i_c), " ", "_"), stage_tic, sprintf("%d it, negative %.1f%%, misfit %.3f", ...
        info.solver.iterations, 100 * mean(field(:) < 0), misfit))
end
clear predicted observed F_active
disp(metrics)
if ~isfolder(c.figure_dir)
    mkdir(c.figure_dir)
end
writetable(metrics, fullfile(c.figure_dir, test_name + ".csv"))

%% Compare maps and the radar section -----------------------------------
mid = ceil(numel(frame_ids) / 2);
k = find(radar_selected.name == section_radar);
[~, section_y] = min(abs(g.lat - radar_selected.lat(k)));
zoom_lon = radar_selected.lon(k) + [-2.5 2.5];
zoom_lat = radar_selected.lat(k) + [-2 2];
land_color = [0.96 0.96 0.94];
water_color = [0.85 0.92 0.97];
background = permute(reshape(land_color, 1, 1, 3) .* ~g.mask_water + reshape(water_color, 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
fig = create_export_figure(3600, 900 * height(configs), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, height(configs), 3, TileSpacing="tight", Padding="tight");
for i_c = 1:height(configs)
    field = nan(g.sz, "single");
    field(active_id) = results{i_c}(:, mid);
    column = max(sum(field, 3, "omitnan"), 0) .* g.dz ./ 1000;
    column(~any(active_mask_3d, 3)) = nan;
    if i_c == 1
        c_lim = [0 quantile(column(:), 0.995)];
        section_lim = quantile(squeeze(field(:, section_y, :)), [0.002 0.998], "all");
        section_lim = [-1 1] * max(abs(section_lim));
    end
    ax = nexttile(t);
    plot_night_map(ax, g, background, column, c_lim, states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
    title(ax, configs.name(i_c) + sprintf(" (prior %.2g, smoothness %.2g, h %.2g)", ...
        configs.prior(i_c), configs.smoothness(i_c), configs.coeff_horver(i_c)), FontSize=12)
    ax = nexttile(t);
    plot_night_map(ax, g, background, column, c_lim, states, radar_selected, [zoom_lon zoom_lat], ones(g.sz(1:2)));
    yline(ax, g.lat(section_y), "w--")
    title(ax, section_radar + " " + string(window_times(mid), "HH:mm") + sprintf(" UTC, negative %.1f%%, misfit %.3f", ...
        100 * metrics.negative_share(i_c), metrics.data_misfit(i_c)), FontSize=12)
    ax = nexttile(t);
    section = squeeze(field(:, section_y, :))';
    imagesc(ax, g.lon, squeeze(g.z), section, AlphaData=~isnan(section))
    set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=11)
    xlim(ax, zoom_lon); ylim(ax, [0 3000]); clim(ax, section_lim)
    colormap(ax, interp1([0 0.5 1], [0.2 0.3 0.9; 1 1 1; 0.85 0.2 0.1], linspace(0, 1, 255)))
    cb = colorbar(ax); cb.Label.String = "Density (birds km^{-3})";
    title(ax, sprintf("Section at %.2f°N (signed)", g.lat(section_y)), FontSize=12)
end
export_fixed_figure(fig, fullfile(c.figure_dir, test_name + ".png"), Resolution=300)
close(fig)

%% Radar-distance curves ---------------------------------------------------
fig = create_export_figure(2400, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
ax = axes(fig); hold(ax, "on"); box(ax, "on")
bin_center = (bias_edges_km(1:end-1) + bias_edges_km(2:end)) / 2;
for i_c = 1:height(configs)
    plot(ax, bin_center, mean(bias_curves{i_c}, 2, "omitnan"), "-o", LineWidth=2, ...
        DisplayName=configs.name(i_c) + sprintf(" (bias %.3f, misfit %.3f)", metrics.radar_bias(i_c), metrics.data_misfit(i_c)))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean (normalized)")
legend(ax, Location="southwest"); set(ax, FontSize=13)
title(ax, "Radar imprint, mean over " + string(window_times(1), "HH:mm") + "-" + string(window_times(end), "HH:mm") + " UTC")
export_fixed_figure(fig, fullfile(c.figure_dir, test_name + "_curves.png"), Resolution=300)
close(fig)
