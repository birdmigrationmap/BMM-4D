%% Sandbox: do observation errors sigma(lambda) recover a known density better?
% Known truth with realistic structure: the sandbox solution (21 Oct 2026,
% 28 radars, 6 frames, unweighted) scaled to a busy (x1), medium (x0.1)
% and quiet (x0.02) night. Observations on the sandbox rows (same rows,
% times and operator F):
%   y = F rho_true x exp(e),  e ~ N(mu(lambda), s(lambda)^2),
% with lambda = (F rho_true) x counting volume (250 m gates) the expected
% birds per gate, mu(lambda) the near-radar deficit of the real sweep pairs
% (median ln upper / lower, SB_sigma_weights: -0.8 below 1 bird per gate,
% -0.1 at 5, 0 above) and s(lambda) the rest of their spread (at most 0.8,
% so the total RMS is about 1.1 below 1 bird per gate, as in the pairs). So the
% synthetic data carry the error measured on real data, and the truth is
% known. Each scenario is solved with perfect data (y = F rho_true: the
% inversion alone), unweighted, and with the weights (0.4 / sigma(lambda))^2
% of SB_sigma_weights (lambda from the radar prior of the synthetic data).
% The radar prior is rebuilt from each synthetic data set, as in a real run.
% Scores against the truth: column / true column by distance to the
% nearest radar, the near-radar column bias (within 30 km against 30-150
% km), relative L2 error of the 3D field and log RMSE of the columns; and
% solve time and PCG iterations.
% Last section: the same weights on the two real nights (sandbox, 12 May
% 2021), weighted against unweighted.
% Outputs: figures/sigma_synthetic.png, sigma_synthetic_maps.png,
% sigma_synthetic.csv, sigma_real.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
S = sandbox_load("sandbox");
g = S.g;
rng(1)

%% Truth and error model ----------------------------------------------------------
base = sandbox_run(S, struct("name", "synthetic truth (unweighted sandbox)", "video", false));
truth_base = max(double(base.field), 0);
scenarios = table(["busy" "medium" "quiet"]', [1 0.1 0.02]', VariableNames=["name" "scale"]);
% Deficit and spread of the real sweep pairs (SB_sigma_weights), in ln
deficit = @(lambda) -0.8 ./ (1 + (lambda / 1.2) .^ 1.5);
spread = @(lambda) min(sqrt(0.2 ^ 2 + 0.2 ^ 2 ./ max(lambda, 1e-6)), 0.8);   % total RMS ~1.1 below 1 bird per gate, as the pairs
sigma_model = @(lambda) sqrt(0.4 ^ 2 + 0.5 ^ 2 ./ max(lambda, 1e-6));   % as SB_sigma_weights
beam = bmm4d.radar.gate.beam();
volume_sample = beam.volume_sample_km3(double(S.observations.row_range_m));
F_active = S.observations.F(:, S.active_id);
observed_mask = S.observations.values ~= 0;
[row, col] = find(observed_mask);
n_time = size(observed_mask, 2);

%% Columns and distance to the nearest radar (for the scores) --------------------
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
air = bmm4d.grid.air_fraction(g, S.active_id);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for i_r = 1:height(S.radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r)))));
end
column_of = @(field) accumarray(column_id, mean(double(field) .* air, 2), [prod(g.sz(1:2)) 1]);
true_column_base = column_of(truth_base);
is_active = true_column_base > 0;
edges_km = 0:5:200;
centre_km = edges_km(1:end - 1)' + 2.5;
column_bin = discretize(distance_km(:), edges_km);

%% Solve every scenario ---------------------------------------------------------------
methods = {
    struct("name", "perfect data")
    struct("name", "unweighted")
    struct("name", "sigma(lambda)", "sigma_model", sigma_model)};
results = table();
curves = cell(height(scenarios), numel(methods));
columns = cell(height(scenarios), numel(methods));
lambda_share = zeros(height(scenarios), 1);
for i_s = 1:height(scenarios)
    truth = scenarios.scale(i_s) * truth_base;
    predicted_true = F_active * truth;
    lambda = predicted_true .* volume_sample;
    noise = deficit(lambda) + spread(lambda) .* randn(size(lambda));
    lambda_share(i_s) = mean(lambda(observed_mask) < 1);
    true_column = scenarios.scale(i_s) * true_column_base;
    for i_m = 1:numel(methods)
        % Synthetic observations on the sandbox rows; the prior is rebuilt from them
        S_syn = S;
        S_syn.prior_cache = containers.Map();
        index = sub2ind(size(observed_mask), row, col);
        value = predicted_true(index);
        if i_m > 1
            value = value .* exp(noise(index));
        end
        S_syn.observations.values = sparse(row, col, max(value, 1e-6), size(observed_mask, 1), n_time);
        option = methods{i_m};
        option.name = "synthetic " + scenarios.name(i_s) + ", " + option.name;
        option.video = false;
        wall = tic;
        result = sandbox_run(S_syn, option);
        wall_seconds = toc(wall);
        column = column_of(result.field);
        ok = is_active & ~isnan(column_bin);
        curves{i_s, i_m} = accumarray(column_bin(ok), column(ok), [numel(centre_km) 1], @sum, NaN) ./ ...
            accumarray(column_bin(ok), true_column(ok), [numel(centre_km) 1], @sum, NaN);
        columns{i_s, i_m} = column;
        near = ok & distance_km(:) < 30;
        mid = ok & distance_km(:) >= 30 & distance_km(:) < 150;
        near_bias = log2((sum(column(near)) / sum(true_column(near))) / (sum(column(mid)) / sum(true_column(mid))));
        field_error = norm(double(result.field(:)) - truth(:)) / norm(truth(:));
        column_log_rmse = sqrt(mean(log((max(column(ok), 0) + 0.1 * mean(true_column(ok))) ./ ...
            (true_column(ok) + 0.1 * mean(true_column(ok)))) .^ 2));
        results = [results; table(scenarios.name(i_s), string(methods{i_m}.name), round(lambda_share(i_s), 3), ...
            round(near_bias, 3), round(log2(sum(column(ok)) / sum(true_column(ok))), 3), round(field_error, 3), ...
            round(column_log_rmse, 3), round(result.seconds, 1), round(wall_seconds, 1), result.iterations, ...
            VariableNames=["scenario" "method" "share_lambda_below_1" "near_bias_log2" "total_bias_log2" ...
            "field_rel_l2" "column_log_rmse" "solve_seconds" "wall_seconds" "iterations"])]; %#ok<AGROW>
        disp(results(end, :))
    end
end
writetable(results, fullfile(figure_dir, "sigma_synthetic.csv"))

%% Plot results: by distance and scores ---------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 1900, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
title(t, "Known truth (sandbox field, scaled), data with the error of real sweep pairs: reconstructed / true column", FontSize=16)
colors = [0.5 0.5 0.5; 0.85 0.33 0.1; 0 0.45 0.74];
for i_s = 1:height(scenarios)
    ax = nexttile(t, i_s); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_m = 1:numel(methods)
        plot(ax, centre_km, curves{i_s, i_m}, "-", Color=colors(i_m, :), LineWidth=2.5, DisplayName=methods{i_m}.name)
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    xlim(ax, [0 200]); ylim(ax, [0.3 1.3]); set(ax, FontSize=13)
    xlabel(ax, "Distance to the nearest radar (km)")
    if i_s == 1, ylabel(ax, "Column / true column"), end
    title(ax, sprintf("%s (x%g), %.0f%% of rows with \\lambda < 1", scenarios.name(i_s), ...
        scenarios.scale(i_s), 100 * lambda_share(i_s)), FontWeight="bold", FontSize=14)
    legend(ax, Location="southeast", FontSize=12)
end
metrics = ["near_bias_log2" "column_log_rmse" "solve_seconds"];
metric_label = ["Near-radar bias (log_2)" "Column log RMSE" "Solve time (s)"];
for i_k = 1:numel(metrics)
    ax = nexttile(t, 3 + i_k); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    values = reshape(results.(metrics(i_k)), numel(methods), [])';
    b = bar(ax, categorical(scenarios.name, scenarios.name), values);
    for i_m = 1:numel(methods)
        b(i_m).FaceColor = colors(i_m, :);
    end
    yline(ax, 0, "k-", HandleVisibility="off")
    set(ax, FontSize=13); title(ax, metric_label(i_k), FontWeight="normal", FontSize=14)
    if i_k == 1, ylim(ax, [-0.35 0.02]), end   % colours as in the legends above
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sigma_synthetic.png"), Resolution=150)
close(fig)

%% Plot results: maps of the quiet scenario -------------------------------------------
% log2(reconstructed / true column), mean over frames: truth, unweighted, weighted
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
fig = bmm4d.plot.create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "Quiet scenario: log_2(reconstructed / true column), mean over frames", FontSize=16)
i_s = find(scenarios.name == "quiet");
true_column = scenarios.scale(i_s) * true_column_base;
for i_m = 1:numel(methods)
    ratio_map = nan(g.sz(1:2));
    ratio_map(is_active) = columns{i_s, i_m}(is_active) ./ true_column(is_active);
    ax = nexttile(t);
    bmm4d.plot.plot_night_map(ax, g, background, log2(max(ratio_map, 1 / 8)), [-1 1], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, signed_map)
    title(ax, methods{i_m}.name, FontWeight="bold", FontSize=14)
end
cb = colorbar(ax); cb.Layout.Tile = "east"; cb.Label.String = "log_2(column / true)"; cb.FontSize = 13;
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sigma_synthetic_maps.png"), Resolution=150)
close(fig)

%% Real nights: weighted against unweighted ----------------------------------------
% Column ratio weighted / unweighted (mean over frames) as a map and by
% distance to the nearest radar, and the expected birds per gate of the
% real observations within 30 km (from the radar prior).
case_names = ["sandbox" "night_20210512"];
fig = bmm4d.plot.create_export_figure(3600, 2200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
title(t, "Real nights: solution with sigma(\lambda) weights against the unweighted solution", FontSize=16)
real_summary = table();
real_curves = nan(numel(centre_km), numel(case_names));
lambda_edges_log10 = -3:0.1:3;
real_lambda_share = nan(numel(lambda_edges_log10) - 1, numel(case_names));
for i_c = 1:numel(case_names)
    S_real = sandbox_load(case_names(i_c));
    unweighted = sandbox_run(S_real, struct("name", "real " + case_names(i_c) + ", unweighted", "video", false));
    weighted = sandbox_run(S_real, struct("name", "real " + case_names(i_c) + ", sigma(lambda)", "sigma_model", sigma_model, "video", false));
    [ix, iy, ~] = ind2sub(S_real.g.sz, S_real.active_id);
    real_column_id = sub2ind(S_real.g.sz(1:2), ix, iy);
    real_air = bmm4d.grid.air_fraction(S_real.g, S_real.active_id);
    column_u = accumarray(real_column_id, mean(double(unweighted.field) .* real_air, 2), [prod(S_real.g.sz(1:2)) 1]);
    column_w = accumarray(real_column_id, mean(double(weighted.field) .* real_air, 2), [prod(S_real.g.sz(1:2)) 1]);
    active = column_u > 0;
    ratio_map = nan(S_real.g.sz(1:2));
    ratio_map(active) = column_w(active) ./ column_u(active);
    ax = nexttile(t, (i_c - 1) * 3 + 1);
    real_background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~S_real.g.mask_water + ...
        reshape([0.85 0.92 0.97], 1, 1, 3) .* S_real.g.mask_water, [2 1 3]);
    bmm4d.plot.plot_night_map(ax, S_real.g, real_background, log2(ratio_map), [-0.5 0.5], states, S_real.radar_selected, ...
        S_real.c.figure_bounds, ones(S_real.g.sz(1:2)));
    colormap(ax, signed_map); cb = colorbar(ax); cb.Label.String = "log_2(weighted / unweighted column)";
    title(ax, case_names(i_c), FontWeight="bold", FontSize=14, Interpreter="none")
    % By distance to the nearest radar
    [lon_r, lat_r] = ndgrid(S_real.g.lon, S_real.g.lat);
    real_distance = inf(S_real.g.sz(1:2));
    for i_r = 1:height(S_real.radar_selected)
        real_distance = min(real_distance, 6371 * acos(min(1, sind(lat_r) .* sind(S_real.radar_selected.lat(i_r)) + ...
            cosd(lat_r) .* cosd(S_real.radar_selected.lat(i_r)) .* cosd(lon_r - S_real.radar_selected.lon(i_r)))));
    end
    real_bin = discretize(real_distance(:), edges_km);
    ok = active & ~isnan(real_bin);
    curve = accumarray(real_bin(ok), column_w(ok), [numel(centre_km) 1], @sum, NaN) ./ ...
        accumarray(real_bin(ok), column_u(ok), [numel(centre_km) 1], @sum, NaN);
    real_curves(:, i_c) = curve;
    ax = nexttile(t, (i_c - 1) * 3 + 2); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, centre_km, curve, "-", Color=[0 0.45 0.74], LineWidth=2.5)
    yline(ax, 1, "k:"); xlim(ax, [0 200]); ylim(ax, [0.8 1.2]); set(ax, FontSize=13)
    xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Weighted / unweighted column")
    title(ax, sprintf("max change %.1f%%, radar-distance bias %.3f -> %.3f", 100 * max(abs(curve - 1)), ...
        unweighted.bias_200km, weighted.bias_200km), FontWeight="normal", FontSize=13)
    % Birds per gate of the real rows within 30 km
    prior_keys = S_real.prior_cache.keys;
    real_prior = S_real.prior_cache(prior_keys{1});
    real_predicted = max(S_real.observations.F(:, S_real.active_id) * double(real_prior), 0);
    real_lambda = real_predicted .* beam.volume_sample_km3(double(S_real.observations.row_range_m));
    near_rows = S_real.observations.values ~= 0 & S_real.observations.row_range_m < 30000;
    real_lambda_share(:, i_c) = histcounts(log10(real_lambda(near_rows)), lambda_edges_log10, Normalization="probability");
    ax = nexttile(t, (i_c - 1) * 3 + 3); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    histogram(ax, log10(real_lambda(near_rows)), -3:0.1:3, Normalization="probability", FaceColor=[0.5 0.5 0.5], EdgeColor="none")
    yyaxis(ax, "right"); plot(ax, -3:0.05:3, (0.4 ./ sigma_model(10 .^ (-3:0.05:3))) .^ 2, "-", Color=[0 0.45 0.74], LineWidth=2.5)
    ylabel(ax, "Weight (0.4 / \sigma)^2"); ax.YAxis(2).Color = [0 0.45 0.74];
    yyaxis(ax, "left"); ylabel(ax, "Share of rows within 30 km"); set(ax, FontSize=13)
    xlabel(ax, "log_{10} expected birds per gate \lambda")
    title(ax, sprintf("%.0f%% of rows within 30 km have \\lambda < 1", 100 * mean(real_lambda(near_rows) < 1)), FontWeight="normal", FontSize=13)
    real_summary = [real_summary; table(case_names(i_c), unweighted.bias_200km, weighted.bias_200km, unweighted.misfit_trusted, ...
        weighted.misfit_trusted, unweighted.seconds, weighted.seconds, max(abs(curve - 1)), ...
        VariableNames=["case" "bias_unweighted" "bias_weighted" "misfit_trusted_unweighted" "misfit_trusted_weighted" ...
        "seconds_unweighted" "seconds_weighted" "max_curve_change"])]; %#ok<AGROW>
end
disp(real_summary)
writetable(real_summary, fullfile(figure_dir, "sigma_real.csv"))
writetable(array2table([centre_km real_curves], VariableNames=["distance_km" "ratio_" + case_names]), ...
    fullfile(figure_dir, "sigma_real_curves.csv"))
writetable(array2table([lambda_edges_log10(1:end - 1)' + 0.05 real_lambda_share], VariableNames=["log10_lambda" "share_" + case_names]), ...
    fullfile(figure_dir, "sigma_real_lambda.csv"))
writetable(array2table([centre_km [curves{:}]], VariableNames=["distance_km" reshape(scenarios.name + "_" + ...
    ["perfect" "unweighted" "weighted"], 1, [])]), fullfile(figure_dir, "sigma_synthetic_curves.csv"))
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sigma_real.png"), Resolution=150)
close(fig)
