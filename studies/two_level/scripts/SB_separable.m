%% Sandbox: separable solver against the current solver
% Separable model (MODEL_REVIEW section 7, step 2): density = column amount
% x_1 times a vertical shape above the smoothed ground, plus local shape
% corrections x_2..x_K (solve_separable). Patterns are the us_east night's
% EOFs (SB_double_layer.m, shape_patterns.mat). Every option uses the
% logistic data weights (near range centred at 28 km, no-echo rows fading
% around 75 km), so only the field model differs. Scores, a table, the
% radar-distance curves, a map gallery and the projection of the current
% solution on K patterns go to studies/two_level/figures/sandbox.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "sandbox"), fullfile(project_root, "studies", "two_level", "functions"))
figure_dir = fullfile(project_root, "studies", "two_level", "figures", "sandbox");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end
patterns = load(fullfile(prod_case("us_east").data_dir, "shape_patterns.mat"));

%% Options ----------------------------------------------------------------
logistic28 = @(r, e, z) (1 ./ (1 + exp(-(r - 28000) / 3000))) .* (~z + z ./ (1 + exp((r - 75000) / 10000)));
sep = {"solver", "separable", "row_weight", logistic28, "smoothness", 1};
regional = {"regional_patterns", 5, "regional_length_m", 75000, "shape_length_m", 50000};
options = {
    struct("name", "current, prior 0.3", "row_weight", logistic28)
    struct("name", "current, prior 0.5", "row_weight", logistic28, "prior", 0.5)
    struct("name", "sep K3, mean shape", sep{1:4}, "smoothness", 0.3)
    struct("name", "sep K1, regional shape", sep{:}, "n_pattern", 1, "regional_patterns", 2)
    struct("name", "sep K1, regional, level 0.1", sep{:}, "n_pattern", 1, "regional_patterns", 2, "column_prior", 0.1)
    struct("name", "sep K1, regional, level 0.3", sep{:}, "n_pattern", 1, "regional_patterns", 2, "column_prior", 0.3)
    struct("name", "sep K3, regional, level 0.3", sep{:}, regional{:}, "n_pattern", 3, "shape_prior", 0.1, "column_prior", 0.3)
    struct("name", "sep K3, regional, level 0.12", sep{:}, regional{:}, "n_pattern", 3, "shape_prior", 0.03, "column_prior", 0.12)
    struct("name", "sep K6, regional, level 0.3", sep{:}, "n_pattern", 6, "regional_patterns", 2, "shape_prior", 0.03, "column_prior", 0.3)};
% Level prior: weight of ||W (x_1 - level)||^2, level = column of the radar
% prior. Per column it is about 6x a lean prior of the same weight, so
% level 0.12 is about the current prior 0.3 (a stand-in for the large-scale
% level V_L of step 3).

%% Run --------------------------------------------------------------------
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = two_level_run(S, options{i_o});
    fprintf("[sandbox] %-30s bias %.3f  misfit %.3f  trusted %.3f  detail %.2f  %4.0f s\n", results{i_o}.name, ...
        results{i_o}.bias_200km, results{i_o}.misfit, results{i_o}.misfit_trusted, results{i_o}.detail, results{i_o}.seconds)
end
summary = table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), ...
    cellfun(@(r) r.misfit, results), cellfun(@(r) r.misfit_trusted, results), cellfun(@(r) r.detail, results), ...
    cellfun(@(r) r.seconds, results), cellfun(@(r) r.iterations, results), ...
    VariableNames=["option" "bias_200km" "misfit" "misfit_trusted" "detail" "seconds" "iterations"]);
disp(summary)
writetable(summary, fullfile(figure_dir, "separable_comparison.csv"))

%% How many patterns does the current solution need? ------------------------
% Project the current solution (prior 0.3) column by column on the first K
% patterns (above the smoothed ground) and recompute its data misfit: the
% misfit the separable model would pay to hold the same field.
g = S.g;
[ix, iy, iz] = ind2sub(g.sz, S.active_id);
[column_xy, ~, column_of_cell] = unique(sub2ind(g.sz(1:2), ix, iy));
n_column = numel(column_xy);
n_active = numel(S.active_id);
ground_smooth = imgaussfilt(max(g.dem, 0), 30000 / mean([g.dx g.dy]));
z = squeeze(g.z);
height = z(iz) - ground_smooth(column_xy(column_of_cell));
phi = interp1(patterns.z_agl, patterns.shape_patterns, min(max(height, patterns.z_agl(1)), patterns.z_agl(end)));
observed = full(S.observations.values);
observed_mask = observed ~= 0;
trusted = observed_mask & S.observations.row_range_m >= 30000 & S.observations.row_range_m <= 120000;
F_active = S.observations.F(:, S.active_id);
current = double(results{1}.field);
projection = table((1:6)', nan(6, 1), nan(6, 1), VariableNames=["K" "misfit_trusted" "field_error"]);
for K = 1:6
    B = sparse(repmat((1:n_active)', K, 1), column_of_cell + n_column * (0:K - 1), reshape(phi(:, 1:K), [], 1), n_active, K * n_column);
    projected = max(B * ((B' * B + 1e-9 * speye(K * n_column)) \ (B' * current)), 0);
    predicted = F_active * projected;
    projection.misfit_trusted(K) = norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted));
    projection.field_error(K) = norm(projected - current, "fro") / norm(current, "fro");
end
disp(projection)
writetable(projection, fullfile(figure_dir, "separable_projection.csv"))

%% Compare the radar-distance curves ------------------------------------
fig = create_export_figure(3200, 1400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
ax = axes(fig); hold(ax, "on"); box(ax, "on")
center = (results{1}.edges_km(1:end-1) + results{1}.edges_km(2:end)) / 2;
colors = lines(numel(results));
for i_o = 1:numel(results)
    style = "-o";
    if options{i_o}.name.startsWith("current")
        style = "-s";
    end
    plot(ax, center, mean(results{i_o}.ratio, 2, "omitnan"), style, Color=colors(i_o, :), LineWidth=2, ...
        DisplayName=sprintf("%s (bias %.3f, trusted %.3f, detail %.2f)", results{i_o}.name, results{i_o}.bias_200km, results{i_o}.misfit_trusted, results{i_o}.detail))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean")
legend(ax, Location="eastoutside", FontSize=11); set(ax, FontSize=13); ylim(ax, [0.6 1.4])
title(ax, "Sandbox: separable solver against the current solver")
export_fixed_figure(fig, fullfile(figure_dir, "separable_curves.png"), Resolution=300)
close(fig)

%% Map gallery at the middle frame ------------------------------------------
t_mid = ceil(size(results{1}.field, 2) / 2);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
fig = create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, 3, TileSpacing="tight", Padding="tight");
for i_o = 1:numel(results)
    field = nan(g.sz, "single");
    field(S.active_id) = results{i_o}.field(:, t_mid);
    column = sum(field, 3, "omitnan") .* g.dz ./ 1000;
    column(~any(S.active_mask_3d, 3)) = nan;
    ax = nexttile(t);
    plot_night_map(ax, g, background, column, [0 200], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    title(ax, {results{i_o}.name, sprintf("bias %.3f, trusted %.3f, detail %.2f", results{i_o}.bias_200km, ...
        results{i_o}.misfit_trusted, results{i_o}.detail)}, FontSize=12, FontWeight="normal")
end
cb = colorbar(ax); cb.Layout.Tile = "east"; cb.Label.String = "Column density (birds km^{-2}), " + string(S.c.times(t_mid), "HH:mm") + " UTC";
export_fixed_figure(fig, fullfile(figure_dir, "gallery_separable.png"), Resolution=150)
close(fig)
