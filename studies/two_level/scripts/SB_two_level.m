%% Sandbox: two-level separable solver (prototype step 3)
% Column amount x_1 = P v + e with the large-scale level V_L (v) and the
% shape coefficients c_k on a 40 km grid, the fine deviation e relative to
% the level (solve_two_level). No radar-profile prior: the level is
% estimated from all radars at once; the radar profiles only start the
% shape (regional_shape.m) and centre its prior. Compared with the current
% solver and the best step-2 run (level = radar prior column). Every option
% uses the logistic data weights. Outputs in studies/two_level/figures/sandbox.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "sandbox"), fullfile(project_root, "studies", "two_level", "functions"))
figure_dir = fullfile(project_root, "studies", "two_level", "figures", "sandbox");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end

%% Options ----------------------------------------------------------------
logistic28 = @(r, e, z) (1 ./ (1 + exp(-(r - 28000) / 3000))) .* (~z + z ./ (1 + exp((r - 75000) / 10000)));
two_level = {"solver", "two_level", "row_weight", logistic28};
options = {
    struct("name", "current, prior 0.3", "row_weight", logistic28)
    struct("name", "current, prior 0.5", "row_weight", logistic28, "prior", 0.5)
    struct("name", "step 2: sep K3, radar level 0.3", "solver", "separable", "row_weight", logistic28, "smoothness", 1, ...
        "n_pattern", 3, "regional_patterns", 5, "regional_length_m", 75000, "shape_prior", 0.1, "shape_length_m", 50000, "column_prior", 0.3)
    struct("name", "2L K3, mean-shape start, fine 0.1", two_level{:}, "n_pattern", 3, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.1, "shape_prior", 0.1, "shape_length_m", 150000))
    struct("name", "2L K3, fine 0.2", two_level{:}, "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000))
    struct("name", "2L K3, fine 0.1", two_level{:}, "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.1, "shape_prior", 0.03, "shape_length_m", 150000))
    struct("name", "2L K6, fine 0.2", two_level{:}, "n_pattern", 6, "regional_patterns", 5, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000))
    struct("name", "2L K6, fine 0.1", two_level{:}, "n_pattern", 6, "regional_patterns", 5, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.1, "shape_prior", 0.03, "shape_length_m", 150000))
    struct("name", "2L K3, fine 0.2, level 0.03", two_level{:}, "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.03))
    struct("name", "2L K3, fine 0.2, level 0.01", two_level{:}, "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.01))
    struct("name", "2L K6, fine 0.2, level 0.03", two_level{:}, "n_pattern", 6, "regional_patterns", 5, ...
        "two_level", struct("near_range", "off", "fine_amplitude", "relative", "advection_on", "density", "fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.03))};
% near_range "off", relative amplitude and advection of the whole density:
% these runs predate steps 4 and the far-range fix.
% Defaults of solve_two_level otherwise: nodes 40 km, level GMRF over
% 600 km (weight 1, which leaves V_L almost constant: the "level" runs
% loosen it), fine deviation correlated over 30 km, 3 rounds.

%% Run --------------------------------------------------------------------
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = two_level_run(S, options{i_o});
    fprintf("[sandbox] %-34s bias %.3f  misfit %.3f  trusted %.3f  detail %.2f  %4.0f s\n", results{i_o}.name, ...
        results{i_o}.bias_200km, results{i_o}.misfit, results{i_o}.misfit_trusted, results{i_o}.detail, results{i_o}.seconds)
end
summary = table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), ...
    cellfun(@(r) r.misfit, results), cellfun(@(r) r.misfit_trusted, results), cellfun(@(r) r.detail, results), ...
    cellfun(@(r) r.seconds, results), VariableNames=["option" "bias_200km" "misfit" "misfit_trusted" "detail" "seconds"]);
disp(summary)
writetable(summary, fullfile(figure_dir, "two_level_comparison.csv"))

%% Trade-off: radar bias against trusted misfit ---------------------------------
fig = create_export_figure(2400, 1600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
ax = axes(fig); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
groups = ["current" "step 2" "2L K3" "2L K6"];
colors = [0.4 0.4 0.4; 0.85 0.33 0.1; 0 0.45 0.74; 0.47 0.67 0.19];
for i_g = 1:numel(groups)
    in_group = cellfun(@(r) startsWith(r.name, groups(i_g)), results);
    plot(ax, summary.misfit_trusted(in_group), summary.bias_200km(in_group), "o", MarkerSize=12, LineWidth=2, ...
        Color=colors(i_g, :), MarkerFaceColor=colors(i_g, :), DisplayName=groups(i_g))
    text(ax, summary.misfit_trusted(in_group) + 0.002, summary.bias_200km(in_group), summary.option(in_group), FontSize=10)
end
xlabel(ax, "Trusted misfit (30-120 km)"); ylabel(ax, "Radar-distance bias (within 200 km)")
legend(ax, Location="northwest"); set(ax, FontSize=13); xlim(ax, [0.22 0.36])
title(ax, "Sandbox: two-level solver against the current solver")
export_fixed_figure(fig, fullfile(figure_dir, "two_level_tradeoff.png"), Resolution=300)
close(fig)

%% Level, deviation and shape maps of one two-level run -----------------------------
% Middle frame of "2L K3, fine 0.2, level 0.03": column amount, level V_L, relative
% deviation x_1 / V_L - 1, and the two shape coefficients (c = 0 is the
% night's mean shape; c_2 > 0 puts more birds aloft, see double_layer_eofs.png).
r = results{find(cellfun(@(r) r.name == "2L K3, fine 0.2, level 0.03", results), 1)};
g = S.g;
t_mid = ceil(size(r.field, 2) / 2);
to_grid = @(values) reshape(accumarray(double(r.info.column_xy), double(values), [prod(g.sz(1:2)) 1], [], nan), g.sz(1:2));
amount = to_grid(r.info.column_amount(:, t_mid));
level = to_grid(r.info.column_level(:, t_mid));
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
fig = create_export_figure(3600, 2200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
panels = {amount, [0 200], "Column amount x_1 (birds km^{-2})", parula(256)
    level, [0 200], "Level V_L (birds km^{-2})", parula(256)
    amount ./ max(level, 10) - 1, [-1 1], "Relative deviation x_1 / V_L - 1", signed_map};
for k = 2:3
    shape_node = double(r.info.shape_nodes(:, k - 1, t_mid));
    shape_grid = griddata(r.info.node_lon, r.info.node_lat, shape_node, lon_grid, lat_grid, "linear");
    shape_grid(isnan(amount)) = nan;
    panels(end + 1, :) = {shape_grid, [-1 1], sprintf("Shape coefficient c_%d", k), signed_map}; %#ok<SAGROW>
end
for i_p = 1:size(panels, 1)
    ax = nexttile(t);
    plot_night_map(ax, g, background, panels{i_p, 1}, panels{i_p, 2}, states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, panels{i_p, 4}); colorbar(ax)
    title(ax, panels{i_p, 3}, FontSize=13, FontWeight="normal")
end
% Mean vertical shape with the domain-mean coefficients, against phi_1
ax = nexttile(t); hold(ax, "on"); box(ax, "on")
patterns = load(fullfile(prod_case("us_east").data_dir, "shape_patterns.mat"));
mean_c = mean(double(r.info.shape_nodes(:, :, t_mid)), 1);
plot(ax, patterns.shape_patterns(:, 1), patterns.z_agl / 1000, "k--", LineWidth=2, DisplayName="night mean \phi_1")
plot(ax, patterns.shape_patterns(:, 1:3) * [1 mean_c]', patterns.z_agl / 1000, LineWidth=2.5, DisplayName="sandbox mean shape")
xlabel(ax, "Shape (km^{-1})"); ylabel(ax, "Height above smoothed ground (km)"); legend(ax); set(ax, FontSize=13)
title(ax, "Domain-mean shape, " + string(S.c.times(t_mid), "HH:mm") + " UTC", FontSize=13, FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "two_level_fields.png"), Resolution=150)
close(fig)
