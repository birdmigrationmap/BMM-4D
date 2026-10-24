%% Sandbox: far-range deficit in the two-level solver
% Columns 150-250 km from the nearest radar have no data at all (no beam
% reaches them), so they are the level V_L plus whatever the deviation
% and advection carry there. Two causes were found:
% 1. the relative amplitude of the fine deviation weights columns by
%    1/level^2, so V_L sits below the dense patches the 600 km level
%    cannot follow (deviation +12% on observed columns), and the far
%    columns, which fall back to V_L, end up low;
% 2. advection of the whole density keeps earlier, lower densities in the
%    data-free areas while the level rises through the night (deviation
%    there 0 at the first frame, about -5 birds km^-2 at the last).
% Fixes (now the solve_two_level defaults): unweighted amplitude of the
% deviation (FINE_AMPLITUDE "uniform") and advection of the deviation
% only (ADVECTION_ON "deviation"). Outputs in prod/experiments/two_level/figures/sandbox.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"), fullfile(project_root, "prod", "sandbox"), fullfile(project_root, "prod", "experiments", "two_level", "functions"))
figure_dir = fullfile(project_root, "prod", "experiments", "two_level", "figures", "sandbox");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end

%% Options ----------------------------------------------------------------
base = {"fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.03};
two_level = @(name, extra, advection) struct("name", name, "solver", "two_level", "n_pattern", 3, ...
    "regional_patterns", 2, "advection", advection, "panel", false, "two_level", struct(base{:}, extra{:}));
options = {
    two_level("before: relative amplitude, advect density", {"fine_amplitude", "relative", "advection_on", "density"}, 1)
    two_level("uniform weights, advect density", {"weight_reference", "uniform", "advection_on", "density"}, 1)
    two_level("uniform amplitude, advect density", {"advection_on", "density"}, 1)
    two_level("uniform amplitude, no advection", {}, 0)
    two_level("uniform amplitude, advect deviation (default)", {}, 1)};

%% Run and diagnose ---------------------------------------------------------
% Data support of each column: total beam weight of the observed rows
g = S.g;
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
[~, ~, column_of_cell] = unique(sub2ind(g.sz(1:2), ix, iy));
observed_row = any(S.observations.values ~= 0, 2);
is_supported = accumarray(column_of_cell, full(sum(S.observations.F(observed_row, S.active_id), 1))') > 0;
results = cell(numel(options), 1);
deviation_by_time = zeros(numel(options), size(S.observations.values, 2));
mean_deviation_observed = zeros(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = two_level_run(S, options{i_o});
    deviation = double(results{i_o}.info.column_amount) - double(results{i_o}.info.column_level);
    deviation_by_time(i_o, :) = mean(deviation(~is_supported, :), 1);
    mean_deviation_observed(i_o) = mean(deviation(is_supported, :), "all");
    fprintf("[far range] %-46s bias %.3f  trusted %.3f  detail %.2f  %4.0f s\n", results{i_o}.name, ...
        results{i_o}.bias_200km, results{i_o}.misfit_trusted, results{i_o}.detail, results{i_o}.seconds)
end
curves = cell2mat(cellfun(@(r) mean(r.ratio, 2, "omitnan")', results, UniformOutput=false));
summary = table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), cellfun(@(r) r.misfit_trusted, results), ...
    cellfun(@(r) r.detail, results), mean_deviation_observed, deviation_by_time(:, end), curves(:, 7), curves(:, 8), ...
    VariableNames=["option" "bias_200km" "misfit_trusted" "detail" "deviation_observed" "deviation_unobserved_last" ...
    "ratio_150_200km" "ratio_200_300km"]);
disp(summary)
writetable(summary, fullfile(figure_dir, "far_range_comparison.csv"))

%% Plot results -----------------------------------------------------------------
fig = create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
colors = lines(numel(results));
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
center = (results{1}.edges_km(1:end - 1) + results{1}.edges_km(2:end)) / 2;
for i_o = 1:numel(results)
    plot(ax, center, curves(i_o, :), "-o", LineWidth=2, Color=colors(i_o, :), ...
        DisplayName=sprintf("%s (bias %.3f)", results{i_o}.name, results{i_o}.bias_200km))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean"); ylim(ax, [0.8 1.2])
legend(ax, Location="southwest", FontSize=11); set(ax, FontSize=13)
title(ax, "Radar-distance curves", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_o = 1:numel(results)
    plot(ax, S.c.times, deviation_by_time(i_o, :), "-o", LineWidth=2, Color=colors(i_o, :), DisplayName=results{i_o}.name)
end
yline(ax, 0, "k--", HandleVisibility="off")
ylabel(ax, "Mean x_1 - V_L on data-free columns (birds km^{-2})"); set(ax, FontSize=13)
title(ax, "Deviation where no beam reaches", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "far_range.png"), Resolution=150)
close(fig)
