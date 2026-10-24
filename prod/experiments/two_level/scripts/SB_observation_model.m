%% Sandbox: observation model in the two-level solver (prototype step 4)
% Replaces the data cuts by an observation model in solve_two_level, each
% part switched on separately: near-range recovery a(r) (fixed, estimated
% globally, or estimated per elevation group from held-out near rows),
% no-echo gates as one-sided penalties below a detection limit, and a
% multiplicative error by reweighting. All observation-model runs use
% every row (no row weights); the logistic-weight run is the step-3
% reference. A held-out diagnostic (solve without rows < 40 km, compare
% their observations with the prediction) shows the near-range deficit by
% range and elevation. Outputs in prod/experiments/two_level/figures/sandbox.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"), fullfile(project_root, "prod", "sandbox"), fullfile(project_root, "prod", "experiments", "two_level", "functions"))
figure_dir = fullfile(project_root, "prod", "experiments", "two_level", "figures", "sandbox");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end

%% Options ----------------------------------------------------------------
% Two-level settings from step 3: K = 3, radar regional shape, fine 0.2,
% level 0.03.
logistic28 = @(r, e, z) (1 ./ (1 + exp(-(r - 28000) / 3000))) .* (~z + z ./ (1 + exp((r - 75000) / 10000)));
held_out_40 = @(r, e, z) double(r >= 40000);
base = {"fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.03, ...
    "fine_amplitude", "relative", "advection_on", "density"};   % settings before the far-range fix
near_off = {"near_range", "off"};
near_elevation = {"near_range", "estimate", "near_range_elev_edges", [0 1.2 2 3 90], "near_range_fit_max_m", 50000};
two_level = @(name, row_weight, extra) struct("name", name, "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
    "row_weight", row_weight, "two_level", struct(base{:}, extra{:}));
options = {
    two_level("logistic weights (step 3)", logistic28, near_off)
    two_level("no cuts, no model", [], near_off)
    two_level("held out < 40 km", held_out_40, near_off)
    two_level("near fixed (10 km, k 3.8)", [], {"near_range", "fixed"})
    two_level("near estimated, one curve", [], {"near_range", "estimate", "near_range_elev_edges", [0 90], "near_range_fit_max_m", 40000})
    two_level("near estimated by elevation", [], near_elevation)
    two_level("no-echo one-sided", [], [near_off {"no_echo", "one_sided"}])
    two_level("error cv 0.4", [], [near_off {"error_cv", 0.4}])
    two_level("near elev + no-echo", [], [near_elevation {"no_echo", "one_sided"}])
    two_level("near elev + error cv 0.4", [], [near_elevation {"error_cv", 0.4}])
    two_level("near elev + error, floor 15", [], [near_elevation {"error_cv", 0.4, "error_floor", 15}])
    two_level("all three", [], [near_elevation {"no_echo", "one_sided", "error_cv", 0.4}])
    struct("name", "current, logistic, prior 0.3", "row_weight", logistic28)};

%% Run --------------------------------------------------------------------
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = two_level_run(S, options{i_o});
    fprintf("[sandbox] %-30s bias %.3f  misfit %.3f  trusted %.3f  detail %.2f  %4.0f s\n", results{i_o}.name, ...
        results{i_o}.bias_200km, results{i_o}.misfit, results{i_o}.misfit_trusted, results{i_o}.detail, results{i_o}.seconds)
end
summary = table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), ...
    cellfun(@(r) r.misfit, results), cellfun(@(r) r.misfit_trusted, results), cellfun(@(r) r.detail, results), ...
    cellfun(@(r) r.seconds, results), VariableNames=["option" "bias_200km" "misfit" "misfit_trusted" "detail" "seconds"]);
disp(summary)
writetable(summary, fullfile(figure_dir, "observation_model_comparison.csv"))
near_curve = results{6}.info.near_curve;
disp(near_curve)
writetable(near_curve, fullfile(figure_dir, "observation_model_near_curve.csv"))

%% Held-out near-range ratio by range and elevation ------------------------------
% Observed / predicted on echo rows within 60 km, prediction from the run
% without rows < 40 km (so it comes from the farther data and the priors).
held_out = results{3};
observed = full(S.observations.values);
predicted = S.observations.F(:, S.active_id) * double(held_out.field);
range_km = double(S.observations.row_range_m) / 1000;
elevation = double(S.observations.row_elev_deg);
range_edges = 0:5:60;
elev_edges = [0 1.2 2 3 90];
ratio = nan(numel(elev_edges) - 1, numel(range_edges) - 1);
for i_e = 1:numel(elev_edges) - 1
    for i_r = 1:numel(range_edges) - 1
        in_bin = range_km >= range_edges(i_r) & range_km < range_edges(i_r + 1) & elevation >= elev_edges(i_e) & elevation < elev_edges(i_e + 1);
        selected = repmat(in_bin, 1, size(observed, 2)) & observed > 1e-6 & predicted > 1;
        ratio(i_e, i_r) = sum(observed(selected)) / sum(predicted(selected));
    end
end
disp(array2table(round(ratio, 2), VariableNames=compose("%d-%d km", range_edges(1:end - 1)', range_edges(2:end)'), ...
    RowNames=compose("%.1f-%.1f deg", elev_edges(1:end - 1)', elev_edges(2:end)')))

fig = create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
colors = parula(numel(elev_edges));
center = (range_edges(1:end - 1) + range_edges(2:end)) / 2;
r_line = linspace(0, 60, 200);
for i_e = 1:numel(elev_edges) - 1
    label = sprintf("%.1f-%.1f°", elev_edges(i_e), elev_edges(i_e + 1));
    plot(ax, center, ratio(i_e, :), "o", MarkerSize=9, LineWidth=2, Color=colors(i_e, :), DisplayName=label + " held out")
    plot(ax, r_line, 1 - exp(-(r_line * 1000 / near_curve.r0_m(i_e)) .^ near_curve.k(i_e)), "-", LineWidth=2, ...
        Color=colors(i_e, :), DisplayName=sprintf("a(r): r_0 %.1f km, k %.2f", near_curve.r0_m(i_e) / 1000, near_curve.k(i_e)))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Range (km)"); ylabel(ax, "Observed / predicted (held out)"); ylim(ax, [0 1.5])
legend(ax, Location="eastoutside", FontSize=11); set(ax, FontSize=13)
title(ax, "Near-range recovery by elevation", FontWeight="normal")
% Radar-distance curves of the main runs
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
shown = ["no cuts, no model" "logistic weights (step 3)" "near estimated, one curve" "near estimated by elevation" "all three" "current, logistic, prior 0.3"];
colors = lines(numel(shown));
bin_center = (results{1}.edges_km(1:end - 1) + results{1}.edges_km(2:end)) / 2;
for i_s = 1:numel(shown)
    r = results{cellfun(@(x) x.name == shown(i_s), results)};
    plot(ax, bin_center, mean(r.ratio, 2, "omitnan"), "-o", LineWidth=2, Color=colors(i_s, :), ...
        DisplayName=sprintf("%s (bias %.3f, trusted %.3f)", r.name, r.bias_200km, r.misfit_trusted))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean"); ylim(ax, [0.6 1.3])
legend(ax, Location="southwest", FontSize=11); set(ax, FontSize=13)
title(ax, "Radar-distance curves", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "observation_model.png"), Resolution=150)
close(fig)
