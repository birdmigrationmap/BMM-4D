%% Near radar D: drop rules in the solve, real and synthetic (dev, ~50 min)
% Which rows to drop near the radars, and can the data around fill the
% gap? Each rule removes observation entries before the prior and the
% weights are rebuilt (nr_apply_drop), then solves with the production
% settings (ar_solve, plain PCG; ~110 s per rule on ar_midwest_shared
% including the prior). Run NR_B first (deficit shape and adaptive R_min
% per radar-night, with a 5% tolerance).
% Real cases:
%   ar_midwest_shared - 28 Midwest radars, 21-22 Oct 2026, 7 frames
%                       04:30-07:30 UTC (us_east cleaning; VCP 35/212/215
%                       and staggered cuts at KMPX, KICT);
%   night_20210512    - the same radars, 13 May 2021, 16 frames 02:30-10:00
%                       UTC (older estimator cleaning, no sd: equal
%                       weights; VCP 31/32/35, where elevation and cut type
%                       differ).
% Validation without a threshold: the split-cut entries within 20 km of
% every second radar are held out of every solve (also "all rows"). They
% see the air below ~500 m near the radar, where the batch cuts at 1.8 deg
% overlap them, so observed / predicted on them is an independent check
% of the near-radar field. The rules are scored on:
%   - column against distance to the radar (relative to 60-100 km);
%   - observed / predicted of the held-out split entries;
%   - observed / predicted of the dropped entries (what was removed) and
%     of the kept batch entries by ring (what was left);
%   - misfit of the kept entries at 30-120 km, solve time.
% Synthetic (ar_midwest_shared): truth = the "batch < 20 km" field (T1),
% T1 with the measured batch deficit put on the batch entries (pooled
% NR_B shape), and T1 times a smooth random field per frame (T2: local
% structure the advection cannot carry). Observations F x truth on the
% observed entries, no noise. Column error against the truth by distance.
% Outputs: <local_data_dir>/near_radar/NR_D_results.mat, figures/NR_D_*.csv,
% figures/NR_D_drop_solves.png.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"), ...
    fullfile(project_root, "studies", "near_radar", "functions"))
figure_dir = fullfile(project_root, "studies", "near_radar", "figures");
out_dir = fullfile(bmm4d.util.local_data_dir(), "near_radar");
radar_night = readtable(fullfile(figure_dir, "NR_B_radar_night.csv"), TextType="string");
pooled = readtable(fullfile(figure_dir, "NR_B_pooled.csv"), TextType="string");

%% Rules
% max_range_km NaN: per radar, from NR_B (adaptive)
rules = struct( ...
    "name", {"all rows", "elevation < 20 km", "batch < 10 km", "batch < 20 km", "batch < 30 km", ...
    "batch+staggered < 20 km", "batch adaptive", "batch all ranges"}, ...
    "cuts", {"none", "elevation", "batch", "batch", "batch", "batch_staggered", "batch", "batch"}, ...
    "max_range_km", {0, 20, 10, 20, 30, 20, NaN, Inf});

%% Real cases
% The spring night (16 frames, ~2x slower) gets the rules that differ there
cases = struct("name", {"ar_midwest_shared", "night_20210512"}, "set", {"us_east_20261021", "night_20210512"}, ...
    "rules", {1:8, [1 2 4 5 7]});
summary = table();
distance = table();
residuals = table();
radars = table();
for i_case = 1:numel(cases)
    case_name = cases(i_case).name;
    if case_name == "ar_midwest_shared"
        c = ar_test_case(case_name, [-100 -85 36 46], ...
            datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30"), "read_options", {"shape", false});
    else
        c = case_settings(case_name);
    end
    S = ar_load(c);
    manifest = load(fullfile(c.data_dir, "manifest.mat"), "radars", "raw_files");
    cut_code = nr_cut_types(S.observations, S.radar_selected, manifest);
    observed = S.observations.values ~= 0;
    fprintf("%s: observed entries by cut code 0-5: %s\n", case_name, mat2str(histcounts(cut_code(observed), -0.5:1:5.5)))

    % Held-out split entries within 20 km of every second radar
    range_km = double(S.observations.row_range_m) / 1000;
    holdout = observed & cut_code == 1 & range_km < 20 & mod(double(S.observations.row_radar_id), 2) == 0;

    % Adaptive R_min of each radar on this night (missing: the night's median)
    mine = radar_night(radar_night.set == cases(i_case).set & isfinite(radar_night.r_min_km), :);
    [found, id] = ismember(S.radar_selected.name, mine.station);
    adaptive_km = repmat(median(mine.r_min_km), height(S.radar_selected), 1);
    adaptive_km(found) = mine.r_min_km(id(found));

    for rule = rules(cases(i_case).rules)
        if isnan(rule.max_range_km)
            rule.max_range_km = adaptive_km;
        end
        drop = nr_drop_rows(S.observations, cut_code, rule);
        S_rule = nr_apply_drop(S, drop | holdout);
        result = ar_solve(S_rule, case_name + ": " + rule.name);
        clear S_rule
        score = nr_score(S, result.field, cut_code, drop);

        % Held-out split entries
        predicted = S.observations.F(:, S.active_id) * double(result.field);
        values = full(S.observations.values);
        ok = holdout & values > 0.5 & predicted > 0.5;
        holdout_ratio = sum(values(ok)) / sum(predicted(ok));
        dropped = score.residuals(score.residuals.dropped, :);
        summary = [summary; table(string(case_name), string(rule.name), nnz(drop) / nnz(observed), ...
            sum(dropped.observed) / sum(dropped.predicted), holdout_ratio, score.misfit_trusted, result.seconds, ...
            score.distance.column_ratio(1), score.distance.column_ratio(2), score.distance.column_ratio(3), ...
            VariableNames=["case" "rule" "dropped_share" "dropped_ratio" "holdout_ratio" "misfit_trusted" "seconds" ...
            "column_0_10" "column_10_20" "column_20_30"])]; %#ok<AGROW>
        tag = @(t) [table(repmat(string(case_name), height(t), 1), repmat(string(rule.name), height(t), 1), ...
            VariableNames=["case" "rule"]) t];
        distance = [distance; tag(score.distance)]; %#ok<AGROW>
        residuals = [residuals; tag(score.residuals)]; %#ok<AGROW>
        radars = [radars; tag(score.radars)]; %#ok<AGROW>
        disp(summary(end, :))
        if case_name == "ar_midwest_shared" && rule.name == "batch < 20 km"
            truth = max(result.field, 0);
        end
        clear result predicted values
    end
    if case_name == "ar_midwest_shared"
        S_midwest = S;
        cut_code_midwest = cut_code;
    end
    clear S cut_code
end
writetable(summary, fullfile(figure_dir, "NR_D_summary.csv"))
writetable(distance, fullfile(figure_dir, "NR_D_distance.csv"))
writetable(residuals, fullfile(figure_dir, "NR_D_residuals.csv"))
writetable(radars, fullfile(figure_dir, "NR_D_radars.csv"))

%% Synthetic: can the data around fill the gap?
S = S_midwest;
cut_code = cut_code_midwest;
clear S_midwest cut_code_midwest
observed = S.observations.values ~= 0;
F = S.observations.F(:, S.active_id);
synthetic_values = @(field) sparse(double(observed) .* (F * double(field) + eps));

% T2: T1 times exp(0.3 z), z smooth (20 km horizontal, 300 m vertical), one draw per frame
rng(1)
dx_km = mean(diff(S.g.x)) / 1000;
truth_2 = truth;
for t = 1:size(truth, 2)
    z = randn(S.g.sz);
    z = smoothdata(smoothdata(z, 1, "gaussian", 5 * 20 / dx_km), 2, "gaussian", 5 * 20 / dx_km);
    z = smoothdata(z, 3, "gaussian", 5 * 300 / S.g.dz);
    z = z(S.active_id) / std(z(S.active_id));
    truth_2(:, t) = truth(:, t) .* exp(0.3 * z);
end

% Measured batch deficit (NR_B pooled, all heights) on the batch entries by ring; 1 beyond 25 km
batch_pooled = groupsummary(pooled(pooled.waveform == "batch", :), "range_lo_km", "sum", ["sum_density" "sum_predicted"]);
far = batch_pooled.range_lo_km >= 25 & batch_pooled.range_lo_km < 50;
far_level = sum(batch_pooled.sum_sum_density(far)) / sum(batch_pooled.sum_sum_predicted(far));
deficit = min(batch_pooled.sum_sum_density ./ batch_pooled.sum_sum_predicted / far_level, 1);
deficit(batch_pooled.range_lo_km >= 25) = 1;
range_km = repmat(double(S.observations.row_range_m) / 1000, 1, size(observed, 2));
entry_deficit = ones(size(observed));
is_batch = cut_code == 4 & range_km < 25;
entry_deficit(is_batch) = interp1(batch_pooled.range_lo_km, deficit, range_km(is_batch), "previous", 1);

truths = struct("name", {"T1", "T1 + batch deficit", "T2 (local structure)"}, ...
    "field", {truth, truth, truth_2}, "values", {synthetic_values(truth), ...
    synthetic_values(truth) .* entry_deficit, synthetic_values(truth_2)});
synthetic_rules = rules([1 4 5]);
synthetic = table();
for truth_case = truths
    for rule = synthetic_rules
        drop = nr_drop_rows(S.observations, cut_code, rule);
        S_rule = nr_apply_drop(S, drop, truth_case.values);
        result = ar_solve(S_rule, truth_case.name + ": " + rule.name);
        score = nr_score(S_rule, result.field, cut_code, drop, truth_case.field);
        clear S_rule
        synthetic = [synthetic; [table(repmat(string(truth_case.name), height(score.distance), 1), ...
            repmat(string(rule.name), height(score.distance), 1), VariableNames=["truth" "rule"]) score.distance]]; %#ok<AGROW>
        disp(synthetic(end - height(score.distance) + 1:end, :))
    end
end
writetable(synthetic, fullfile(figure_dir, "NR_D_synthetic.csv"))
save(fullfile(out_dir, "NR_D_results.mat"), "summary", "distance", "residuals", "radars", "synthetic", "rules")

%% Plot results
fig = bmm4d.plot.create_export_figure(2400, 1600);
set(fig, "DefaultAxesFontName", "Helvetica", "DefaultTextFontName", "Helvetica")
layout = tiledlayout(fig, 2, 3, "TileSpacing", "tight", "Padding", "tight");
rule_names = string({rules.name});
rule_colors = [0.3 0.3 0.3; lines(numel(rules) - 1)];
for i_case = 1:numel(cases)
    % Column against distance per rule (row i_case: tiles 1-2, synthetic in 3)
    ax = nexttile(layout, (i_case - 1) * 3 + 1);
    hold(ax, "on")
    for i_rule = 1:numel(rules)
        t = distance(distance.case == cases(i_case).name & distance.rule == rule_names(i_rule), :);
        plot(ax, (t.lo_km + t.hi_km) / 2, t.column_ratio, "-o", Color=rule_colors(i_rule, :), MarkerSize=3, ...
            LineWidth=1, DisplayName=rule_names(i_rule))
    end
    yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
    xlabel(ax, "Distance to the nearest radar (km)")
    ylabel(ax, "Column / column at 60-100 km")
    title(ax, cases(i_case).name, FontWeight="bold", Interpreter="none")
    if i_case == 1, legend(ax, Location="southeast", Box="off", FontSize=5), end
    grid(ax, "on")

    % Held-out split entries and kept batch entries at 20-30 km
    ax = nexttile(layout, (i_case - 1) * 3 + 2);
    t = summary(summary.case == cases(i_case).name, :);
    kept = residuals(residuals.case == cases(i_case).name & residuals.cuts == "batch" & ~residuals.dropped & ...
        residuals.lo_km == 20, :);
    [~, order] = ismember(t.rule, kept.rule);
    kept_ratio = nan(height(t), 1);
    kept_ratio(order > 0) = kept.ratio(order(order > 0));
    bar(ax, categorical(t.rule, t.rule), [t.holdout_ratio kept_ratio])
    yline(ax, 1, ":", Color=[0.4 0.4 0.4])
    ylabel(ax, "Observed / predicted")
    legend(ax, ["held-out split < 20 km" "kept batch 20-30 km"], Location="southoutside", Box="off")
    title(ax, "Independent check near the radars", FontWeight="bold")
    grid(ax, "on")
end
% Synthetic: column error against the truth
ax = nexttile(layout, 3);
hold(ax, "on")
markers = ["o" "s" "^"];
synthetic_names = string({synthetic_rules.name});
for i_truth = 1:numel(truths)
    for i_rule = 1:numel(synthetic_rules)
        t = synthetic(synthetic.truth == truths(i_truth).name & synthetic.rule == synthetic_names(i_rule), :);
        plot(ax, (t.lo_km + t.hi_km) / 2, t.truth_ratio, "-" + markers(i_truth), MarkerSize=3, LineWidth=1, ...
            Color=rule_colors(find(rule_names == synthetic_names(i_rule), 1), :), ...
            DisplayName=truths(i_truth).name + ", " + synthetic_names(i_rule))
    end
end
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)")
ylabel(ax, "Reconstructed / true column")
legend(ax, Location="southeast", Box="off", FontSize=4)
title(ax, "Synthetic: does the gap fill?", FontWeight="bold")
grid(ax, "on")
ax = nexttile(layout, 6);
hold(ax, "on")
for i_truth = 1:numel(truths)
    for i_rule = 1:numel(synthetic_rules)
        t = synthetic(synthetic.truth == truths(i_truth).name & synthetic.rule == synthetic_names(i_rule), :);
        plot(ax, (t.lo_km + t.hi_km) / 2, t.truth_error, "-" + markers(i_truth), MarkerSize=3, LineWidth=1, ...
            Color=rule_colors(find(rule_names == synthetic_names(i_rule), 1), :))
    end
end
xlabel(ax, "Distance to the nearest radar (km)")
ylabel(ax, "Median |log(reconstructed / true)|")
title(ax, "Synthetic: column error", FontWeight="bold")
grid(ax, "on")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "NR_D_drop_solves.png"), Resolution=300)
close(fig)
