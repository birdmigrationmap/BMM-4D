%% Near radar E: the drop rule over many radars and nights (dev, long run)
% On the shared night set (studies/nights/PLAN.md): 18 full nights
% (night_us_*: 17 Dark Ecology nights 2021-2025 in six migration phases,
% weak and strong, plus the 21-22 Oct 2026 reference), the 107 us_east
% radars every 30 min, cleaned with production settings without the shape
% correction (cache 493fd773), raw files kept. No download, no cleaning.
% The heavy parts are functions meant for cluster jobs; each skips what is
% already done (files in products/near_radar/NR_E), so this script either
% runs everything in one job or, after the jobs, only collects:
% 1. nr_diagnose_nights: the per-volume diagnostic of every volume
%    (~38,000 volumes, ~0.4 s each; ~1 h on 8 workers);
% 2. consistency: per radar-night deficit and adaptive R_min
%    (nr_radar_night); radar or night (two-way model of the 5-15 km
%    deficit); R_min by season, phase, year and VCP, and against track
%    speed and density (seconds);
% 3. nr_gap_fill: per night and region (Midwest, East), a regional case
%    built from the cleaned volumes and solved with four rules (all rows,
%    cut_set within fixed_range_km, cut_set within the night's adaptive
%    R_min, the whole cut_set), scored as NR_D (~10 min per night-region,
%    36 in all);
% 4. decision table (decision rule: studies/near_radar/PLAN.md).
% Set cut_set and fixed_range_km from the 1 h run (NR_B-NR_D) first.
% Outputs: products/near_radar/NR_E/ (rows_*.mat, gap_*.mat);
% figures/NR_E_*.csv and NR_E_multi_night.png in the study folder, or in
% products/near_radar/NR_E/figures when the study folder is read-only (a
% cluster job's snapshot).

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"), ...
    fullfile(project_root, "studies", "near_radar", "functions"))
products_dir = nr_shared_path("/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/BMM-4D/products");
out_dir = fullfile(products_dir, "near_radar", "NR_E");
figure_dir = fullfile(project_root, "studies", "near_radar", "figures");
probe = fopen(fullfile(figure_dir, ".write_test"), "w");   % a job's snapshot is read-only
if probe < 0
    figure_dir = fullfile(out_dir, "figures");
else
    fclose(probe);
    delete(fullfile(figure_dir, ".write_test"))
end
if ~isfolder(figure_dir), mkdir(figure_dir), end
nights = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string");
nights = nights(startsWith(nights.name, "night_us_"), :);
nights.night = datetime(nights.date);
cut_set = "batch";          % from the 1 h run: "batch" or ["batch" "staggered"]
fixed_range_km = 20;        % from the 1 h run
region_names = ["midwest" "east"];
n_workers = 6;
fprintf("%d nights\n", height(nights))

%% 1. The diagnostic of every volume
row_files = nr_diagnose_nights(nights.name, n_workers);
parts = cell(numel(row_files), 1);
for i_n = 1:numel(row_files)
    loaded = load(row_files(i_n), "rows");
    % The night of every row is its night's date (frames run past midnight UTC)
    loaded.rows.night = repmat(nights.night(i_n), height(loaded.rows), 1);
    parts{i_n} = loaded.rows;
end
rows = vertcat(parts{:});
clear parts loaded

%% 2. Consistency over radars and nights
[radar_night, deficit, ~, ring_edges] = nr_radar_night(rows, cut_set);
radar_night = outerjoin(radar_night, nights(:, ["night" "season" "phase" "speed_ms" "median_refl" "track_deg"]), ...
    Keys="night", Type="left", MergeKeys=true);
radar_night.year = year(radar_night.night);
writetable(radar_night, fullfile(figure_dir, "NR_E_radar_night.csv"))

% Pooled by cut type, and the deficit's spread over radar-nights
by_type = groupsummary(rows, ["waveform" "range_lo_km"], "sum", ["density" "predicted"]);
by_type.ratio = round(by_type.sum_density ./ by_type.sum_predicted, 2);
disp(unstack(by_type(:, ["range_lo_km" "waveform" "ratio"]), "ratio", "waveform"))
writetable(by_type, fullfile(figure_dir, "NR_E_by_type.csv"))
spread = table(ring_edges(1:end-1)', ring_edges(2:end)', sum(isfinite(deficit))', ...
    prctile(deficit, 10)', median(deficit, "omitnan")', prctile(deficit, 90)', ...
    VariableNames=["lo_km" "hi_km" "radar_nights" "p10" "median" "p90"]);
disp(spread(spread.radar_nights > 0, :))
writetable(spread, fullfile(figure_dir, "NR_E_deficit_spread.csv"))
staggered_night = nr_radar_night(rows, "staggered");
staggered_night = staggered_night(staggered_night.has_cuts > 0, :);
writetable(staggered_night, fullfile(figure_dir, "NR_E_radar_night_staggered.csv"))

% Radar or night? Two-way additive model of the log deficit at 5-15 km
near_rings = ring_edges(1:end-1) >= 5 & ring_edges(2:end) <= 15;
radar_night.log_deficit_5_15 = log(mean(deficit(:, near_rings), 2, "omitnan"));
use = isfinite(radar_night.log_deficit_5_15);
model_table = table(categorical(radar_night.station(use)), categorical(string(radar_night.night(use))), ...
    radar_night.log_deficit_5_15(use), VariableNames=["radar" "night" "y"]);
variance = anova(fitlm(model_table, "y ~ radar + night"));
disp(variance)
writetable(variance, fullfile(figure_dir, "NR_E_variance.csv"), WriteRowNames=true)

% R_min against season, phase, year, VCP, track speed and density
fprintf("R_min: median %.1f km, 10-90%% %.1f-%.1f km (%d radar-nights)\n", median(radar_night.r_min_km, "omitnan"), ...
    prctile(radar_night.r_min_km, 10), prctile(radar_night.r_min_km, 90), nnz(isfinite(radar_night.r_min_km)))
for factor = ["season" "phase" "year" "vcp"]
    disp(groupsummary(radar_night(isfinite(radar_night.r_min_km), :), factor, ["median" "mean"], ...
        ["r_min_km" "log_deficit_5_15"]))
end
covariates = radar_night(isfinite(radar_night.r_min_km) & isfinite(radar_night.speed_ms), :);
fprintf("Spearman with R_min: speed %.2f, density %.2f\n", corr(covariates.r_min_km, covariates.speed_ms, Type="Spearman"), ...
    corr(covariates.r_min_km, covariates.median_refl, Type="Spearman"))
fprintf("Radar-nights with R_min > %g km: %.0f%%\n", fixed_range_km, ...
    100 * mean(radar_night.r_min_km(isfinite(radar_night.r_min_km)) > fixed_range_km))

%% 3. Gap filling on regional cases
summary = table();
distance = table();
for region_name = region_names
    gap_files = nr_gap_fill(nights.name', region_name, cut_set, fixed_range_km);
    for gap_file = gap_files
        loaded = load(gap_file, "summary", "distance");
        summary = [summary; loaded.summary]; %#ok<AGROW>
        distance = [distance; loaded.distance]; %#ok<AGROW>
    end
end
writetable(summary, fullfile(figure_dir, "NR_E_summary.csv"))
writetable(distance, fullfile(figure_dir, "NR_E_distance.csv"))

%% 4. Decision table
decision = groupsummary(summary, "rule", {"median", @(x) prctile(x, 10), @(x) prctile(x, 90)}, ...
    ["column_0_10" "column_10_20" "holdout_ratio" "kept_batch_20_40" "dropped_share" "misfit_trusted"]);
disp(decision)
writetable(decision, fullfile(figure_dir, "NR_E_decision.csv"))

%% Plot results
fig = bmm4d.plot.create_export_figure(2400, 900);
set(fig, "DefaultAxesFontName", "Helvetica", "DefaultTextFontName", "Helvetica")
layout = tiledlayout(fig, 1, 3, "TileSpacing", "tight", "Padding", "tight");

% Deficit by ring: median and 10-90% over radar-nights
ax = nexttile(layout);
hold(ax, "on")
shown = spread(spread.radar_nights >= 20, :);   % rings most radar-nights reach
x = shown.lo_km;
fill(ax, [x; flipud(x)], [shown.p10; flipud(shown.p90)], [0.85 0.33 0.10], FaceAlpha=0.2, EdgeColor="none", ...
    DisplayName="10-90% of radar-nights")
plot(ax, x, shown.median, "-o", Color=[0.85 0.33 0.10], MarkerSize=4, LineWidth=1.2, DisplayName="median")
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Range, inner edge of the ring (km)")
ylabel(ax, "Deficit (ring / own level at 25-50 km)")
legend(ax, Location="southeast", Box="off")
title(ax, sprintf("%d radar-nights", height(radar_night)), FontWeight="bold")
grid(ax, "on")

% Adaptive R_min by night
ax = nexttile(layout);
night_labels = categorical(string(radar_night.night, "yyyy-MM-dd"));
boxchart(ax, night_labels(isfinite(radar_night.r_min_km)), radar_night.r_min_km(isfinite(radar_night.r_min_km)), ...
    BoxFaceColor=[0.85 0.33 0.10])
yline(ax, fixed_range_km, ":", "fixed rule", Color=[0.4 0.4 0.4])
ylabel(ax, "Adaptive R_{min} (km)")
title(ax, "Per radar, by night", FontWeight="bold")
grid(ax, "on")

% Column near the radars by rule, over night-regions
ax = nexttile(layout);
hold(ax, "on")
rule_names = unique(summary.rule, "stable");
rule_colors = [0.3 0.3 0.3; lines(numel(rule_names) - 1)];
for i_rule = 1:numel(rule_names)
    t = summary(summary.rule == rule_names(i_rule), :);
    swarmchart(ax, repmat(i_rule, height(t), 1), t.column_0_10, 12, rule_colors(i_rule, :), "filled", ...
        MarkerFaceAlpha=0.6)
end
yline(ax, 1, ":", Color=[0.4 0.4 0.4])
xticks(ax, 1:numel(rule_names))
xticklabels(ax, rule_names)
ylabel(ax, "Column 0-10 km / 60-100 km")
title(ax, "Gap filling, per night and region", FontWeight="bold")
grid(ax, "on")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "NR_E_multi_night.png"), Resolution=300)
close(fig)
