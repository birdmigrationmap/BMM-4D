%% Near radar B: steep cuts against the split-cut profile, every cleaned volume (dev)
% The cheap diagnostic, no solve: in each cleaned volume, every sweep and
% range ring against the beam-weighted vertical profile fitted to the
% split cuts at 20-60 km of the same volume (nr_volume_ratio). The cut
% type comes from the Level II header (nr_read_vcp). Answers, over many
% radars and nights:
% - which cuts read low near the radar (batch, Doppler, staggered; by VCP);
% - how the deficit depends on range and height, and how far out it goes;
% - whether it is the same at every radar and night (fixed rule) or not
%   (adaptive rule): per radar-night, R_min = the outer edge of the
%   outermost range ring where the batch cuts read below the split cuts
%   of that radar-night by more than two standard errors over volumes.
% Sets (estimator cleanings on disk):
%   us_east_20261021 - cache 9c2b4ef2 (us_east production cleaning, no
%                      aspect correction), ~101 radars, 24 frames;
%   night_20210512   - cache fad5e1e1, 28 Midwest radars, 13 May 2021
%                      02-10 UTC (older estimator cleaning with the
%                      3000 birds km^-3 cap and patch filter);
%   sandbox_nights   - cache fad5e1e1, 9 nights x 28 radars, one volume each.
% every_kth thins the volume lists of the two full nights (1 = all, ~2800
% volumes, ~25 min serial; light first look: 6, ~7 min).
% Outputs: <local_data_dir>/near_radar/NR_B_rows.mat (all rows),
% figures/NR_B_pooled.csv, NR_B_radar_night.csv, NR_B_volume_ratio.png.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "near_radar", "functions"))
figure_dir = fullfile(project_root, "studies", "near_radar", "figures");
out_dir = fullfile(bmm4d.util.local_data_dir(), "near_radar");
if ~isfolder(out_dir), mkdir(out_dir), end
volume_root = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes");
sandbox_raw_dir = fullfile(bmm4d.util.local_data_dir(), "sandbox_nights", "raw");
every_kth = 1;

%% List the cleaned volumes
sets = struct("name", {"us_east_20261021", "night_20210512", "sandbox_nights"}, ...
    "folder", {"9c2b4ef2", "fad5e1e1", "fad5e1e1"}, "pattern", {"K*_V06.mat", "K*20210513_*_V06.mat", "K*_V06.mat"});
volumes = table();
for volume_set = sets
    listed = dir(fullfile(volume_root, volume_set.folder, volume_set.pattern));
    names = string({listed.name})';
    if volume_set.name == "sandbox_nights"
        % One volume per radar and night: those read from sandbox_nights/raw
        names = names(isfile(fullfile(sandbox_raw_dir, erase(names, ".mat"))));
    else
        names = names(1:every_kth:end);
    end
    volumes = [volumes; table(fullfile(volume_root, volume_set.folder, names), repmat(string(volume_set.name), numel(names), 1), ...
        VariableNames=["file" "set"])]; %#ok<AGROW>
end
fprintf("%d volumes\n", height(volumes))

%% Diagnostic of every volume
parts = cell(height(volumes), 1);
run_tic = tic;
for i_v = 1:height(volumes)
    volume = load(volumes.file(i_v), "volume").volume;
    raw_file = volume.source_file;
    if ~isfile(raw_file)
        [~, raw_name] = fileparts(volumes.file(i_v));
        raw_file = fullfile(sandbox_raw_dir, raw_name);
    end
    if ~isfile(raw_file), continue, end
    vcp = nr_read_vcp(raw_file);
    if isempty(vcp.cuts), continue, end
    rows = nr_volume_ratio(volume, vcp.cuts);
    if isempty(rows), continue, end
    rows.set = repmat(volumes.set(i_v), height(rows), 1);
    rows.vcp = repmat(vcp.number, height(rows), 1);
    parts{i_v} = rows;
    if mod(i_v, 100) == 0
        fprintf("%d / %d volumes, %.0f s\n", i_v, height(volumes), toc(run_tic))
    end
end
rows = vertcat(parts{:});
% Night: the date of the evening (UTC minus 16 h: US night frames run
% from ~22 to ~13 UTC)
rows.night = dateshift(rows.datetime - hours(16), "start", "day");
save(fullfile(out_dir, "NR_B_rows.mat"), "rows")
fprintf("%d rows from %d volumes in %.0f s\n", height(rows), numel(unique(rows.station + string(rows.datetime))), toc(run_tic))

%% Pooled: cut type, range ring and height band
height_edges = [0 300 600 900 1200 1500 2000 3000];
rows.height_band = discretize(rows.height_m, height_edges);
pooled = groupsummary(rows, ["waveform" "range_lo_km" "height_band"], "sum", ["n" "density" "predicted"]);
pooled.ratio = pooled.sum_density ./ pooled.sum_predicted;
pooled = pooled(pooled.sum_n >= 200, :);
writetable(pooled, fullfile(figure_dir, "NR_B_pooled.csv"))
by_type = groupsummary(rows, ["waveform" "range_lo_km"], "sum", ["density" "predicted"]);
by_type.ratio = round(by_type.sum_density ./ by_type.sum_predicted, 2);
disp(unstack(by_type(:, ["range_lo_km" "waveform" "ratio"]), "ratio", "waveform"))
batch = rows(rows.waveform == "batch", :);
by_height = groupsummary(batch, ["range_lo_km" "height_band"], "sum", ["density" "predicted"]);
by_height.ratio = round(by_height.sum_density ./ by_height.sum_predicted, 2);
disp("Batch cuts by range ring (rows) and height band (columns: 0-300, 300-600, ... m)")
disp(unstack(by_height(:, ["range_lo_km" "height_band" "ratio"]), "ratio", "height_band", VariableNamingRule="preserve"))

%% Per radar-night: near-radar deficit of the batch cuts and adaptive R_min
% nr_radar_night: each ring of the batch cuts over the radar-night's own
% batch level at 25-50 km (deficit, jackknife se over volumes), the far
% batch level over the split cuts (far_level, a cut-type offset), and
% R_min = outer edge of the outermost ring within 25 km that is low by more
% than two standard errors. The staggered cuts get the same table.
[radar_night, deficit, ~, ring_edges] = nr_radar_night(rows, "batch");
writetable(radar_night, fullfile(figure_dir, "NR_B_radar_night.csv"))
far_level = radar_night.far_level;
r_min = radar_night.r_min_km;
staggered_night = nr_radar_night(rows, "staggered");
staggered_night = staggered_night(staggered_night.has_cuts > 0, :);
writetable(staggered_night, fullfile(figure_dir, "NR_B_radar_night_staggered.csv"))
fprintf("Staggered cuts: %d radar-nights, R_min median %.1f km\n", height(staggered_night), ...
    median(staggered_night.r_min_km, "omitnan"))

% Spread over radar-nights of the deficit by ring
spread = table(ring_edges(1:end-1)', ring_edges(2:end)', sum(isfinite(deficit))', ...
    prctile(deficit, 10)', median(deficit, "omitnan")', prctile(deficit, 90)', ...
    VariableNames=["lo_km" "hi_km" "radar_nights" "p10" "median" "p90"]);
disp("Batch deficit (ring / the radar-night's batch level at 25-50 km), spread over radar-nights")
disp(spread(spread.radar_nights > 0, :))
fprintf("Far batch level (batch 25-50 km / split 10-40 km): median %.2f, 10-90%% %.2f-%.2f\n", ...
    median(far_level, "omitnan"), prctile(far_level, 10), prctile(far_level, 90))
fprintf("Adaptive R_min (radar-nights with >= 3 volumes): median %.1f km, 10-90%% %.1f-%.1f km, n = %d\n", ...
    median(r_min, "omitnan"), prctile(r_min, 10), prctile(r_min, 90), nnz(isfinite(r_min)))
disp(groupsummary(radar_night(isfinite(radar_night.r_min_km), :), ["set" "vcp"], ["median" "min" "max"], "r_min_km"))

%% Plot results
fig = bmm4d.plot.create_export_figure(2400, 900);
set(fig, "DefaultAxesFontName", "Helvetica", "DefaultTextFontName", "Helvetica")
layout = tiledlayout(fig, 1, 3, "TileSpacing", "tight", "Padding", "tight");

% Pooled ratio against range by cut type
ax = nexttile(layout);
hold(ax, "on")
types = ["split" "batch" "staggered" "doppler"];
colors = [0.3 0.3 0.3; 0.85 0.33 0.10; 0.1 0.45 0.6; 0.45 0.65 0.2];
for i_type = 1:numel(types)
    t = by_type(by_type.waveform == types(i_type), :);
    if isempty(t), continue, end
    plot(ax, t.range_lo_km, t.sum_density ./ t.sum_predicted, "-o", Color=colors(i_type, :), MarkerSize=4, ...
        LineWidth=1.2, DisplayName=types(i_type))
end
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Range, inner edge of the ring (km)")
ylabel(ax, "Observed / split-cut profile")
legend(ax, Location="southeast", Box="off")
title(ax, "By cut type, all volumes", FontWeight="bold")
grid(ax, "on")

% Batch cuts by height band
ax = nexttile(layout);
hold(ax, "on")
bands = unique(by_height.height_band)';
band_colors = copper(numel(bands) + 1);
for i_b = 1:numel(bands)
    t = by_height(by_height.height_band == bands(i_b), :);
    plot(ax, t.range_lo_km, t.sum_density ./ t.sum_predicted, "-o", Color=band_colors(i_b, :), MarkerSize=4, ...
        LineWidth=1.2, DisplayName=sprintf("%d-%d m", height_edges(bands(i_b)), height_edges(bands(i_b) + 1)))
end
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Range, inner edge of the ring (km)")
ylabel(ax, "Observed / split-cut profile")
legend(ax, Location="southeast", Box="off")
title(ax, "Batch cuts by height above the antenna", FontWeight="bold")
grid(ax, "on")

% Spread over radar-nights
ax = nexttile(layout);
hold(ax, "on")
shown = spread(spread.radar_nights >= 20, :);   % rings most radar-nights reach
x = shown.lo_km;
fill(ax, [x; flipud(x)], [shown.p10; flipud(shown.p90)], [0.85 0.33 0.10], FaceAlpha=0.2, EdgeColor="none", ...
    DisplayName="10-90% of radar-nights")
plot(ax, x, shown.median, "-o", Color=[0.85 0.33 0.10], MarkerSize=4, LineWidth=1.2, DisplayName="median")
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Range, inner edge of the ring (km)")
ylabel(ax, "Batch / own batch level at 25-50 km")
legend(ax, Location="southeast", Box="off")
title(ax, sprintf("Per radar-night (%d)", height(radar_night)), FontWeight="bold")
grid(ax, "on")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "NR_B_volume_ratio.png"), Resolution=300)
close(fig)
