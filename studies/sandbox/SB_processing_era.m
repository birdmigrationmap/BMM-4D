%% Sandbox: near-range deficit before and after the radar's point clutter filter
% Within one volume the radar's own reflectivity reads lower the closer a
% gate is to the radar at the same height (SB_near_raw, SB_sweep_offset;
% the gate-density report). Leading explanation (H1): near the
% radar a gate holds less than one bird, so each bird is an isolated
% strong spike, which the radar's point clutter filter (RDA Build 19,
% 2019-2020, with the CFP moment) replaces by the background. Test: the
% same measurement on volumes from 2015-2018, before that processing.
% Per volume, sweep (one pass per elevation: the one with the most REF
% gates, i.e. the surveillance pass of split cuts; 0.5 deg left out for its
% ground clutter) and 100 m of height (300-2500 m above the antenna):
% mean linear Z (below threshold = 0) by 1 km range bin, 3-60 km, and the
% share of gates above 20 dBZ (a single bird in a small gate). Range
% factor fitted on pairs of adjacent sweeps at the same height, as in
% SB_near_raw (1 at 30-60 km).
% Volumes: the 100 rain-free gate-study scans of 2015-2021 (GD_E) and the
% 124 of 2021-2025 (GD_J), grouped by year.
% Outputs: figures/processing_era.png, processing_era.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
range_edges_m = 3000:1000:60000;
height_edges = 300:100:2500;
fit_edges_km = [3 5 7.5 10 15 20 30 45 60];

%% Volumes ------------------------------------------------------------------
selected = dir("/Users/rafnuss/Data/BMM-within/gate_density/selected_scans_volumes/K*.mat");
raw_files = strings(numel(selected), 1);
for i_f = 1:numel(selected)
    raw_files(i_f) = load(fullfile(selected(i_f).folder, selected(i_f).name)).volume.source_file;
end
post2020 = readtable(fullfile(project_root, "studies", "gate_density", "figures", "GD_J_selected_scans.csv"), TextType="string");
raw_files = [raw_files; post2020.raw_file];
[~, names] = fileparts(raw_files);
volumes = table(raw_files, extractBefore(names, 5), double(extractBetween(names, 5, 8)), ...
    VariableNames=["file" "radar" "year"]);
volumes = volumes(isfile(volumes.file), :);
volumes.era = discretize(volumes.year, [2015 2019 2020 2026], "categorical", ["2015-2018" "2019" "2020-2025"]);
disp(groupcounts(volumes, "era"))

%% Statistics per volume ------------------------------------------------------
stats = cell(height(volumes), 1);
parfor i_v = 1:height(volumes)
    try
        stats{i_v} = volume_stats(volumes.file(i_v), range_edges_m, height_edges);
    catch
        stats{i_v} = [];
    end
end

%% Pairs of adjacent sweeps at the same height ---------------------------------
pair_rows = cell(height(volumes), 1);
for i_v = 1:height(volumes)
    s = stats{i_v};
    if numel(s) < 2, continue, end
    r_lower = [s(1:end - 1).range_m]; r_upper = [s(2:end).range_m];
    z_lower = [s(1:end - 1).mean_z]; z_upper = [s(2:end).mean_z];
    ok = ~isnan(r_lower) & ~isnan(r_upper) & z_lower > 0 & z_upper > 0;
    pair_rows{i_v} = table(repmat(volumes.era(i_v), nnz(ok), 1), r_lower(ok), r_upper(ok), log2(z_upper(ok) ./ z_lower(ok)), ...
        VariableNames=["era" "r_lower" "r_upper" "log_ratio"]);
end
pairs = vertcat(pair_rows{:});
fprintf("%d pairs\n", height(pairs))

%% Range factor and spikes by era -------------------------------------------------
eras = categories(volumes.era);
range_factor = nan(numel(fit_edges_km) - 1, numel(eras));
range_centre_km = (range_edges_m(1:end - 1) + range_edges_m(2:end)) / 2000;
spike_share = nan(numel(range_centre_km), numel(eras));
for i_e = 1:numel(eras)
    in_era = pairs.era == eras{i_e};
    if nnz(in_era) < 50, continue, end
    range_factor(:, i_e) = fit_range(pairs(in_era, :), fit_edges_km);
    % Share of gates above 20 dBZ by range, 300-1500 m, all sweeps
    n_spike = zeros(numel(range_centre_km), 1); n_gate = n_spike;
    for i_v = find(volumes.era == eras{i_e})'
        for sweep = stats{i_v}
            n_spike = n_spike + sweep.n_spike;
            n_gate = n_gate + sweep.n_gate;
        end
    end
    spike_share(:, i_e) = n_spike ./ n_gate;
end
fit_centre = (fit_edges_km(1:end - 1) + fit_edges_km(2:end)) / 2;
result = array2table(round(range_factor, 2), VariableNames=eras);
result.range_km = fit_centre';
disp(result)
writetable(result, fullfile(figure_dir, "processing_era.csv"))

%% Beam smoothing alone: synthetic profiles through the same pair fit ------------
% A far sweep's beam is wider (0.95 deg two-way Gaussian, sd ~ r), so at
% the same nominal height it averages a thicker layer. Five profile shapes,
% sampled every 250 m along each sweep (0.9-6.4 deg), 100 m height bins.
elevations = [0.9 1.3 1.8 2.4 3.1 4.0 5.1 6.4];
range_m = (3000:250:60000)';
z_fine = (0:10:5000)';
sd_two_way = deg2rad(0.95) / sqrt(2) / (2 * sqrt(2 * log(2)));
profiles = {@(z) exp(-z / 600), @(z) exp(-(z - 1000) .^ 2 / (2 * 300 ^ 2)), @(z) exp(-(z - 1600) .^ 2 / (2 * 400 ^ 2)), ...
    @(z) double(z < 1500), @(z) double(z >= 500 & z < 2000)};
profile_names = ["exponential, 600 m" "layer at 1000 m" "layer at 1600 m" "uniform to 1500 m" "empty below 500 m"];
synthetic_factor = nan(numel(fit_edges_km) - 1, numel(profiles));
for i_p = 1:numel(profiles)
    beam_mean = nan(numel(height_edges) - 1, numel(elevations));
    beam_range = beam_mean;
    for i_e = 1:numel(elevations)
        [~, centre_m] = slant2ground(range_m, elevations(i_e));
        weight = exp(-(z_fine' - centre_m) .^ 2 ./ (2 * (range_m * sd_two_way) .^ 2));
        observed = (weight * profiles{i_p}(z_fine)) ./ sum(weight, 2);
        height_bin = discretize(centre_m, height_edges);
        ok = ~isnan(height_bin);
        beam_mean(:, i_e) = accumarray(height_bin(ok), observed(ok), [numel(height_edges) - 1 1], @mean, NaN);
        beam_range(:, i_e) = accumarray(height_bin(ok), range_m(ok), [numel(height_edges) - 1 1], @mean, NaN);
    end
    lower_mean = beam_mean(:, 1:end - 1); upper_mean = beam_mean(:, 2:end);
    ok = lower_mean > 1e-3 & upper_mean > 1e-3;
    lower_range = beam_range(:, 1:end - 1); upper_range = beam_range(:, 2:end);
    synthetic_pairs = table(lower_range(ok), upper_range(ok), log2(upper_mean(ok) ./ lower_mean(ok)), ...
        VariableNames=["r_lower" "r_upper" "log_ratio"]);
    synthetic_factor(:, i_p) = fit_range(synthetic_pairs, fit_edges_km);
end
disp(array2table(round(synthetic_factor, 2), VariableNames=profile_names))

%% Plot results -----------------------------------------------------------
era_colors = [0.47 0.67 0.19; 0.6 0.6 0.6; 0.85 0.33 0.10];
fig = bmm4d.plot.create_export_figure(3000, 1150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Radar output near the radar, before and after the point clutter filter (RDA Build 19, 2019-2020)", ...
    FontSize=16, Color="k")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
fill(ax, [fit_centre fliplr(fit_centre)], [min(synthetic_factor, [], 2)' fliplr(max(synthetic_factor, [], 2)')], ...
    [0.85 0.85 0.85], EdgeColor="none", DisplayName="beam smoothing alone (5 profile shapes)")
for i_e = 1:numel(eras)
    n_volume = nnz(volumes.era == eras{i_e} & ~cellfun(@isempty, stats));
    plot(ax, fit_centre, range_factor(:, i_e), "-o", Color=era_colors(i_e, :), LineWidth=2, ...
        MarkerFaceColor=era_colors(i_e, :), DisplayName=sprintf("%s (%d volumes)", eras{i_e}, n_volume))
end
yline(ax, 1, "k:", HandleVisibility="off")
xlim(ax, [0 60]); ylim(ax, [0 1.3]); set(ax, FontSize=13)
xlabel(ax, "Range (km)"); ylabel(ax, "Radar output at the same height (1 at 30-60 km)")
legend(ax, Location="southeast", FontSize=12)
title(ax, "Range factor from sweep pairs at the same height", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_e = 1:numel(eras)
    plot(ax, range_centre_km, 100 * spike_share(:, i_e), "-", Color=era_colors(i_e, :), LineWidth=2, DisplayName=eras{i_e})
end
xlim(ax, [0 60]); set(ax, FontSize=13, YScale="log")
xlabel(ax, "Range (km)"); ylabel(ax, "Gates above 20 dBZ (%), 300-1500 m")
legend(ax, Location="northeast", FontSize=12)
title(ax, "Bright single gates (one bird in a small gate)", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "processing_era.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function s = volume_stats(raw_file, range_edges_m, height_edges)
% Per sweep (one pass per elevation, 0.5 deg left out), by height bin: mean
% linear Z and mean range of the gates; by range bin: gates and gates above
% 20 dBZ at 300-1500 m above the antenna.
sweeps = bmm4d.radar.read_level2_moments(raw_file, "REF", 7);
sweeps = sweeps(arrayfun(@(sw) ~isempty(sw.REF), sweeps));
elev = round([sweeps.elev] * 10) / 10;
n_gates = arrayfun(@(sw) numel(sw.REF_range), sweeps);
s = struct("elev", {}, "mean_z", {}, "range_m", {}, "n_spike", {}, "n_gate", {});
for e = unique(elev(elev >= 0.8))
    candidates = find(elev == e);
    [~, k] = max(n_gates(candidates));
    sweep = sweeps(candidates(k));
    range_m = sweep.REF_range;
    z = 10 .^ (sweep.REF / 10);
    z(sweep.REF_code == 0) = 0;
    valid = sweep.REF_code ~= 1;
    [~, height_m] = slant2ground(range_m, sweep.elev);
    in_range = range_m >= range_edges_m(1) & range_m < range_edges_m(end);
    height_bin = discretize(height_m, height_edges);
    mean_z = nan(numel(height_edges) - 1, 1);
    mean_range = nan(numel(height_edges) - 1, 1);
    for i_h = 1:numel(height_edges) - 1
        rows = in_range & height_bin == i_h;
        g = valid & rows;
        if nnz(g) >= 300
            mean_z(i_h) = mean(z(g));
            mean_range(i_h) = mean(range_m(rows));
        end
    end
    range_bin = discretize(range_m, range_edges_m);
    low = height_m >= 300 & height_m <= 1500 & ~isnan(range_bin);
    n_range = numel(range_edges_m) - 1;
    n_gate = accumarray(range_bin(low), sum(valid(low, :), 2), [n_range 1]);
    n_spike = accumarray(range_bin(low), sum(valid(low, :) & sweep.REF(low, :) > 20, 2), [n_range 1]);
    s(end + 1) = struct("elev", e, "mean_z", mean_z, "range_m", mean_range, "n_spike", n_spike, "n_gate", n_gate); %#ok<AGROW>
end
end

function range_factor = fit_range(pairs, range_edges_km)
% Robust fit of the pairs' log2 ratio as f(r_upper) - f(r_lower), f with
% mean 0 at 30-60 km; returns 2^f per range bin.
n = height(pairs);
n_bin = numel(range_edges_km) - 1;
upper_bin = discretize(pairs.r_upper / 1000, range_edges_km);
lower_bin = discretize(pairs.r_lower / 1000, range_edges_km);
A = sparse([1:n 1:n]', [upper_bin; lower_bin], [ones(n, 1); -ones(n, 1)], n, n_bin);
anchor = zeros(1, n_bin);
anchor(range_edges_km(1:end - 1) >= 30) = 1 / nnz(range_edges_km(1:end - 1) >= 30);
y = pairs.log_ratio;
weight = ones(n, 1);
for i_pass = 1:6
    x = [A .* sqrt(weight); 100 * anchor] \ [y .* sqrt(weight); 0];
    residual = y - A * x;
    scale = 1.4826 * median(abs(residual));
    weight = min(1, 1.345 * scale ./ abs(residual));
end
range_factor = 2 .^ x;
end
