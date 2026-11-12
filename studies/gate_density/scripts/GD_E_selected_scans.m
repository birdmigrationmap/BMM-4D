%% Gate density E: the near-range deficit on 100 selected rain-free scans
% The 100 scans selected for the BMM-within simulator (03_simulator_fit,
% S03_A_select_scans: 31 radars, 2015-2021, spring and autumn nights with
% migration, checked by eye for rain). All elevation cuts up to 10 deg,
% gates to 80 km. For every gate, densities by:
% - output mean: the radar output, gate by gate;
% - estimator: the homogeneous-density estimator (beam from the radar
%   specifications, assumed cross-section log-sd 1.25, BMM4D.RADAR.GATE.BEAM), +-1 km window;
% - estimator with log-sd 1 and 1.75: the sensitivity to the assumed spread;
% - old filter cleaning: bmm4d.radar.clean_volume with the sandbox settings
%   (BMM-4D 0.1.0 cleaning; cut at 7 deg), cached.
% Densities are averaged per scan, elevation cut, range bin and 100 m height
% layer (beam centre above ground). The deficit D compares each bin with the
% same scan's density at 30-60 km in the same 100 m layer, pooled over cuts;
% thin layers keep the vertical profile from leaking into D.
% Two views:
% 1. D(r) per layer: median and interquartile range over scans.
% 2. D(r) per elevation cut, over all layers: the sampling deficit depends on
%    range only (gate volume), while the radar's clutter filter acts mostly
%    at the lowest cuts, so differences between cuts at the same range
%    point to the clutter filter.
% 3. Sensitivity: estimate with log-sd 1 and 1.75 over the estimate with
%    1.25, per range bin and elevation class, median over scans.
% Outputs: studies/gate_density/figures/GD_E_*.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
% The scan selection lives in the BMM-within project (data/derived/scans.mat)
scans_project_dir = "/Users/rafnuss/Library/CloudStorage/OneDrive-Vogelwarte/BMM-within";
% Assumed cross-section log-sd 1.25 (BMM4D.RADAR.GATE.BEAM), and 1 and 1.75 for the sensitivity
spread_values = [1.25 1 1.75];
n_spread = numel(spread_values);
beams = arrayfun(@(v) bmm4d.radar.gate.beam(sigma_log_sd=v), spread_values, UniformOutput=false);
tables = cellfun(@(b) bmm4d.radar.gate.bird_sum_table(b, 10000, 20000), beams, UniformOutput=false);   % covers 2000 birds km^-3 at 80 km
s_floor_cm2 = 0.1;
radius_m = 1000;
max_range_m = 80000;
max_elevation_deg = 10;
max_density = 3000;
bin_edges_km = [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80];
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end)) / 2;
n_bin = numel(bin_km);
layer_edges_m = 100:100:3000;
n_layer = numel(layer_edges_m) - 1;
% Thicker layers for the figures and tables
display_edges_m = [300 500; 500 800; 800 1200; 1200 1700];
display_label = compose("%d-%d m", display_edges_m(:, 1), display_edges_m(:, 2));
n_display = numel(display_label);
% Elevation cuts, grouped in five classes (VCP 212/215 cuts; older VCP 21
% cuts 1.5, 3.4, 4.3, 6.0 and 9.9 deg fall in the same classes)
elevation_edges = [0 1.1 2.1 3.6 5.75 max_elevation_deg + 0.5];
elevation_label = ["0.5-0.9°" "1.3-1.8°" "2.4-3.4°" "4.0-5.1°" "6.0-10°"];
n_elevation = numel(elevation_label);
method_name = ["output mean" "estimator" "old filter cleaning" "estimator, s_\sigma 1" "estimator, s_\sigma 1.75"];
estimate_method = [2 4 5];                       % method of each spread value
n_method = numel(method_name);
earth_radius_m = 4 / 3 * 6371000;

%% Selected scans --------------------------------------------------------------
scans = load(fullfile(scans_project_dir, "data", "derived", "scans.mat")).scans;
radar_table = readtable(fullfile(scans_project_dir, "data", "reference", "radar.csv"), TextType="string");
[~, base, ext] = fileparts(scans.filename);
scans.raw_file = fullfile(scans_project_dir, "data", "raw_scans", "nexrad", ...
    extractBefore(base, 5) + "_" + extractBetween(base, 5, 12), base + ext);
[~, i_radar] = ismember(scans.station, radar_table.name);
scans.antenna_agl = radar_table.height(i_radar) - radar_table.dem(i_radar);
n_scan = height(scans);
assert(all(isfile(scans.raw_file)), "Missing raw scans in " + scans_project_dir)

%% Old filter cleaning of the selected scans (cached) ---------------------------
% Same cleaning settings as the sandbox case (read_options, vp, filter)
c = case_settings("sandbox");
volume_dir = fullfile(bmm4d.util.local_data_dir(), "gate_density", "selected_scans_volumes");
if ~isfolder(volume_dir)
    mkdir(volume_dir)
end
scans.volume_file = fullfile(volume_dir, base + ".mat");
todo = find(~isfile(scans.volume_file));
read_options = c.read_options;
vp = c.vp;
filter = c.filter;
status = strings(numel(todo), 1);
todo_raw = scans.raw_file(todo);
todo_out = scans.volume_file(todo);
% MistNet needs about 2.6 GB per worker (bmm4d.clean_volumes)
if ~isempty(todo) && isempty(gcp("nocreate"))
    parpool("Processes", 4);
end
parfor k = 1:numel(todo)
    status(k) = bmm4d.radar.clean_volume(todo_raw(k), read_options, vp, filter, todo_out(k));
end
fprintf("[selected scans] old cleaning: %d done, %d failed\n", nnz(status == "ok"), nnz(status ~= "ok"))

%% Densities per range bin, layer, elevation cut and scan ------------------------
density_sum = zeros(n_bin, n_layer, n_elevation, n_scan, n_method);
gate_count = zeros(n_bin, n_layer, n_elevation, n_scan, n_method);
for i_s = 1:n_scan
    level2 = bmm4d.radar.read_level2_moments(scans.raw_file(i_s), "REF", max_elevation_deg);
    level2 = level2(arrayfun(@(s) ~isempty(s.REF), level2));
    % Split cuts repeat an elevation: keep the first (surveillance) cut
    [~, keep] = unique(round([level2.elev] * 2.5) / 2.5, "stable");
    level2 = level2(keep);
    cleaned = [];
    if isfile(scans.volume_file(i_s))
        cleaned = load(scans.volume_file(i_s)).volume.dens.sweeps;
    end
    for i_e = 1:numel(level2)
        i_cut = discretize(level2(i_e).elev, elevation_edges);
        n_gate = min(size(level2(i_e).REF, 1), find(level2(i_e).REF_range <= max_range_m, 1, "last"));
        range_m = level2(i_e).REF_range(1:n_gate);
        dbz = level2(i_e).REF(1:n_gate, :);
        az_step_deg = 360 / size(dbz, 2);
        % Clutter and rain (> 3000 birds km^-3): left out of every method
        excluded = z_to_refl(idb(dbz)) / 11 > max_density;
        output = z_to_refl(idb(dbz)) / 11;
        output(isnan(dbz)) = 0;
        output(excluded) = nan;
        detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        per_gate = cell(1, n_method);
        per_gate{1} = output;
        for i_v = 1:n_spread
            [lo, hi] = bmm4d.radar.gate.bounds(dbz, isnan(dbz), range_m, beams{i_v}, detection_dbz, s_floor_cm2);
            lo(excluded) = -Inf; hi(excluded) = Inf;
            estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, az_step_deg, beams{i_v}, tables{i_v}, radius_m=radius_m);
            per_gate{estimate_method(i_v)} = estimate.density;
        end
        % Average per range bin, then assign the bin to its layer
        beam_agl = scans.antenna_agl(i_s) + bin_km' * 1000 * sind(level2(i_e).elev) + (bin_km' * 1000) .^ 2 / (2 * earth_radius_m);
        layer = discretize(beam_agl, layer_edges_m);
        bin_of_gate = discretize(range_m / 1000, bin_edges_km);
        has_clean = ~isempty(cleaned) && min(abs([cleaned.elev] - level2(i_e).elev)) < 0.2;
        if has_clean
            [~, i_clean] = min(abs([cleaned.elev] - level2(i_e).elev));
            bin_of_clean = discretize(cleaned(i_clean).range(:) / 1000, bin_edges_km);
        end
        for k = find(~isnan(layer))'
            rows = bin_of_gate == k;
            if ~any(rows)
                continue
            end
            for i_m = [1 estimate_method]
                values = per_gate{i_m}(rows, :);
                density_sum(k, layer(k), i_cut, i_s, i_m) = density_sum(k, layer(k), i_cut, i_s, i_m) + sum(values, "all", "omitnan");
                gate_count(k, layer(k), i_cut, i_s, i_m) = gate_count(k, layer(k), i_cut, i_s, i_m) + nnz(~isnan(values));
            end
            if has_clean
                values = cleaned(i_clean).data(bin_of_clean == k, :);
                density_sum(k, layer(k), i_cut, i_s, 3) = density_sum(k, layer(k), i_cut, i_s, 3) + sum(values, "all", "omitnan");
                gate_count(k, layer(k), i_cut, i_s, 3) = gate_count(k, layer(k), i_cut, i_s, 3) + nnz(~isnan(values));
            end
        end
    end
    fprintf("[selected scans] %d/%d %s done\n", i_s, n_scan, scans.station(i_s))
end
save(fullfile(figure_dir, "GD_E_selected_scans.mat"), "density_sum", "gate_count", "bin_km", "layer_edges_m", ...
    "elevation_label", "method_name", "scans")

%% Deficit relative to 30-60 km in the same layer -------------------------------
% Reference: all cuts pooled, per scan, 100 m layer and method
reference = bin_km >= 30 & bin_km <= 60;
reference_density = sum(density_sum(reference, :, :, :, :), [1 3]) ./ sum(gate_count(reference, :, :, :, :), [1 3]);
reference_density(reference_density < 5) = nan;                  % too few birds to compare
relative = density_sum ./ gate_count ./ reference_density;        % [bin layer cut scan method]
weight = gate_count .* ~isnan(relative);
relative(isnan(relative)) = 0;
% Per display layer, all cuts pooled: gate-weighted mean of D over its 100 m layers
layer_mid_m = (layer_edges_m(1:end - 1) + layer_edges_m(2:end)) / 2;
layer_relative = nan(n_bin, n_display, n_scan, n_method);
for i_d = 1:n_display
    in_display = layer_mid_m > display_edges_m(i_d, 1) & layer_mid_m < display_edges_m(i_d, 2);
    layer_relative(:, i_d, :, :) = sum(relative(:, in_display, :, :, :) .* weight(:, in_display, :, :, :), [2 3]) ./ ...
        sum(weight(:, in_display, :, :, :), [2 3]);
end
% Per elevation class, all layers pooled
cut_relative = squeeze(sum(relative .* weight, 2) ./ sum(weight, 2));   % [bin cut scan method]

summary = table(bin_km', VariableNames="range_km");
for i_d = [1 3]
    for i_m = 1:n_method
        summary.(matlab.lang.makeValidName(display_label(i_d) + "_" + method_name(i_m))) = ...
            round(median(layer_relative(:, i_d, :, i_m), 3, "omitnan"), 2);
    end
end
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_E_deficit.csv"))
cut_summary = table(bin_km', VariableNames="range_km");
for i_m = 1:n_method
    for i_c = 1:n_elevation
        cut_summary.(matlab.lang.makeValidName(method_name(i_m) + "_" + elevation_label(i_c))) = ...
            round(median(cut_relative(:, i_c, :, i_m), 3, "omitnan"), 2);
    end
end
disp(cut_summary)
writetable(cut_summary, fullfile(figure_dir, "GD_E_deficit_by_cut.csv"))
% Reference density per display layer, median over scans
reference_display = nan(n_display, n_method);
for i_d = 1:n_display
    in_display = layer_mid_m > display_edges_m(i_d, 1) & layer_mid_m < display_edges_m(i_d, 2);
    reference_display(i_d, :) = median(squeeze(mean(reference_density(1, in_display, 1, :, :), 2, "omitnan")), 1, "omitnan");
end
reference_table = array2table(round(reference_display, 1), VariableNames=cellstr(method_name), RowNames=cellstr(display_label));
disp("Density at 30-60 km (birds km^-3), median over scans:"); disp(reference_table)
writetable(reference_table, fullfile(figure_dir, "GD_E_reference.csv"), WriteRowNames=true)

%% Plot: deficit by layer, median and interquartile range over scans -------------------
method_colors = [0.1 0.4 0.8; 0 0 0; 0.8 0.2 0.15];
fig = bmm4d.plot.create_export_figure(3600, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, n_display, TileSpacing="tight", Padding="tight");
for i_d = 1:n_display
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_m = [3 1 2]
        values = squeeze(layer_relative(:, i_d, :, i_m));
        q = prctile(values, [25 50 75], 2);
        q(sum(~isnan(values), 2) < 10, :) = nan;
        valid = ~isnan(q(:, 2));
        fill(ax, [bin_km(valid) fliplr(bin_km(valid))], [q(valid, 1)' fliplr(q(valid, 3)')], method_colors(i_m, :), ...
            FaceAlpha=0.15, EdgeColor="none", HandleVisibility="off")
        plot(ax, bin_km(valid), q(valid, 2), "-o", Color=method_colors(i_m, :), LineWidth=2, MarkerSize=4, DisplayName=method_name(i_m))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.01 3])
    xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.01 0.03 0.1 0.3 1 3])
    xlabel(ax, "Range (km)"); ylabel(ax, "D(r), median and IQR over scans")
    title(ax, "Beam centre " + display_label(i_d) + " AGL", FontWeight="normal")
    if i_d == 1
        legend(ax, Location="southeast", FontSize=10)
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_E_deficit.png"), Resolution=150)
close(fig)

%% Plot: deficit by elevation class, all layers ------------------------------------------
cut_colors = parula(n_elevation + 1);
fig = bmm4d.plot.create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
for i_m = [3 1 2]
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_c = 1:n_elevation
        values = squeeze(cut_relative(:, i_c, :, i_m));
        med = median(values, 2, "omitnan");
        valid = sum(~isnan(values), 2) >= 10;
        plot(ax, bin_km(valid), med(valid), "-o", Color=cut_colors(i_c, :), LineWidth=2, MarkerSize=4, ...
            DisplayName=elevation_label(i_c))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.01 3])
    xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.01 0.03 0.1 0.3 1 3])
    xlabel(ax, "Range (km)"); ylabel(ax, "D(r), median over scans")
    title(ax, method_name(i_m), FontWeight="normal")
    if i_m == 2
        legend(ax, Location="southeast", FontSize=10)
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_E_deficit_by_cut.png"), Resolution=150)
close(fig)

%% Sensitivity to the cross-section spread ---------------------------------------------
% Estimate with log-sd 1 and 1.75 over the estimate with 1.25, per scan, range
% bin and elevation class (all layers), then the median over scans
cut_density = squeeze(sum(density_sum, 2) ./ sum(gate_count, 2));         % [bin cut scan method]
spread_ratio = cut_density(:, :, :, estimate_method(2:end)) ./ cut_density(:, :, :, estimate_method(1));
spread_ratio(isinf(spread_ratio)) = nan;
spread_summary = table(bin_km', VariableNames="range_km");
for i_v = 2:n_spread
    for i_c = 1:n_elevation
        spread_summary.(matlab.lang.makeValidName(sprintf("sd%g_over_sd1.25_", spread_values(i_v)) + elevation_label(i_c))) = ...
            round(median(spread_ratio(:, i_c, :, i_v - 1), 3, "omitnan"), 3);
    end
end
disp(spread_summary)
writetable(spread_summary, fullfile(figure_dir, "GD_E_spread_sensitivity.csv"))
% Same ratio for the 30-60 km reference at 300-500 m
reference_layers = layer_mid_m > 300 & layer_mid_m < 500;
disp("Reference density at 30-60 km, 300-500 m, ratio to log-sd 1.25 (median over scans):")
disp(round(median(squeeze(mean(reference_density(1, reference_layers, 1, :, estimate_method(2:end)), 2, "omitnan") ./ ...
    mean(reference_density(1, reference_layers, 1, :, estimate_method(1)), 2, "omitnan")), 1, "omitnan"), 3))

%% Plot: sensitivity to the cross-section spread ------------------------------------------
spread_colors = [0.35 0.35 0.35; 0.55 0.75 0.95; 0.05 0.2 0.5];
fig = bmm4d.plot.create_export_figure(3600, 1000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, n_elevation, TileSpacing="tight", Padding="tight");
for i_c = 1:n_elevation
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_v = 2:n_spread
        values = squeeze(spread_ratio(:, i_c, :, i_v - 1));
        q = prctile(values, [25 50 75], 2);
        q(sum(~isnan(values), 2) < 10, :) = nan;
        valid = ~isnan(q(:, 2));
        fill(ax, [bin_km(valid) fliplr(bin_km(valid))], [q(valid, 1)' fliplr(q(valid, 3)')], spread_colors(i_v, :), ...
            FaceAlpha=0.2, EdgeColor="none", HandleVisibility="off")
        plot(ax, bin_km(valid), q(valid, 2), "-o", Color=spread_colors(i_v, :), LineWidth=2, MarkerSize=4, ...
            DisplayName=sprintf("s_\\sigma = %g", spread_values(i_v)))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.5 2])
    xticks(ax, [2 5 10 20 60]); yticks(ax, [0.5 0.7 1 1.4 2])
    xlabel(ax, "Range (km)"); title(ax, elevation_label(i_c), FontWeight="normal")
    if i_c == 1
        ylabel(ax, "Estimate / estimate with s_\sigma = 1.25")
        legend(ax, Location="northeast", FontSize=10)
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_E_spread_sensitivity.png"), Resolution=150)
close(fig)
