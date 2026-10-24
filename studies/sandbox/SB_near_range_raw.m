%% Sandbox: near-range echo deficit in the radar data themselves
% Model-free definition of the inconsistency: within one volume, gates of
% different elevations that sample the same height above ground at
% different ranges should give about the same bird density if birds are
% horizontally uniform over ~50 km. For each radar, gates are binned by
% beam-centre height above the radar's ground (layers) and by range, the
% 6 frames pooled; each layer is normalized by its own mean at 30-60 km.
% A curve below 1 at short range is the near-range echo deficit.
% Three stages are compared:
% - radar output: REF as distributed (after the radar's clutter filter);
% - unfiltered: REF + CFP, the clutter filter power removed (dB), i.e. the
%   reflectivity before the radar's clutter filter (ground clutter
%   included where there is some);
% - cleaned: our cleaned volume cache (MistNet, outliers, azimuthal shape
%   correction, smoothing, max density, patch filter, thinning).
% Only the surveillance cuts (those with CFP) are used for the first two.
% No-echo gates count as 0 density; gates above 3000 birds km^-3 (~40 dBZ:
% clutter, rain) are dropped. Gates with no echo after the filter but
% power removed have an unknown unfiltered value and are left out of the
% unfiltered stage. A second split is by azimuth relative to the bird
% flight direction (mean EP-VVP wind of the cleaned volume): along the
% flight (|cos| > 0.7) or across it (|cos| < 0.3), where the radial
% velocity is near 0 and a clutter filter removes bird echo first.
% Outputs: studies/sandbox/figures/near_range_raw*.png and .csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
manifest = load(fullfile(c.data_dir, "manifest.mat"));
load(fullfile(c.data_dir, "case.mat"), "radar_selected")
[n_radar, n_time] = size(manifest.raw_files);

%% Bins ------------------------------------------------------------------
% 2 km bins below 20 km: cleaned volumes are thinned to 1 km gates
range_edges_km = [2:2:20, 22.5:2.5:50, 55:5:80];
range_center_km = (range_edges_km(1:end - 1) + range_edges_km(2:end)) / 2;
layer_edges_m = [100 300 500 800 1200 1700];
layer_label = compose("%d-%d m", layer_edges_m(1:end - 1)', layer_edges_m(2:end)');
n_range = numel(range_center_km);
n_layer = numel(layer_label);
reference = range_center_km >= 30 & range_center_km <= 60;
max_density = 3000;
earth_radius_m = 4 / 3 * 6371000;
stage_name = ["Radar output (REF)" "Unfiltered (REF + CFP)" "Cleaned volumes"];
n_stage = numel(stage_name);
direction_name = ["along flight" "across flight"];

%% Accumulate per radar: density sums and counts per range, layer, stage ---
% and, for the radar output, per flight-direction sector; CFP statistics
% per range and elevation group.
density_sum = zeros(n_range, n_layer, n_radar, n_stage);
gate_count = zeros(n_range, n_layer, n_radar, n_stage);
direction_sum = zeros(n_range, n_layer, n_radar, 2, 2);     % [.. direction stage(output, unfiltered)]
direction_count = zeros(n_range, n_layer, n_radar, 2, 2);
elevation_edges = [0 0.7 1.1 1.5 2.1 2.8 3.5 4.5 5.5 7];
elevation_label = ["0.5°" "0.9°" "1.3°" "1.8°" "2.4°" "3.1°" "4.0°" "5.1°" "6.4°"];
n_elevation = numel(elevation_label);
cfp_filtered = zeros(n_range, n_elevation, n_radar);         % gates with power removed (code >= 8)
cfp_count = zeros(n_range, n_elevation, n_radar);
cfp_db_sum = zeros(n_range, n_elevation, n_radar);           % sum of dB removed over filtered gates
% Reflectivity distribution of the radar output in the 300-500 m layer,
% near (5-10 km) and in the reference ring (30-60 km), and the share of
% gates flagged by the point clutter filter (CFP code 1) per range
dbz_edges = -30:2.5:50;
ring_from_m = [5000 30000];
ring_to_m = [10000 60000];
dbz_histogram = zeros(numel(dbz_edges) - 1, 2);
ring_gate_count = zeros(1, 2);
point_filtered = zeros(n_range, 1);
point_count = zeros(n_range, 1);
vcp = nan(n_radar, n_time);
for i_r = 1:n_radar
    antenna_agl = radar_selected.height(i_r) - radar_selected.dem(i_r);
    for i_t = 1:n_time
        raw_file = manifest.raw_files(i_r, i_t);
        volume_file = manifest.volume_files(i_r, i_t);
        if ismissing(raw_file) || ~isfile(raw_file) || ismissing(volume_file) || ~isfile(volume_file)
            continue
        end
        volume = load(volume_file).volume;
        uv = mean(volume.uv, 1, "omitnan");
        flight_az = mod(atan2d(uv(1), uv(2)), 360);           % direction of travel, clockwise from north
        level2 = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], 7);
        level2 = level2(arrayfun(@(s) ~isempty(s.CFP) && ~isempty(s.REF), level2));
        % Sweeps of each stage as {elev, range (m), az, density}
        sweeps = {};
        for i_e = 1:numel(level2)
            n_gate = min(size(level2(i_e).REF, 1), size(level2(i_e).CFP, 1));
            ref_dbz = level2(i_e).REF(1:n_gate, :);
            code = double(level2(i_e).CFP_code(1:n_gate, :));
            range_m = level2(i_e).REF_range(1:n_gate);
            output = z_to_refl(idb(ref_dbz)) / 11;
            output(isnan(ref_dbz)) = 0;
            unfiltered_dbz = ref_dbz + max(code - 8, 0);
            unfiltered = z_to_refl(idb(unfiltered_dbz)) / 11;
            unfiltered(isnan(ref_dbz) & code < 8) = 0;
            unfiltered(isnan(ref_dbz) & code >= 8) = nan;
            sweeps(end + 1, :) = {level2(i_e).elev, range_m, level2(i_e).az, output, 1}; %#ok<SAGROW>
            beam_agl = antenna_agl + range_m * sind(level2(i_e).elev) + range_m .^ 2 / (2 * earth_radius_m);
            in_layer = beam_agl >= 300 & beam_agl < 500;
            for i_h = 1:2
                ring = in_layer & range_m >= ring_from_m(i_h) & range_m < ring_to_m(i_h);
                values = ref_dbz(ring, :);
                dbz_histogram(:, i_h) = dbz_histogram(:, i_h) + histcounts(values(:), dbz_edges)';
                ring_gate_count(i_h) = ring_gate_count(i_h) + numel(values);
            end
            point_bin = discretize(range_m / 1000, range_edges_km);
            rows_point = ~isnan(point_bin);
            point_filtered = point_filtered + accumarray(point_bin(rows_point), sum(code(rows_point, :) == 1, 2), [n_range 1]);
            point_count = point_count + accumarray(point_bin(rows_point), size(code, 2), [n_range 1]);
            sweeps(end + 1, :) = {level2(i_e).elev, range_m, level2(i_e).az, unfiltered, 2}; %#ok<SAGROW>
            % CFP statistics by range and elevation
            range_bin = discretize(range_m / 1000, range_edges_km);
            elevation_bin = discretize(level2(i_e).elev, elevation_edges);
            rows = ~isnan(range_bin);
            if ~isnan(elevation_bin)
                filtered = code(rows, :) >= 8;
                cfp_filtered(:, elevation_bin, i_r) = cfp_filtered(:, elevation_bin, i_r) + accumarray(range_bin(rows), sum(filtered, 2), [n_range 1]);
                cfp_count(:, elevation_bin, i_r) = cfp_count(:, elevation_bin, i_r) + accumarray(range_bin(rows), size(code, 2), [n_range 1]);
                cfp_db_sum(:, elevation_bin, i_r) = cfp_db_sum(:, elevation_bin, i_r) + ...
                    accumarray(range_bin(rows), sum((code(rows, :) - 8) .* filtered, 2), [n_range 1]);
            end
        end
        for i_e = 1:numel(volume.dens.sweeps)
            cleaned = volume.dens.sweeps(i_e);
            sweeps(end + 1, :) = {cleaned.elev, cleaned.range(:), cleaned.az(:)', cleaned.data, 3}; %#ok<SAGROW>
        end
        % Bin every sweep by range and beam-centre height above ground
        for i_s = 1:size(sweeps, 1)
            [elev, range_m, az, density, stage] = sweeps{i_s, :};
            beam_agl = antenna_agl + range_m * sind(elev) + range_m .^ 2 / (2 * earth_radius_m);
            range_bin = discretize(range_m / 1000, range_edges_km);
            layer = discretize(beam_agl, layer_edges_m);
            rows = ~isnan(range_bin) & ~isnan(layer);
            data = density(rows, :);
            valid = ~isnan(data) & data <= max_density;
            data(~valid) = 0;
            key = [range_bin(rows) layer(rows)];
            density_sum(:, :, i_r, stage) = density_sum(:, :, i_r, stage) + accumarray(key, sum(data, 2), [n_range n_layer]);
            gate_count(:, :, i_r, stage) = gate_count(:, :, i_r, stage) + accumarray(key, sum(valid, 2), [n_range n_layer]);
            if stage <= 2
                alignment = abs(cosd(az - flight_az));
                sectors = {alignment > 0.7, alignment < 0.3};
                for i_d = 1:2
                    in_sector = valid & sectors{i_d};
                    direction_sum(:, :, i_r, i_d, stage) = direction_sum(:, :, i_r, i_d, stage) + ...
                        accumarray(key, sum(data .* in_sector, 2), [n_range n_layer]);
                    direction_count(:, :, i_r, i_d, stage) = direction_count(:, :, i_r, i_d, stage) + ...
                        accumarray(key, sum(in_sector, 2), [n_range n_layer]);
                end
            end
        end
    end
    fprintf("[raw near range] %s done\n", radar_selected.name(i_r))
end

%% Relative density per radar: each layer over its own 30-60 km mean ------
% Layers with < 5 birds km^-3 in the reference ring are left out.
mean_density = density_sum ./ gate_count;                      % [range layer radar stage]
reference_density = mean(mean_density(reference, :, :, :), 1, "omitnan");
relative = mean_density ./ reference_density;
relative(repmat(reference_density < 5, n_range, 1, 1, 1)) = nan;
direction_density = direction_sum ./ direction_count;          % [range layer radar direction stage]
direction_relative = direction_density ./ mean(direction_density(reference, :, :, :, :), 1, "omitnan");
direction_relative(repmat(reference_density(:, :, :, 1) < 5, n_range, 1, 1, 2, 2)) = nan;
cfp_share = cfp_filtered ./ cfp_count;
cfp_mean_db = cfp_db_sum ./ cfp_filtered;

%% Summary tables: median over radars ------------------------------------------
summary = table(range_center_km(:), VariableNames="range_km");
for stage = 1:n_stage
    for i_l = [2 4]
        summary.(sprintf("stage%d_%s", stage, matlab.lang.makeValidName(layer_label(i_l)))) = ...
            round(median(relative(:, i_l, :, stage), 3, "omitnan"), 2);
    end
end
for i_d = 1:2
    summary.("output_" + matlab.lang.makeValidName(direction_name(i_d)) + "_300_500m") = ...
        round(median(direction_relative(:, 2, :, i_d, 1), 3, "omitnan"), 2);
end
disp(summary)
writetable(summary, fullfile(figure_dir, "near_range_raw.csv"))
cfp_table = array2table(round(median(cfp_mean_db, 3, "omitnan"), 1), VariableNames=cellstr(elevation_label));
cfp_table = addvars(cfp_table, range_center_km(:), Before=1, NewVariableNames="range_km");
disp("Mean CFP (dB removed, filtered gates), median over radars:"); disp(cfp_table)
writetable(cfp_table, fullfile(figure_dir, "near_range_cfp.csv"))

%% Plot: relative density per radar, three stages ---------------------------------
radar_colors = parula(n_radar + 4);
plot_layers = [1 2 3 4];
fig = bmm4d.plot.create_export_figure(3600, 2700, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, n_stage, numel(plot_layers), TileSpacing="tight", Padding="tight");
title(t, "Density at the same height above ground, relative to 30-60 km, per radar (6 frames pooled; black: median)", ...
    FontSize=14, Color="k")
for stage = 1:n_stage
    for i_l = plot_layers
        ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
        for i_r = 1:n_radar
            plot(ax, range_center_km, relative(:, i_l, i_r, stage), Color=[radar_colors(i_r, :) 0.6], LineWidth=1)
        end
        plot(ax, range_center_km, median(relative(:, i_l, :, stage), 3, "omitnan"), "k-", LineWidth=3)
        yline(ax, 1, "k--")
        set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.03 10])
        xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.03 0.1 0.3 1 3 10])
        xlabel(ax, "Range (km)")
        if i_l == plot_layers(1)
            ylabel(ax, stage_name(stage) + newline + "density / 30-60 km density")
        end
        title(ax, "Beam centre " + layer_label(i_l) + " AGL", FontWeight="normal")
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "near_range_raw.png"), Resolution=150)
close(fig)

%% Plot: clutter filter power removed, and the flight-direction split -------------
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 2, TileSpacing="tight", Padding="tight");
elevation_colors = parula(n_elevation + 2);
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_e = 1:n_elevation
    plot(ax, range_center_km, median(cfp_mean_db(:, i_e, :), 3, "omitnan"), "-", Color=elevation_colors(i_e, :), ...
        LineWidth=2, DisplayName=elevation_label(i_e))
end
set(ax, XScale="log", FontSize=12); xlim(ax, [2 80]); xticks(ax, [2 5 10 20 30 60])
xlabel(ax, "Range (km)"); ylabel(ax, "Mean power removed on filtered gates (dB)")
legend(ax, Location="northeast", NumColumns=2); title(ax, "Clutter filter power removed (median over radars)", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_e = 1:n_elevation
    plot(ax, range_center_km, median(cfp_share(:, i_e, :), 3, "omitnan"), "-", Color=elevation_colors(i_e, :), LineWidth=2)
end
set(ax, XScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0 1]); xticks(ax, [2 5 10 20 30 60])
plot(ax, range_center_km, point_filtered ./ point_count, "k--", LineWidth=2)
xlabel(ax, "Range (km)"); ylabel(ax, "Share of gates")
title(ax, "Gates with power removed, per elevation (dashed: point clutter filter, all cuts)", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
line_style = ["-" "--"];
direction_colors = [0.8 0.3 0.1; 0.1 0.4 0.8];
for stage = 1:2
    for i_d = 1:2
        plot(ax, range_center_km, median(direction_relative(:, 2, :, i_d, stage), 3, "omitnan"), line_style(stage), ...
            Color=direction_colors(i_d, :), LineWidth=2, DisplayName=stage_name(stage) + ", " + direction_name(i_d))
    end
end
yline(ax, 1, "k--", HandleVisibility="off")
set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.03 10]); xticks(ax, [2 5 10 20 30 60])
xlabel(ax, "Range (km)"); ylabel(ax, "Density / 30-60 km density (300-500 m AGL)")
legend(ax, Location="southeast"); title(ax, "Along vs across the flight direction", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
dbz_center = (dbz_edges(1:end - 1) + dbz_edges(2:end)) / 2;
ring_names = ["5-10 km" "30-60 km"];
for i_h = 1:2
    plot(ax, dbz_center, dbz_histogram(:, i_h) / ring_gate_count(i_h), "-o", Color=direction_colors(i_h, :), LineWidth=2, ...
        DisplayName=sprintf("%s (no echo %.0f%%)", ring_names(i_h), 100 * (1 - sum(dbz_histogram(:, i_h)) / ring_gate_count(i_h))))
end
set(ax, FontSize=12); xlabel(ax, "Reflectivity, radar output (dBZ)"); ylabel(ax, "Share of gates")
legend(ax, Location="northwest"); title(ax, "Reflectivity distribution, beam centre 300-500 m AGL", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "near_range_cfp.png"), Resolution=150)
close(fig)
