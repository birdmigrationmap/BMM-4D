%% Gate density D: the near-range deficit on the sandbox's real sweeps, by method
% Level II sweeps of the sandbox case (28 radars, frames 2, 4, 6, cuts
% below 4.5 deg). For every gate, densities from a 2 km window
% (GD_WINDOW_SUM) by:
% - output mean: the radar output, gate by gate (no window);
% - window mean: linear mean of the output over the window;
% - step 1, sampling: the homogeneous-density estimator, output read as
%   exact (calibrated cross-section spread and beam volume, GD_C);
% - step 1, uncalibrated: same with the default model (log-sd 1, volume x1);
% - step 2, + CFP: step 1 with the radar's clutter filter information:
%   filtered gates between the output and output + CFP, point-filtered
%   gates at least the output;
% - our cleaned volumes: the current read_rd/clean_volume output (cache).
% The deficit D(r, layer) = density at r / density at 30-60 km in the same
% layer (beam-centre height above ground), per radar, then the median over
% radars (NEAR_RANGE_DEFICIT.md). Also maps of one sweep.
% Outputs: studies/gate_density/figures/GD_D_*.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");
c = prod_case("sandbox");
load(fullfile(c.data_dir, "case.mat"), "radar_selected")
manifest = load(fullfile(c.data_dir, "manifest.mat"));
n_radar = height(radar_selected);
frames = [2 4 6];

%% Settings ---------------------------------------------------------------
calibration = load(fullfile(figure_dir, "GD_C_calibration.mat"));
beam_calibrated = gd_beam(sigma_log_sd=calibration.calibrated_sigma_log_sd, ...
    width_deg=0.95 * sqrt(calibration.calibrated_volume_factor));
beam_default = gd_beam();
table_calibrated = gd_bird_sum_table(beam_calibrated, 600, 20000);
table_default = gd_bird_sum_table(beam_default, 600, 20000);
s_floor_cm2 = 0.1;
radius_m = 2000;
max_range_m = 80000;
max_density = 3000;
bin_edges_km = [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80];
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end)) / 2;
n_bin = numel(bin_km);
layer_edges_m = [100 300 500 800 1200 1700];
layer_label = compose("%d-%d m", layer_edges_m(1:end - 1)', layer_edges_m(2:end)');
n_layer = numel(layer_label);
method_name = ["output mean" "window mean" "step 1, sampling" "step 1, uncalibrated" "step 2, + CFP" "our cleaned volumes"];
n_method = numel(method_name);
earth_radius_m = 4 / 3 * 6371000;
example = struct("radar", "KDMX", "frame", 4, "elevation", 2.4);

%% Densities per gate, averaged per range bin, layer and radar ------------------------
example_sweep = [];
density_sum = zeros(n_bin, n_layer, n_radar, n_method);
gate_count = zeros(n_bin, n_layer, n_radar, n_method);
for i_r = 1:n_radar
    antenna_agl = radar_selected.height(i_r) - radar_selected.dem(i_r);
    for i_t = frames
        raw_file = manifest.raw_files(i_r, i_t);
        volume_file = manifest.volume_files(i_r, i_t);
        if ismissing(raw_file) || ~isfile(raw_file) || ismissing(volume_file) || ~isfile(volume_file)
            continue
        end
        level2 = read_level2_moments(raw_file, ["REF" "CFP"], 4.5);
        level2 = level2(arrayfun(@(s) ~isempty(s.CFP) && ~isempty(s.REF), level2));
        cleaned = load(volume_file).volume.dens.sweeps;
        for i_e = 1:numel(level2)
            n_gate = min([size(level2(i_e).REF, 1), size(level2(i_e).CFP, 1), find(level2(i_e).REF_range <= max_range_m, 1, "last")]);
            range_m = level2(i_e).REF_range(1:n_gate);
            dbz = level2(i_e).REF(1:n_gate, :);
            code = level2(i_e).CFP_code(1:n_gate, :);
            az_step_deg = 360 / size(dbz, 2);
            % Clutter and rain (> 3000 birds km^-3): left out of every method
            excluded = z_to_refl(idb(dbz)) / 11 > max_density;
            output = z_to_refl(idb(dbz)) / 11;
            output(isnan(dbz)) = 0;
            output(excluded) = nan;
            detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
            % The estimator variants
            [lo, hi] = gd_gate_bounds(dbz, zeros(size(code)), range_m, beam_calibrated, detection_dbz, s_floor_cm2 * calibration.calibrated_volume_factor);
            lo(excluded) = -Inf; hi(excluded) = Inf;
            step1 = gd_estimate_sweep(lo, hi, range_m, az_step_deg, beam_calibrated, table_calibrated, radius_m=radius_m);
            [lo, hi] = gd_gate_bounds(dbz, code, range_m, beam_calibrated, detection_dbz, s_floor_cm2 * calibration.calibrated_volume_factor);
            lo(excluded) = -Inf; hi(excluded) = Inf;
            step2 = gd_estimate_sweep(lo, hi, range_m, az_step_deg, beam_calibrated, table_calibrated, radius_m=radius_m);
            [lo, hi] = gd_gate_bounds(dbz, zeros(size(code)), range_m, beam_default, detection_dbz, s_floor_cm2);
            lo(excluded) = -Inf; hi(excluded) = Inf;
            step1_default = gd_estimate_sweep(lo, hi, range_m, az_step_deg, beam_default, table_default, radius_m=radius_m);
            window_mean = gd_window_mean(output, range_m, az_step_deg, radius_m);
            [~, i_clean] = min(abs([cleaned.elev] - level2(i_e).elev));
            clean_range = cleaned(i_clean).range(:);
            per_gate = {output, window_mean, step1.density, step1_default.density, step2.density};
            % Average per range bin, then assign the bin to its layer
            beam_agl = antenna_agl + bin_km' * 1000 * sind(level2(i_e).elev) + (bin_km' * 1000) .^ 2 / (2 * earth_radius_m);
            layer = discretize(beam_agl, layer_edges_m);
            bin_of_gate = discretize(range_m / 1000, bin_edges_km);
            bin_of_clean = discretize(clean_range / 1000, bin_edges_km);
            for k = find(~isnan(layer))'
                rows = bin_of_gate == k;
                if ~any(rows)
                    continue
                end
                for i_m = 1:5
                    values = per_gate{i_m}(rows, :);
                    density_sum(k, layer(k), i_r, i_m) = density_sum(k, layer(k), i_r, i_m) + sum(values, "all", "omitnan");
                    gate_count(k, layer(k), i_r, i_m) = gate_count(k, layer(k), i_r, i_m) + nnz(~isnan(values));
                end
                values = cleaned(i_clean).data(bin_of_clean == k, :);
                density_sum(k, layer(k), i_r, 6) = density_sum(k, layer(k), i_r, 6) + sum(values, "all", "omitnan");
                gate_count(k, layer(k), i_r, 6) = gate_count(k, layer(k), i_r, 6) + nnz(~isnan(values));
            end
            % Keep the example sweep for the maps
            if radar_selected.name(i_r) == example.radar && i_t == example.frame && abs(level2(i_e).elev - example.elevation) < 0.2
                example_sweep = struct("range_m", range_m, "az_step_deg", az_step_deg, "output", output, ...
                    "window_mean", window_mean, "step1", step1.density, "step2", step2.density, ...
                    "cfp_db", double(code) - 8, "code", code, "elev", level2(i_e).elev);
            end
        end
    end
    fprintf("[gate density] %s done\n", radar_selected.name(i_r))
end
save(fullfile(figure_dir, "GD_D_sandbox.mat"), "density_sum", "gate_count", "bin_km", "layer_label", "method_name", "example_sweep")

%% Deficit per radar: each layer over its own 30-60 km value ---------------------------
density_mean = density_sum ./ gate_count;                        % [bin layer radar method]
reference = bin_km >= 30 & bin_km <= 60;
reference_density = mean(density_mean(reference, :, :, :), 1, "omitnan");
relative = density_mean ./ reference_density;
relative(repmat(reference_density < 5, n_bin, 1, 1, 1)) = nan;
summary = table(bin_km', VariableNames="range_km");
for i_l = [2 4]
    for i_m = 1:n_method
        summary.(matlab.lang.makeValidName(layer_label(i_l) + "_" + method_name(i_m))) = ...
            round(median(relative(:, i_l, :, i_m), 3, "omitnan"), 2);
    end
end
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_D_deficit.csv"))
reference_table = array2table(round(squeeze(median(reference_density(1, [2 4], :, :), 3, "omitnan")), 1), ...
    VariableNames=cellstr(method_name), RowNames=cellstr(layer_label([2 4])));
disp("Density at 30-60 km (birds km^-3), median over radars:"); disp(reference_table)
writetable(reference_table, fullfile(figure_dir, "GD_D_reference.csv"), WriteRowNames=true)

%% Plot: deficit by method ---------------------------------------------------------------
method_colors = [0.1 0.4 0.8; 0.5 0.5 0.5; 0 0 0; 0.65 0.65 0.65; 0.55 0.3 0.75; 0.8 0.2 0.15];
method_style = ["-" "-" "-" "--" "-" "-"];
plot_layers = [2 3 4];
fig = create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, numel(plot_layers), TileSpacing="tight", Padding="tight");
title(t, "Near-range echo deficit by method (sandbox, 28 radars, 3 frames): density / 30-60 km density in the same layer", ...
    FontSize=14, Color="k")
for i_l = plot_layers
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_m = 1:n_method
        plot(ax, bin_km, median(relative(:, i_l, :, i_m), 3, "omitnan"), method_style(i_m) + "o", Color=method_colors(i_m, :), ...
            LineWidth=2, DisplayName=method_name(i_m))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.01 3])
    xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.01 0.03 0.1 0.3 1 3])
    xlabel(ax, "Range (km)"); ylabel(ax, "Median over radars")
    title(ax, "Beam centre " + layer_label(i_l) + " AGL", FontWeight="normal")
    if i_l == plot_layers(1)
        legend(ax, Location="southeast", FontSize=10)
    end
end
radar_colors = parula(n_radar + 4);
for i_m = [1 3 5]
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_r = 1:n_radar
        plot(ax, bin_km, relative(:, 2, i_r, i_m), Color=[radar_colors(i_r, :) 0.6], LineWidth=1)
    end
    plot(ax, bin_km, median(relative(:, 2, :, i_m), 3, "omitnan"), "k-", LineWidth=3)
    yline(ax, 1, "k:")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.01 10])
    xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.01 0.03 0.1 0.3 1 3 10])
    xlabel(ax, "Range (km)"); ylabel(ax, "Per radar, 300-500 m")
    title(ax, method_name(i_m) + " (black: median)", FontWeight="normal")
end
export_fixed_figure(fig, fullfile(figure_dir, "GD_D_deficit.png"), Resolution=150)
close(fig)

%% Plot: maps of the example sweep -----------------------------------------------------
assert(~isempty(example_sweep), "The example sweep was not found; change EXAMPLE.")
in_map = example_sweep.range_m <= 40000;
az_edges = 0:example_sweep.az_step_deg:360;
[az_grid, r_grid] = meshgrid(az_edges, example_sweep.range_m(in_map) / 1000);
x = r_grid .* sind(az_grid);
y = r_grid .* cosd(az_grid);
cfp_map = example_sweep.cfp_db;
cfp_map(example_sweep.code < 8) = nan;
panels = {example_sweep.output, example_sweep.window_mean, example_sweep.step1, example_sweep.step2, cfp_map};
panel_titles = ["Radar output (each gate)" "Window mean (2 km)" "Step 1: sampling estimator" "Step 2: + CFP bounds" ...
    "CFP: power removed by the radar (dB)"];
fig = create_export_figure(3600, 820, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(panels), TileSpacing="tight", Padding="tight");
title(t, sprintf("%s, %s UTC, %.1f deg, within 40 km", example.radar, string(manifest.times(example.frame), "HH:mm"), example_sweep.elev), ...
    FontSize=13, Color="k")
for i_p = 1:numel(panels)
    ax = nexttile(t);
    values = panels{i_p}(in_map, [1:end 1]);
    if i_p < numel(panels)
        values = log10(max(values, 0.1));
    end
    surf(ax, x, y, zeros(size(x)), values, EdgeColor="none")
    view(ax, 2); axis(ax, "equal", "tight"); set(ax, FontSize=11, Color=[0.92 0.92 0.92])
    if i_p < numel(panels)
        colormap(ax, parula); clim(ax, [0 2.5])
    else
        colormap(ax, flipud(gray(256))); clim(ax, [0 30])
        cb = colorbar(ax, Color="k"); cb.Label.String = "dB";
    end
    xlabel(ax, "x (km)"); ylabel(ax, "y (km)"); title(ax, panel_titles(i_p), FontWeight="normal")
end
ax = nexttile(t, 4);
cb = colorbar(ax, Color="k"); cb.Label.String = "log_{10} birds km^{-3}";
export_fixed_figure(fig, fullfile(figure_dir, "GD_D_maps.png"), Resolution=150)
close(fig)
