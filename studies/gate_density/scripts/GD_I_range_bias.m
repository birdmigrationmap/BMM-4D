%% Gate density I: range bias of the radar output, the earlier cleaning, a moving mean and the estimator
% On the reference scans of GD_J (124 scans, 2021-2025, 49 radars):
% 1. mean density per scan, sweep and range bin for the radar output (gates
%    above MAX_DENSITY left out, without which its bin means are spiky), the
%    cleaning of BMM-4D before the estimator (GD_OLD_CLEANING: outlier
%    filters, 2 km moving median, 3000 birds km^-3 cap and 25 km^2 patch
%    filter), the moving mean of the capped output over the estimator's
%    window (BMM4D.RADAR.GATE.WINDOW_SUM) and the estimator
%    (BMM4D.RADAR.GATE.ESTIMATE_SWEEP, no density cap), on all cuts up to 10
%    deg and gates to 80 km (scans in parallel, 4 workers);
% 2. range bias: a vertical profile per scan and one bias per range bin,
%    fitted to all sweeps with the beam's vertical extent (GD_FIT_RANGE_BIAS),
%    for all scans and per density class of the scan, with a bootstrap band
%    over scans for the earlier cleaning and the estimator. An unbiased estimator
%    has b(d) = 1 at every distance d.
% Outputs: studies/gate_density/figures/GD_I_range_bias.png, .csv,
% GD_I_range_bias_bootstrap.csv and GD_I_selected_scans.mat.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
beam = bmm4d.radar.gate.beam();                      % cross-section log-sd 1.25 (assumed)
table_s = bmm4d.radar.gate.bird_sum_table(beam, 2000, 20000);
s_floor_cm2 = 0.1;
radius_m = 1000;
max_range_m = 80000;
max_elevation_deg = 10;
max_density = 3000;
bin_edges_km = [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80];
% Density classes of the scans: the scan's own mean output at 20-60 km on the lowest
% cut (GD_J), which sets how many birds the near-range gates hold
class_edges = [0 20 100 Inf];                        % birds km^-3
class_name = ["below 20" "20-100" "above 100"];
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end)) / 2;
n_bin = numel(bin_km);
% Radar output: gates above MAX_DENSITY left out, without which its bin means are spiky;
% moving mean of the same capped output over the estimator's window
method_name = ["radar output" "BMM-4D before (outlier filters)" "moving mean" "estimator"];
n_method = numel(method_name);

%% 1. The reference scans (GD_J): density per scan, sweep and range bin --------------------
scans = load(fullfile(bmm4d.util.local_data_dir(), "gate_density", "post2020_scans", "scans.mat")).scans;
n_scan = height(scans);
n_sweep_max = 20;
per_scan = cell(n_scan, 1);
pool = gcp("nocreate");
if isempty(pool) || pool.NumWorkers > 4
    delete(pool);
    parpool("Processes", 4);
end
parfor i_s = 1:n_scan
    level2 = bmm4d.radar.read_level2_moments(scans.raw_file(i_s), "REF", max_elevation_deg);
    level2 = level2(arrayfun(@(s) ~isempty(s.REF), level2));
    % Split cuts repeat an elevation: keep the first (surveillance) cut
    [~, keep] = unique(round([level2.elev] * 2.5) / 2.5, "stable");
    level2 = level2(keep);
    out = struct("density_sum", zeros(n_bin, n_sweep_max, n_method), "gate_count", zeros(n_bin, n_sweep_max, n_method), ...
        "elevation", nan(n_sweep_max, 1));
    for i_e = 1:min(numel(level2), n_sweep_max)
        n_gate = min(size(level2(i_e).REF, 1), find(level2(i_e).REF_range <= max_range_m, 1, "last"));
        range_m = level2(i_e).REF_range(1:n_gate);
        dbz = level2(i_e).REF(1:n_gate, :);
        az_step_deg = 360 / size(dbz, 2);
        output = z_to_refl(idb(dbz)) / 11;
        output(isnan(dbz)) = 0;
        output_capped = output;
        output_capped(output > max_density) = nan;
        old_cleaning = gd_old_cleaning(dbz, range_m, double(level2(i_e).az), level2(i_e).elev);
        window_mean = gd_window_mean(output_capped, range_m, az_step_deg, radius_m);
        detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        [lo, hi] = bmm4d.radar.gate.bounds(dbz, isnan(dbz), range_m, beam, detection_dbz, s_floor_cm2);
        estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, az_step_deg, beam, table_s, radius_m=radius_m);
        per_gate = {output_capped, old_cleaning, window_mean, estimate.density};
        bin_of_gate = discretize(range_m / 1000, bin_edges_km);
        for i_m = 1:n_method
            values = per_gate{i_m};
            valid = ~isnan(values) & ~isnan(bin_of_gate);
            rows = any(valid, 2);
            out.density_sum(:, i_e, i_m) = accumarray(bin_of_gate(rows), sum(values(rows, :) .* valid(rows, :), 2, "omitnan"), [n_bin 1]);
            out.gate_count(:, i_e, i_m) = accumarray(bin_of_gate(rows), sum(valid(rows, :), 2), [n_bin 1]);
        end
        out.elevation(i_e) = level2(i_e).elev;
    end
    per_scan{i_s} = out;
    fprintf("[range bias] %d/%d %s done\n", i_s, n_scan, scans.station(i_s))
end
parts = cellfun(@(o) o.density_sum, per_scan, UniformOutput=false);
density_sum = permute(cat(4, parts{:}), [1 2 4 3]);
parts = cellfun(@(o) o.gate_count, per_scan, UniformOutput=false);
gate_count = permute(cat(4, parts{:}), [1 2 4 3]);
parts = cellfun(@(o) o.elevation, per_scan, UniformOutput=false);
sweep_elevation = cat(2, parts{:});
save(fullfile(figure_dir, "GD_I_selected_scans.mat"), "density_sum", "gate_count", "sweep_elevation", "bin_edges_km", ...
    "bin_km", "method_name", "scans")

%% 2. Range bias with the beam's vertical extent ---------------------------------------
% Four groups: all scans and the three density classes
scans.density_class = discretize(scans.bird_density, class_edges);
[i_bin, i_sweep, i_scan] = ndgrid(1:n_bin, 1:n_sweep_max, 1:n_scan);
range_obs_m = bin_km(i_bin(:))' * 1000;
elevation_obs = sweep_elevation(sub2ind(size(sweep_elevation), i_sweep(:), i_scan(:)));
antenna_obs = scans.antenna_agl(i_scan(:));
subsets = [{true(n_scan, 1)} arrayfun(@(k) scans.density_class == k, 1:numel(class_name), UniformOutput=false)];
subset_name = [sprintf("All %d scans", n_scan) arrayfun(@(k) sprintf("%d scans, %s birds km^{-3}", ...
    nnz(scans.density_class == k), class_name(k)), 1:numel(class_name))];
bias = nan(n_bin, n_method, numel(subsets));
fits = cell(n_method, numel(subsets));
for i_sub = 1:numel(subsets)
    in_subset = subsets{i_sub}(i_scan(:));
    for i_m = 1:n_method
        mean_density = reshape(density_sum(:, :, :, i_m) ./ gate_count(:, :, :, i_m), [], 1);
        count = reshape(gate_count(:, :, :, i_m), [], 1) .* in_subset;
        fits{i_m, i_sub} = gd_fit_range_bias(mean_density, count, range_obs_m, elevation_obs, i_scan(:), antenna_obs, beam, ...
            bin_edges_km=bin_edges_km);
        bias(:, i_m, i_sub) = fits{i_m, i_sub}.bias;
    end
end
bias_table = table(bin_km', VariableNames="range_km");
for i_sub = 1:numel(subsets)
    for i_m = 1:n_method
        bias_table.(matlab.lang.makeValidName(subset_name(i_sub) + "_" + method_name(i_m))) = round(bias(:, i_m, i_sub), 2);
    end
end
disp(bias_table)
writetable(bias_table, fullfile(figure_dir, "GD_I_range_bias.csv"))
% Spread over scans: refit after resampling whole scans with replacement (moving
% mean and estimator, all scans; the same resamples for both)
n_boot = 200;
boot_methods = [2 4];
rng(3)
times_picked = zeros(n_scan, n_boot);
for i_b = 1:n_boot
    times_picked(:, i_b) = accumarray(randi(n_scan, n_scan, 1), 1, [n_scan 1]);
end
boot_table = table(bin_km', VariableNames="range_km");
boot_band = nan(n_bin, 2, n_method);
for i_m = boot_methods
    mean_density = reshape(density_sum(:, :, :, i_m) ./ gate_count(:, :, :, i_m), [], 1);
    count = reshape(gate_count(:, :, :, i_m), [], 1);
    bias_boot = nan(n_bin, n_boot);
    for i_b = 1:n_boot
        fit = gd_fit_range_bias(mean_density, count .* times_picked(i_scan(:), i_b), range_obs_m, elevation_obs, i_scan(:), ...
            antenna_obs, beam, bin_edges_km=bin_edges_km);
        bias_boot(:, i_b) = fit.bias;
    end
    boot_band(:, :, i_m) = prctile(bias_boot, [5 95], 2);
    name = matlab.lang.makeValidName(method_name(i_m));
    boot_table.(name) = round(bias(:, i_m, 1), 2);
    boot_table.(name + "_q05") = round(boot_band(:, 1, i_m), 2);
    boot_table.(name + "_q95") = round(boot_band(:, 2, i_m), 2);
end
writetable(boot_table, fullfile(figure_dir, "GD_I_range_bias_bootstrap.csv"))
% Check of the profile model: observed over predicted density without the bias, per
% elevation class (the GD_E classes), estimator, all scans. Flat in range and
% equal across classes if the beam's vertical extent explains the elevation offset of GD_E.
% Beams centred above 3500 m are left out: the profile is not constrained there.
elevation_edges = [0 1.1 2.1 3.6 5.75 max_elevation_deg + 0.5];
elevation_label = ["0.5-0.9°" "1.3-1.8°" "2.4-3.4°" "4.0-5.1°" "6.0-10°"];
fit = fits{4, 1};
i_class = discretize(fit.elevation_deg, elevation_edges);
cut_ratio = nan(n_bin, numel(elevation_label));
for i_c = 1:numel(elevation_label)
    in_class = i_class == i_c & fit.height_m < 3500;
    observed = accumarray(fit.i_bin(in_class), fit.weight(in_class) .* fit.density(in_class), [n_bin 1]);
    predicted = accumarray(fit.i_bin(in_class), fit.weight(in_class) .* fit.predicted(in_class) ./ fit.bias(fit.i_bin(in_class)), [n_bin 1]);
    enough = accumarray(fit.i_bin(in_class), 1, [n_bin 1]) >= 20;
    cut_ratio(enough, i_c) = observed(enough) ./ predicted(enough);
end
disp(array2table(round([bin_km' cut_ratio], 2), VariableNames=["range_km" elevation_label]))

%% Plot: range bias ----------------------------------------------------------------------
method_colors = [0.55 0.55 0.55; 0.75 0.15 0.15; 0.85 0.5 0.1; 0 0 0];   % as GD_A: output grey, earlier cleaning red
fig = bmm4d.plot.create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(subsets), TileSpacing="tight", Padding="tight");
for i_sub = 1:numel(subsets)
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    if i_sub == 1
        for i_m = boot_methods
            fill(ax, [bin_km fliplr(bin_km)], [boot_band(:, 1, i_m)' fliplr(boot_band(:, 2, i_m)')], method_colors(i_m, :), ...
                FaceAlpha=0.2, EdgeColor="none", DisplayName=method_name(i_m) + ", 5-95% over scans")
        end
    end
    for i_m = 1:n_method
        valid = ~isnan(bias(:, i_m, i_sub));
        plot(ax, bin_km(valid), bias(valid, i_m, i_sub), "-o", Color=method_colors(i_m, :), LineWidth=2, MarkerSize=4, DisplayName=method_name(i_m))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0.03 1.5]); xticks(ax, [2 5 10 20 40 80])
    yticks(ax, [0.03 0.1 0.3 1])
    xlabel(ax, "Distance to the radar (km)"); title(ax, subset_name(i_sub), FontWeight="normal")
    if i_sub == 1
        ylabel(ax, "Bias with distance b(d)")
        lgd = legend(ax, FontSize=12, Orientation="horizontal");
        lgd.Layout.Tile = "south";
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_I_range_bias.png"), Resolution=150)
close(fig)
