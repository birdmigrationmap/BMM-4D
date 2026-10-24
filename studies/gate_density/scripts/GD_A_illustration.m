%% Gate density A: what happens near a radar, on simulated sweeps
% Synthetic illustration of the near-range problem and of the fix:
% 1. why gates near a radar are spiky: birds per gate, empty gates, and the
%    reflectivity of one bird, by range;
% 2. what the gates look like (value distribution by range);
% 3. maps of one simulated sweep: true density, radar output, our current
%    cleaning (read_rd), the window mean and the homogeneous-density
%    estimator, all with the same 2 km window;
% 4. recovered / true density by range for 5, 20, 50 and 200 birds km^-3,
%    and the per-gate error of the window mean and of the estimator.
% Outputs: studies/gate_density/figures/GD_A_*.png and .csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
beam = gd_beam();                                    % 0.95 deg beam, 250 m gates, sigma 11 cm^2, log-sd 1
table_s = gd_bird_sum_table(beam, 300, 20000);
range_m = (2:0.25:60)' * 1000;
az_step_deg = 0.5;
radius_m = 2000;                                     % window half-width (range) and half-arc
detection_dbz = -28 + 20 * log10(range_m / 5000);    % about the clear-air sensitivity seen near range
background_dbz = -15;
s_floor_cm2 = 0.1;                                   % 1% of a mean bird
bin_edges_km = [2 3 4 5 6 7 8 10 12 15 20 25 30 40 50 60];
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end))' / 2;
bin_of = discretize(range_m / 1000, bin_edges_km);
bin_mean = @(x) accumarray(bin_of, mean(x, 2, "omitnan"), [numel(bin_km) 1], @mean, nan);
bin_rmse = @(x) sqrt(accumarray(bin_of, mean((x - 1) .^ 2, 2, "omitnan"), [numel(bin_km) 1], @mean, nan));
density_from_dbz = @(dbz) z_to_refl(idb(dbz)) / 11;

%% 1. Why the gates are spiky ---------------------------------------------------
plot_range_m = logspace(log10(1000), log10(100000), 200)';
example_density = [5 20 50 200];
birds_per_gate = example_density .* beam.volume_sample_km3(plot_range_m);
one_bird_density = beam.sigma_mean_cm2 ./ beam.volume_eff_km3(plot_range_m) / 11;   % one bird on the beam axis
fig = create_export_figure(3600, 1150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
density_colors = [0.55 0.75 0.95; 0.25 0.5 0.85; 0.1 0.3 0.6; 0 0.1 0.3];
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_d = 1:numel(example_density)
    plot(ax, plot_range_m / 1000, birds_per_gate(:, i_d), Color=density_colors(i_d, :), LineWidth=2, ...
        DisplayName=example_density(i_d) + " birds km^{-3}")
end
yline(ax, 1, "k--", HandleVisibility="off")
set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [1 100]); ylim(ax, [1e-3 1e3]); xticks(ax, [1 2 5 10 20 50 100])
xlabel(ax, "Range (km)"); ylabel(ax, "Birds per gate (mean)"); legend(ax, Location="northwest")
title(ax, "Gate volume grows with range^2", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_d = 1:numel(example_density)
    plot(ax, plot_range_m / 1000, exp(-birds_per_gate(:, i_d)), Color=density_colors(i_d, :), LineWidth=2)
end
set(ax, XScale="log", FontSize=12); xlim(ax, [1 100]); ylim(ax, [0 1]); xticks(ax, [1 2 5 10 20 50 100])
xlabel(ax, "Range (km)"); ylabel(ax, "Share of gates without a bird"); title(ax, "Most near-range gates are empty", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
plot(ax, plot_range_m / 1000, one_bird_density, "k-", LineWidth=2, DisplayName="one bird on the beam axis")
for i_d = 1:numel(example_density)
    yline(ax, example_density(i_d), "--", Color=density_colors(i_d, :), LineWidth=1.5, HandleVisibility="off")
end
set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [1 100]); ylim(ax, [1 1e5]); xticks(ax, [1 2 5 10 20 50 100])
xlabel(ax, "Range (km)"); ylabel(ax, "Gate density it produces (birds km^{-3})"); legend(ax, Location="northeast")
title(ax, "... and the few occupied gates are spikes", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "GD_A_why_spikes.png"), Resolution=150)
close(fig)

%% 2 and 3. One simulated sweep: values by range, and maps ----------------------------
rng(11)
true_density = gd_random_field(range_m, az_step_deg, 20, 0.4, 15000);
eta = gd_simulate_sweep(true_density, range_m, az_step_deg, beam);
dbz = gd_observe(eta, detection_dbz, background_dbz);
output = density_from_dbz(dbz);
output(isnan(dbz)) = 0;
cleaned = gd_read_rd_filters(dbz, range_m, az_step_deg);
window_mean = gd_window_mean(output, range_m, az_step_deg, radius_m);
[lo, hi] = gd_gate_bounds(dbz, zeros(size(dbz)), range_m, beam, detection_dbz, s_floor_cm2);
estimate = gd_estimate_sweep(lo, hi, range_m, az_step_deg, beam, table_s, radius_m=radius_m);
truth_window = gd_window_mean(true_density, range_m, az_step_deg, radius_m);   % the target: true density averaged over the window

dbz_edges = -40:2.5:50;
dbz_center = (dbz_edges(1:end - 1) + dbz_edges(2:end)) / 2;
fig = create_export_figure(2400, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
ax = axes(fig); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
ring_plot = [3 5; 5 10; 15 20; 30 60];
ring_colors = [0.8 0.3 0.1; 0.9 0.6 0.1; 0.3 0.6 0.3; 0.1 0.4 0.8];
for i_p = 1:size(ring_plot, 1)
    rows = range_m >= ring_plot(i_p, 1) * 1000 & range_m < ring_plot(i_p, 2) * 1000;
    values = dbz(rows, :);
    plot(ax, dbz_center, histcounts(values(:), dbz_edges) / numel(values), "-o", Color=ring_colors(i_p, :), LineWidth=2, ...
        DisplayName=sprintf("%g-%g km", ring_plot(i_p, :)))
end
set(ax, FontSize=13); xlabel(ax, "Gate reflectivity (dBZ)"); ylabel(ax, "Share of gates"); legend(ax, Location="northwest")
title(ax, "Simulated sweep, ~20 birds km^{-3}: background (left mode) and gates holding birds (right)", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "GD_A_gate_values.png"), Resolution=150)
close(fig)

map_range = range_m <= 30000;
% Last azimuth repeated at 360 deg so the map closes at north
[az_grid, r_grid] = meshgrid(0:az_step_deg:360, range_m(map_range) / 1000);
x = r_grid .* sind(az_grid);
y = r_grid .* cosd(az_grid);
panels = {true_density, output, cleaned, window_mean, estimate.density};
panel_titles = ["True density" "Radar output (each gate)" "Our current cleaning (read_rd)" ...
    "Window mean (2 km)" "Estimator (2 km)"];
fig = create_export_figure(3600, 820, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(panels), TileSpacing="tight", Padding="tight");
for i_p = 1:numel(panels)
    ax = nexttile(t);
    values = panels{i_p}(map_range, [1:end 1]);
    surf(ax, x, y, zeros(size(x)), log10(max(values, 0.1)), EdgeColor="none")
    view(ax, 2); axis(ax, "equal", "tight"); set(ax, FontSize=11, Color=[0.92 0.92 0.92]); colormap(ax, parula); clim(ax, [0 2.5])
    xlabel(ax, "x (km)"); ylabel(ax, "y (km)"); title(ax, panel_titles(i_p), FontWeight="normal", Interpreter="none")
end
cb = colorbar(ax, Color="k"); cb.Layout.Tile = "east"; cb.Label.String = "log_{10} birds km^{-3}";
export_fixed_figure(fig, fullfile(figure_dir, "GD_A_maps.png"), Resolution=150)
close(fig)

%% 4. Recovery by range, four densities --------------------------------------------
test_density = [5 20 50 200];
n_replicate = 3;
method_name = ["output, ring mean" "our cleaning" "window mean" "estimator"];
recovered = nan(numel(bin_km), numel(method_name), n_replicate, numel(test_density));
gate_rmse = nan(numel(bin_km), 2, numel(test_density));          % window mean, estimator
for i_d = 1:numel(test_density)
    squared = zeros(numel(bin_km), 2);
    for i_rep = 1:n_replicate
        rng(100 * i_d + i_rep)
        eta = gd_simulate_sweep(test_density(i_d), range_m, az_step_deg, beam);
        dbz = gd_observe(eta, detection_dbz, background_dbz);
        output = density_from_dbz(dbz);
        output(isnan(dbz)) = 0;
        window_mean = gd_window_mean(output, range_m, az_step_deg, radius_m);
        [lo, hi] = gd_gate_bounds(dbz, zeros(size(dbz)), range_m, beam, detection_dbz, s_floor_cm2);
        estimate = gd_estimate_sweep(lo, hi, range_m, az_step_deg, beam, table_s, radius_m=radius_m);
        values = {output, gd_read_rd_filters(dbz, range_m, az_step_deg), window_mean, estimate.density};
        for i_m = 1:numel(method_name)
            recovered(:, i_m, i_rep, i_d) = bin_mean(values{i_m} / test_density(i_d));
        end
        squared = squared + [bin_rmse(window_mean / test_density(i_d)) bin_rmse(estimate.density / test_density(i_d))] .^ 2;
    end
    gate_rmse(:, :, i_d) = sqrt(squared / n_replicate);
end
summary = table(bin_km);
for i_d = 1:numel(test_density)
    for i_m = 1:numel(method_name)
        summary.(matlab.lang.makeValidName(sprintf("d%d_%s", test_density(i_d), method_name(i_m)))) = ...
            round(mean(recovered(:, i_m, :, i_d), 3), 2);
    end
    summary.(sprintf("d%d_gate_rmse_window_mean", test_density(i_d))) = round(gate_rmse(:, 1, i_d), 2);
    summary.(sprintf("d%d_gate_rmse_estimator", test_density(i_d))) = round(gate_rmse(:, 2, i_d), 2);
end
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_A_recovery.csv"))

method_colors = [0.1 0.4 0.8; 0.8 0.2 0.15; 0.5 0.5 0.5; 0 0 0];
fig = create_export_figure(3600, 2100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, numel(test_density), TileSpacing="tight", Padding="tight");
for i_d = 1:numel(test_density)
    ax = nexttile(t, i_d); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_m = 1:numel(method_name)
        gd_plot_band(ax, bin_km, squeeze(recovered(:, i_m, :, i_d)), method_colors(i_m, :), method_name(i_m))
    end
    yline(ax, 1, "k--", HandleVisibility="off")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 60]); ylim(ax, [0.005 3])
    xticks(ax, [2 5 10 20 30 60]); yticks(ax, [0.01 0.03 0.1 0.3 1 3])
    xlabel(ax, "Range (km)"); ylabel(ax, "Recovered / true (ring mean)")
    title(ax, sprintf("%d birds km^{-3}", test_density(i_d)), FontWeight="normal")
    if i_d == 1
        legend(ax, Location="southeast", FontSize=10)
    end
    ax = nexttile(t, numel(test_density) + i_d); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, bin_km, gate_rmse(:, 1, i_d), "-o", Color=method_colors(3, :), LineWidth=2, DisplayName="window mean")
    plot(ax, bin_km, gate_rmse(:, 2, i_d), "-o", Color=method_colors(4, :), LineWidth=2, DisplayName="estimator")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [2 60]); ylim(ax, [0.02 5])
    xticks(ax, [2 5 10 20 30 60]); xlabel(ax, "Range (km)"); ylabel(ax, "Per-gate relative RMS error")
    if i_d == 1
        legend(ax, Location="northeast", FontSize=10)
    end
end
export_fixed_figure(fig, fullfile(figure_dir, "GD_A_recovery.png"), Resolution=150)
close(fig)
