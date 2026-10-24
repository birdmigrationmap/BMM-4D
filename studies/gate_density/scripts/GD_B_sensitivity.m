%% Gate density B: sensitivity of the estimator to its settings and assumptions
% Simulated sweeps (base: 20 birds km^-3 with a smooth lognormal pattern,
% log-sd 0.4 over 15 km; 2 km window), one factor changed at a time:
% data side (what the world is like)
% - bird cross-section spread (log-sd; the model assumes 1)
% - beam volume (true volume / the model's: beam width, gate length)
% - background echo level (insects, noise)
% - small-scale heterogeneity inside the window (pattern log-sd, 5 km)
% - density
% model side (settings of the estimator)
% - detection limit error (dB)
% - floor below which a gate is "no bird" (cm^2)
% - window radius
% - largest number of birds per gate in the table (N_max)
% Scores per range band, against the true density averaged over the same
% window: bias (median ratio estimate / truth) and per-gate relative RMS
% error, for the estimator and, as a reference, the window mean.
% Outputs: studies/gate_density/figures/GD_B_sensitivity.png and .csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Base settings ------------------------------------------------------------
range_m = (2:0.25:60)' * 1000;
az_step_deg = 0.5;
detection_dbz = -28 + 20 * log10(range_m / 5000);
band_edges_km = [2 5 10 20 40 60];
band_label = compose("%d-%d km", band_edges_km(1:end - 1)', band_edges_km(2:end)');
band_of = discretize(range_m / 1000, band_edges_km);
n_band = numel(band_label);
n_replicate = 2;
base = struct("sigma_log_sd", 1, "volume_factor", 1, "background_dbz", -15, "pattern_log_sd", 0.4, ...
    "density", 20, "detection_error_db", 0, "s_floor_cm2", 0.1, "radius_m", 2000, "n_max", 300);
factors = {
    "sigma_log_sd",       [0.5 0.75 1 1.25 1.5],        "cross-section spread (log-sd, model 1)"
    "volume_factor",      [0.5 0.71 1 1.41 2],          "true / model beam volume"
    "background_dbz",     [-Inf -25 -15 -5 5],          "background echo (dBZ)"
    "pattern_log_sd",     [0 0.4 0.8 1.2],              "small-scale pattern (log-sd, 5 km)"
    "density",            [3 10 20 100 300],            "density (birds km^{-3})"
    "detection_error_db", [-6 -3 0 3 6],                "detection limit error (dB)"
    "s_floor_cm2",        [0.01 0.03 0.1 0.3 1],        "no-bird floor (cm^2)"
    "radius_m",           [1000 2000 3000 5000],        "window radius (m)"
    "n_max",              [100 200 300],                "largest number of birds per gate in the table"};
n_factor = size(factors, 1);
tables = containers.Map();                           % bird-sum tables by (log-sd, N_max)
get_table = @(beam, n_max) gd_bird_sum_table(beam, n_max, 20000);

%% Run --------------------------------------------------------------------------
rows = {};
for i_f = 1:n_factor
    [name, levels] = factors{i_f, 1:2};
    for level = levels
        setting = base;
        setting.(name) = level;
        model_beam = gd_beam();
        key = sprintf("%g", setting.n_max);
        if ~isKey(tables, key)
            tables(key) = get_table(model_beam, setting.n_max);
        end
        table_s = tables(key);
        % Data: true beam width scaled so its volume is volume_factor times the model's
        data_beam = gd_beam(sigma_log_sd=setting.sigma_log_sd, width_deg=0.95 * sqrt(setting.volume_factor));
        [ratio_estimate, ratio_mean] = deal(cell(n_band, 1));
        for i_rep = 1:n_replicate
            rng(1000 * i_f + 10 * find(levels == level) + i_rep)
            if strcmp(name, "pattern_log_sd")
                truth = gd_random_field(range_m, az_step_deg, setting.density, setting.pattern_log_sd, 5000);
            else
                truth = gd_random_field(range_m, az_step_deg, setting.density, setting.pattern_log_sd, 15000);
            end
            eta = gd_simulate_sweep(truth, range_m, az_step_deg, data_beam);
            dbz = gd_observe(eta, detection_dbz, setting.background_dbz);
            output = z_to_refl(idb(dbz)) / 11;
            output(isnan(dbz)) = 0;
            [lo, hi] = gd_gate_bounds(dbz, zeros(size(dbz)), range_m, model_beam, detection_dbz + setting.detection_error_db, ...
                setting.s_floor_cm2);
            estimate = gd_estimate_sweep(lo, hi, range_m, az_step_deg, model_beam, table_s, radius_m=setting.radius_m);
            target = gd_window_mean(truth, range_m, az_step_deg, setting.radius_m);
            window_mean = gd_window_mean(output, range_m, az_step_deg, setting.radius_m);
            for b = 1:n_band
                in_band = band_of == b;
                ratio_estimate{b} = [ratio_estimate{b}; reshape(estimate.density(in_band, :) ./ target(in_band, :), [], 1)];
                ratio_mean{b} = [ratio_mean{b}; reshape(window_mean(in_band, :) ./ target(in_band, :), [], 1)];
            end
        end
        for b = 1:n_band
            rows(end + 1, :) = {string(name), level, band_label(b), median(ratio_estimate{b}, "omitnan"), ...
                sqrt(mean((ratio_estimate{b} - 1) .^ 2, "omitnan")), median(ratio_mean{b}, "omitnan"), ...
                sqrt(mean((ratio_mean{b} - 1) .^ 2, "omitnan"))}; %#ok<SAGROW>
        end
        fprintf("[sensitivity] %s = %g done\n", name, level)
    end
end
sensitivity = cell2table(rows, VariableNames=["factor" "level" "band" "estimator_bias" "estimator_rmse" ...
    "window_mean_bias" "window_mean_rmse"]);
writetable(sensitivity, fullfile(figure_dir, "GD_B_sensitivity.csv"))

%% Plot: bias and error by range band, one panel per factor ---------------------------------
fig = create_export_figure(3600, 3600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, 3, TileSpacing="tight", Padding="tight");
title(t, "Estimator sensitivity on simulated sweeps (solid: estimator, dashed: window mean; legend: per-gate RMS error of the estimator over the bands)", ...
    FontSize=14, Color="k")
band_center = 1:n_band;                              % bands as categories
for i_f = 1:n_factor
    [name, levels, label] = factors{i_f, :};
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    level_colors = interp1([0 1], [0.85 0.55 0.2; 0.1 0.25 0.6], linspace(0, 1, numel(levels)));
    for i_l = 1:numel(levels)
        rows_l = sensitivity.factor == name & sensitivity.level == levels(i_l);
        is_base = levels(i_l) == base.(name);
        plot(ax, band_center, sensitivity.estimator_bias(rows_l), "-o", Color=level_colors(i_l, :), LineWidth=2 + is_base, ...
            DisplayName=sprintf("%g (RMS %.2f-%.2f)", levels(i_l), min(sensitivity.estimator_rmse(rows_l)), max(sensitivity.estimator_rmse(rows_l))))
        plot(ax, band_center, sensitivity.window_mean_bias(rows_l), "--", Color=level_colors(i_l, :), LineWidth=1, HandleVisibility="off")
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, FontSize=12); xlim(ax, [0.7 n_band + 0.3]); ylim(ax, [0.4 1.8]); xticks(ax, band_center)
    xticklabels(ax, band_label); xtickangle(ax, 0)
    ylabel(ax, "Estimate / truth (median)"); legend(ax, Location="northwest", FontSize=9, Interpreter="none")
    title(ax, label, FontWeight="normal")
end
export_fixed_figure(fig, fullfile(figure_dir, "GD_B_sensitivity.png"), Resolution=150)
close(fig)
