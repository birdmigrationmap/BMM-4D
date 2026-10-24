%% Gate density C: calibrating the cross-section spread and the beam volume
% The two settings the estimator depends on most (GD_B) are the spread of
% bird cross-sections (log-sd) and the gate volume (beam width, gate
% length, the radar's own volume convention). Both change the shape of the
% gate distribution, so they can be estimated from the data by profile
% likelihood: for each pair on a grid, the density of every range ring of
% every sweep is estimated, and the maximum log likelihoods are summed.
% The gate floor is kept at the same reflectivity for every pair, so the
% likelihoods compare the same data.
% 1. Recovery test on simulated sweeps with known values.
% 2. The sandbox's real sweeps (10 radars, frame 3, cuts below 3.5 deg),
%    separately for rings at 15-45 km (little clutter filtering) and at
%    4-15 km (most gates filtered by the radar).
% Outputs: studies/gate_density/figures/GD_C_calibration.png, .csv, .mat
% (the chosen pair, used by GD_D_sandbox).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Grid of model settings ------------------------------------------------------
sigma_grid = [0.75 1 1.25 1.5 1.75 2 2.25 2.5];
volume_grid = 2 .^ (-0.5:0.25:2.5);
n_sigma = numel(sigma_grid);
n_volume = numel(volume_grid);
density_grid = logspace(-1, 3.7, 236);
tables = cell(n_sigma, 1);
for i_s = 1:n_sigma
    tables{i_s} = gd_bird_sum_table(gd_beam(sigma_log_sd=sigma_grid(i_s)), 300, 20000);
end
model_beam = @(i_s, i_v) gd_beam(sigma_log_sd=sigma_grid(i_s), width_deg=0.95 * sqrt(volume_grid(i_v)));
s_floor_cm2 = 0.1;                                  % at volume factor 1; scaled with the volume
ring_edges_near_m = (4:1:15) * 1000;
ring_edges_far_m = (15:2.5:45) * 1000;

%% 1. Recovery test on simulated sweeps -------------------------------------------
range_m = (2:0.25:60)' * 1000;
az_step_deg = 0.5;
detection_dbz = -28 + 20 * log10(range_m / 5000);
truths = [1 1; 1.25 1.41; 0.75 0.71];               % [log-sd, volume factor] used to simulate
synthetic_ll = nan(n_sigma, n_volume, size(truths, 1));
for i_truth = 1:size(truths, 1)
    data_beam = gd_beam(sigma_log_sd=truths(i_truth, 1), width_deg=0.95 * sqrt(truths(i_truth, 2)));
    sweeps = cell(2, 1);
    for i_rep = 1:2
        rng(70 + 10 * i_truth + i_rep)
        eta = gd_simulate_sweep(gd_random_field(range_m, az_step_deg, 20, 0.4, 15000), range_m, az_step_deg, data_beam);
        sweeps{i_rep} = gd_observe(eta, detection_dbz, -15);
    end
    for i_s = 1:n_sigma
        for i_v = 1:n_volume
            beam = model_beam(i_s, i_v);
            total = 0;
            for i_rep = 1:2
                [lo, hi] = gd_gate_bounds(sweeps{i_rep}, zeros(size(sweeps{i_rep})), range_m, beam, detection_dbz, ...
                    s_floor_cm2 * volume_grid(i_v));
                rings = gd_estimate_rings(lo, hi, range_m, beam, tables{i_s}, [ring_edges_near_m ring_edges_far_m(2:end)], density_grid);
                total = total + sum(rings.log_likelihood, "omitnan");
            end
            synthetic_ll(i_s, i_v, i_truth) = total;
        end
    end
    fprintf("[calibration] synthetic truth %g / %g done\n", truths(i_truth, :))
end

%% 2. Real sweeps of the sandbox --------------------------------------------------
c = case_settings("sandbox");
manifest = load(fullfile(c.data_dir, "manifest.mat"));
radar_ids = 1:3:28;
real_sweeps = {};
for i_r = radar_ids
    level2 = bmm4d.radar.read_level2_moments(manifest.raw_files(i_r, 3), ["REF" "CFP"], 3.5);
    level2 = level2(arrayfun(@(s) ~isempty(s.CFP) && ~isempty(s.REF), level2));
    for i_e = 1:numel(level2)
        n_gate = min(size(level2(i_e).REF, 1), size(level2(i_e).CFP, 1));
        dbz = level2(i_e).REF(1:n_gate, :);
        dbz(z_to_refl(idb(dbz)) / 11 > 3000) = nan;      % clutter and rain: below are counted as no echo
        detection = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        real_sweeps(end + 1, :) = {dbz, level2(i_e).REF_range(1:n_gate), detection}; %#ok<SAGROW>
    end
end
real_ll = nan(n_sigma, n_volume, 2);                % [.., near / far]
ring_sets = {ring_edges_near_m, ring_edges_far_m};
for i_s = 1:n_sigma
    for i_v = 1:n_volume
        beam = model_beam(i_s, i_v);
        totals = [0 0];
        for i_w = 1:size(real_sweeps, 1)
            [dbz, sweep_range, detection] = real_sweeps{i_w, :};
            [lo, hi] = gd_gate_bounds(dbz, zeros(size(dbz)), sweep_range, beam, detection, s_floor_cm2 * volume_grid(i_v));
            for i_set = 1:2
                rings = gd_estimate_rings(lo, hi, sweep_range, beam, tables{i_s}, ring_sets{i_set}, density_grid);
                totals(i_set) = totals(i_set) + sum(rings.log_likelihood, "omitnan");
            end
        end
        real_ll(i_s, i_v, :) = totals;
    end
    fprintf("[calibration] real data, log-sd %g done\n", sigma_grid(i_s))
end

%% Best pairs -------------------------------------------------------------------
best_of = @(ll) deal(sigma_grid(mod(find(ll == max(ll(:)), 1) - 1, n_sigma) + 1), ...
    volume_grid(ceil(find(ll == max(ll(:)), 1) / n_sigma)));
rows = {};
for i_truth = 1:size(truths, 1)
    [best_sigma, best_volume] = best_of(synthetic_ll(:, :, i_truth));
    rows(end + 1, :) = {sprintf("simulated, truth %g / %.2f", truths(i_truth, :)), best_sigma, best_volume}; %#ok<SAGROW>
end
set_names = ["real, rings 4-15 km" "real, rings 15-45 km"];
for i_set = 1:2
    [best_sigma, best_volume] = best_of(real_ll(:, :, i_set));
    rows(end + 1, :) = {set_names(i_set), best_sigma, best_volume}; %#ok<SAGROW>
end
calibration = cell2table(rows, VariableNames=["data" "best_sigma_log_sd" "best_volume_factor"]);
disp(calibration)
writetable(calibration, fullfile(figure_dir, "GD_C_calibration.csv"))
[calibrated_sigma_log_sd, calibrated_volume_factor] = best_of(real_ll(:, :, 2));
save(fullfile(figure_dir, "GD_C_calibration.mat"), "calibrated_sigma_log_sd", "calibrated_volume_factor", ...
    "sigma_grid", "volume_grid", "synthetic_ll", "real_ll", "truths")

%% Plot: log-likelihood surfaces --------------------------------------------------
% Colour: log10(1 + max - log likelihood), so 0 is the best pair
surfaces = {synthetic_ll(:, :, 1), synthetic_ll(:, :, 2), synthetic_ll(:, :, 3), real_ll(:, :, 1), real_ll(:, :, 2)};
surface_titles = [compose("Simulated %g / x%.2f", truths)' "Real, 4-15 km" "Real, 15-45 km"];
fig = bmm4d.plot.create_export_figure(3600, 1050, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(surfaces), TileSpacing="tight", Padding="tight");
title(t, "Profile log likelihood over cross-section log-sd and beam volume factor (white cross: best; black circle: simulated truth; white line: 3 below the best)", ...
    FontSize=13, Color="k")
for i_p = 1:numel(surfaces)
    ax = nexttile(t); hold(ax, "on"); box(ax, "on")
    drop = max(surfaces{i_p}(:)) - surfaces{i_p};
    imagesc(ax, log2(volume_grid), sigma_grid, log10(1 + drop))
    contour(ax, log2(volume_grid), sigma_grid, drop, [3 3], "w-", LineWidth=1.5)
    [i_best_s, i_best_v] = find(drop == 0, 1);
    plot(ax, log2(volume_grid(i_best_v)), sigma_grid(i_best_s), "wx", MarkerSize=14, LineWidth=3)
    if i_p <= size(truths, 1)
        plot(ax, log2(truths(i_p, 2)), truths(i_p, 1), "ko", MarkerSize=12, LineWidth=2)
    end
    axis(ax, "tight"); set(ax, YDir="normal", FontSize=11); colormap(ax, flipud(parula)); clim(ax, [0 6])
    xticks(ax, -0.5:1:2.5); xticklabels(ax, compose("%.2f", 2 .^ (-0.5:1:2.5))); xtickangle(ax, 0)
    xlabel(ax, "Beam volume factor"); ylabel(ax, "Cross-section log-sd")
    title(ax, surface_titles(i_p), FontWeight="normal", FontSize=12)
end
cb = colorbar(ax, Color="k"); cb.Layout.Tile = "east"; cb.Label.String = "log_{10}(1 + max - log likelihood)";
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_C_calibration.png"), Resolution=150)
close(fig)
