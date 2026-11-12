%% Sandbox: do simulated sparse birds make the sweeps disagree?
% On quiet nights the upper sweeps read lower than the lower ones at the
% same height (SB_obs_sections score, down to -0.8 below 3 birds km^-3
% aloft), and neither the beam geometry nor refraction produces that for
% any vertical profile (README). Here the observation chain alone: birds
% placed at random (Poisson) in every sweep of a VCP 212 volume with the
% density of a vertical profile at each gate's beam centre (so the true
% field is consistent between sweeps by construction), lognormal cross-
% sections and the two-way beam weighting of the gate model
% (studies/gate_density gd_simulate_sweep; 0.5 deg radials up to 1.3 deg,
% 1 deg above), 0.5 dB coding and a detection floor (gd_observe), then the
% production gate estimator (bmm4d.radar.filter.estimate_density). Scored
% like the real data: median over azimuth per sweep and 0.5 km, then the
% median log2(upper / lower sweep) at the same height, upper within 20 km.
% Also the raw mean (no echo = 0) and the truth.
% Profiles (birds km^-3, h above the antenna): quiet to busy, with and
% without a dense layer near the ground.
% Outputs: figures/sparse_sim.png, sparse_sim.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
elevations = [0.5 0.9 1.3 1.8 2.4 3.1 4.0 5.1 6.4];
range_m = (2125:250:45000)';
detection_dbz = -28 + 20 * log10(range_m / 5000);       % as in GD_B
parms = struct("smoothing", 2000, "minRange", 3000, "maxRange", 45000, "maxAlt", 3500, "noEchoRange", 150000);
range_edges_m = 3000:500:45000;
scenarios = table( ...
    ["uniform 2" "uniform 10" "uniform 50" "uniform 200" "quiet, low layer" "medium, low layer" "busy, low layer"]', ...
    {@(h) 2 + 0 * h; @(h) 10 + 0 * h; @(h) 50 + 0 * h; @(h) 200 + 0 * h; ...
    @(h) 2 + 30 * exp(-h / 250); @(h) 15 + 60 * exp(-h / 400); @(h) 50 + 200 * exp(-h / 800)}, ...
    VariableNames=["name" "profile"]);

%% Simulate, estimate and score -------------------------------------------------
results = table();
sections = cell(height(scenarios), 1);
for i_p = 1:height(scenarios)
    profile = scenarios.profile{i_p};
    sweep_sections = cell(numel(elevations), 1);
    parfor i_e = 1:numel(elevations)
        rng(100 * i_p + i_e)
        e = elevations(i_e);
        az_step = 0.5 + 0.5 * (e > 1.5);
        [~, height_m] = slant2ground(range_m, e);
        truth = profile(height_m);
        beam = bmm4d.radar.gate.beam();
        eta = gd_simulate_sweep(truth, range_m, az_step, beam);
        dbz = gd_observe(eta, detection_dbz, -Inf);
        gate = struct("elev", e, "range", range_m, "az", (0:az_step:360 - az_step)', "data", dbz, "unit", "dBZ", "no_echo", isnan(dbz));
        estimated = bmm4d.radar.filter.estimate_density(struct("gates", gate), parms);
        density = estimated.gates.data;
        raw = z_to_refl(idb(dbz)) / 11;
        raw(isnan(dbz)) = 0;
        range_bin = discretize(range_m, range_edges_m);
        ok = ~isnan(range_bin);
        n_bin = numel(range_edges_m) - 1;
        centre = (range_edges_m(1:end - 1) + range_edges_m(2:end))' / 2;
        [~, centre_height] = slant2ground(centre, e);
        sweep_sections{i_e} = table(repmat(e, n_bin, 1), centre / 1000, centre_height, ...
            accumarray(range_bin(ok), median(density(ok, :), 2, "omitnan"), [n_bin 1], @median, NaN), ...
            accumarray(range_bin(ok), mean(raw(ok, :), 2), [n_bin 1], @mean, NaN), ...
            profile(centre_height), ...
            beam.volume_sample_km3(centre) .* profile(centre_height), ...
            VariableNames=["elev" "range_km" "height_m" "estimator" "raw" "truth" "birds_per_gate"]);
    end
    section = vertcat(sweep_sections{:});
    sections{i_p} = section;
    row = table(scenarios.name(i_p), VariableNames="scenario");
    for quantity = ["truth" "raw" "estimator"]
        row.(quantity) = round(pair_score(section, quantity, elevations), 2);
    end
    aloft = section.height_m >= 500 & section.height_m < 1500;
    row.density_aloft = round(median(section.truth(aloft)), 1);
    row.birds_per_gate_upper_near = round(median(section.birds_per_gate(section.elev >= 1.8 & section.range_km <= 20)), 2);
    results = [results; row]; %#ok<AGROW>
    disp(row)
end
writetable(results, fullfile(figure_dir, "sparse_sim.csv"))

%% Plot results -----------------------------------------------------------
% Left: score against density aloft, simulated (estimator, raw) over the
% real radar-nights. Right: estimator / truth by expected birds per gate.
observed = readtable(fullfile(figure_dir, "obs_sections_score.csv"), ReadRowNames=true, ReadVariableNames=true, ...
    VariableNamingRule="preserve");
fig = bmm4d.plot.create_export_figure(3200, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Simulated birds through a VCP 212 volume and the gate estimator, scored like the real data", FontSize=16, Color="k")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
if isfile(fullfile(figure_dir, "obs_sections_density.csv"))
    real = readtable(fullfile(figure_dir, "obs_sections_density.csv"));
    scatter(ax, real.density_aloft, real.score, 20, [0.7 0.7 0.7], "filled", DisplayName="real radar-nights")
end
scatter(ax, results.density_aloft, results.estimator, 120, [0 0.45 0.74], "filled", DisplayName="simulated, estimator")
scatter(ax, results.density_aloft, results.raw, 120, [0.85 0.33 0.1], "d", "filled", DisplayName="simulated, raw mean")
scatter(ax, results.density_aloft, results.truth, 60, "k", "x", LineWidth=2, DisplayName="simulated, truth")
yline(ax, 0, "k:", HandleVisibility="off")
set(ax, XScale="log", FontSize=13); xlim(ax, [0.5 500]); ylim(ax, [-2 0.6])
xlabel(ax, "Density at 500-1500 m (birds km^{-3})"); ylabel(ax, "Score: log_2 upper / lower sweep")
legend(ax, Location="southeast", FontSize=12)
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
all_sections = vertcat(sections{:});
lambda_edges = 10 .^ (-2:0.25:3);
lambda_bin = discretize(all_sections.birds_per_gate, lambda_edges);
ok = ~isnan(lambda_bin) & all_sections.range_km >= 3;
lambda_centre = sqrt(lambda_edges(1:end - 1) .* lambda_edges(2:end));
for quantity = ["estimator" "raw"]
    ratio = accumarray(lambda_bin(ok), all_sections.(quantity)(ok) ./ all_sections.truth(ok), [numel(lambda_centre) 1], @median, NaN);
    count = accumarray(lambda_bin(ok), 1, [numel(lambda_centre) 1]);
    ratio(count < 5) = NaN;
    plot(ax, lambda_centre, ratio, "-o", LineWidth=2, DisplayName=quantity + " / truth (median over azimuth)")
end
yline(ax, 1, "k:", HandleVisibility="off")
set(ax, XScale="log", FontSize=13); ylim(ax, [0 2])
xlabel(ax, "Expected birds per gate \lambda"); ylabel(ax, "Value / truth")
legend(ax, Location="southeast", FontSize=12)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sparse_sim.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function score = pair_score(section, quantity, elevations)
% Median log2(upper / lower sweep) at the same height, upper sweep within
% 20 km, both above 0.5 birds km^-3 (as SB_obs_sections)
log_ratio = [];
for i_e = 2:numel(elevations)
    upper = section(section.elev == elevations(i_e) & section.range_km <= 20 & section.(quantity) > 0.5, :);
    lower = section(section.elev == elevations(i_e - 1) & section.(quantity) > 0.5, :);
    if height(upper) == 0 || height(lower) < 2, continue, end
    log_ratio = [log_ratio; log2(upper.(quantity)) - interp1(lower.height_m, log2(lower.(quantity)), upper.height_m)]; %#ok<AGROW>
end
score = median(log_ratio, "omitnan");
end
