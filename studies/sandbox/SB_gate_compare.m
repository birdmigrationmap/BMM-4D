%% Sandbox: real gates near the radar against simulated birds
% Simulated sparse birds keep the sweeps consistent (SB_sparse_sim), the
% real data do not on quiet nights. What is different in the real gates?
% For each radar-night, the density profile is taken from the lower sweeps
% far from the radar (0.9-1.8 deg at 25-45 km: many birds per gate, mean of
% the radar output with no echo = 0, by 100 m of height); birds are then
% simulated with that profile through one upper sweep (nearest 2.4 deg,
% its own radials and gates; gd_simulate_sweep, gd_observe with the
% detection floor of GD_B). Real and simulated gates are compared in range
% bands near the radar:
%   - dBZ distribution (no echo at the left edge);
%   - mean density relative to the profile (raw mean, no echo = 0);
%   - the radar's clutter-filter flags of the real gates (CFP: 0 filter
%     not applied, 1 point clutter filter, >= 8 power removed).
% Two quiet radar-nights (worst scores) and two busy ones (best).
% Outputs: figures/gate_compare.png, gate_compare.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
cache_dir = bmm4d.radar.volume_cache_dir(c);
cases = table(["KGRB" "KABR" "KSGF" "KDMX"]', ["2023-05-20" "2025-09-27" "2025-10-07" "2023-05-16"]', ...
    VariableNames=["radar" "night"]);
band_edges_km = [3 6 10 15 25];
dbz_edges = -30:2:40;
height_edges = 0:100:3000;

%% Read, simulate, compare ------------------------------------------------------
summary = table();
panels = cell(height(cases), 1);
for i_c = 1:height(cases)
    files = dir(fullfile(cache_dir, cases.radar(i_c) + "*.mat"));
    time = datetime(extractBetween(string({files.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
    in_night = find(string(dateshift(time - hours(12), "start", "day"), "yyyy-MM-dd") == cases.night(i_c));
    [~, k] = min(abs(time(in_night) - (datetime(cases.night(i_c)) + days(1) + hours(4))));
    raw_file = load(fullfile(files(in_night(k)).folder, files(in_night(k)).name)).volume.source_file;
    sweeps = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], 7);
    sweeps = sweeps(arrayfun(@(s) ~isempty(s.REF) && ~isempty(s.CFP), sweeps));   % surveillance and batch passes
    % Profile from far, low sweeps
    heights = []; values = [];
    for sweep = sweeps(abs([sweeps.elev] - 1.3) < 0.6)
        far = sweep.REF_range >= 25000 & sweep.REF_range <= 45000;
        z = 10 .^ (sweep.REF(far, :) / 10); z(sweep.REF_code(far, :) == 0) = 0; z(sweep.REF_code(far, :) == 1) = NaN;
        [~, height_m] = slant2ground(sweep.REF_range(far), sweep.elev);
        heights = [heights; height_m]; values = [values; z_to_refl(mean(z, 2, "omitnan")) / 11]; %#ok<AGROW>
    end
    bin = discretize(heights, height_edges); ok = ~isnan(bin);
    profile = accumarray(bin(ok), values(ok), [numel(height_edges) - 1 1], @mean, NaN);
    profile = fillmissing(profile, "nearest");
    centre = (height_edges(1:end - 1) + height_edges(2:end))' / 2;
    % Upper sweep nearest 2.4 deg: real gates and simulated birds
    [~, i_s] = min(abs([sweeps.elev] - 2.4));
    sweep = sweeps(i_s);
    range_m = sweep.REF_range(sweep.REF_range <= 26000);
    n_gate = numel(range_m);
    real_dbz = sweep.REF(1:n_gate, :); real_dbz(sweep.REF_code(1:n_gate, :) == 1) = NaN;   % range folded: drop
    real_code = sweep.REF_code(1:n_gate, :);
    cfp_code = sweep.CFP_code(1:n_gate, :);
    [~, height_m] = slant2ground(range_m, sweep.elev);
    truth = interp1(centre, profile, height_m, "linear", "extrap");
    az_step = 360 / numel(sweep.az);
    rng(i_c)
    eta = gd_simulate_sweep(truth, range_m, az_step, bmm4d.radar.gate.beam());
    sim_dbz = gd_observe(eta, -28 + 20 * log10(range_m / 5000), -Inf);
    beam = bmm4d.radar.gate.beam(gate_length_m=median(diff(range_m)));
    panel = struct("elev", sweep.elev, "profile", profile, "bands", []);
    for i_b = 1:numel(band_edges_km) - 1
        rows = range_m >= band_edges_km(i_b) * 1000 & range_m < band_edges_km(i_b + 1) * 1000;
        real_band = real_dbz(rows, :); real_band = real_band(~isnan(real_band) | real_code(rows, :) == 0);
        real_no_echo = mean(real_code(rows, :) == 0, "all");
        sim_band = sim_dbz(rows, :);
        sim_no_echo = mean(isnan(sim_band), "all");
        real_mean = mean(z_to_refl(idb(fillmissing(real_dbz(rows, :), "constant", -Inf))) / 11, "all", "omitnan");
        sim_mean = mean(z_to_refl(idb(fillmissing(sim_band, "constant", -Inf))) / 11, "all");
        truth_mean = mean(truth(rows));
        codes = double(cfp_code(rows, :));
        row = table(cases.radar(i_c), cases.night(i_c), round(sweep.elev, 1), string(band_edges_km(i_b)) + "-" + ...
            string(band_edges_km(i_b + 1)) + " km", round(truth_mean, 1), round(mean(beam.volume_sample_km3(range_m(rows))) * truth_mean, 2), ...
            round(real_mean / truth_mean, 2), round(sim_mean / truth_mean, 2), round(real_no_echo, 2), round(sim_no_echo, 2), ...
            round(mean(real_band > 15), 3), round(mean(sim_band(~isnan(sim_band)) > 15) * (1 - sim_no_echo), 3), ...
            round(mean(codes == 0, "all"), 2), round(mean(codes == 1, "all"), 2), round(mean(codes >= 11, "all"), 2), ...
            VariableNames=["radar" "night" "elev" "band" "profile_density" "birds_per_gate" "real_mean_over_profile" ...
            "sim_mean_over_profile" "real_no_echo" "sim_no_echo" "real_above_15dBZ" "sim_above_15dBZ" ...
            "cfp_not_applied" "cfp_point" "cfp_removed_3dB"]);
        summary = [summary; row]; %#ok<AGROW>
        panel.bands(i_b).real = histcounts(min(max(real_dbz(rows, :), dbz_edges(1)), dbz_edges(end) - 0.1), dbz_edges) / numel(real_dbz(rows, :));
        panel.bands(i_b).sim = histcounts(min(max(sim_band, dbz_edges(1)), dbz_edges(end) - 0.1), dbz_edges) / numel(sim_band);
        panel.bands(i_b).real_no_echo = real_no_echo;
        panel.bands(i_b).sim_no_echo = sim_no_echo;
    end
    panels{i_c} = panel;
end
disp(summary)
writetable(summary, fullfile(figure_dir, "gate_compare.csv"))

%% Plot results -----------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 750 * height(cases), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, height(cases), numel(band_edges_km) - 1, TileSpacing="compact", Padding="tight");
title(t, "Gate values of an upper sweep near the radar: real (orange) against simulated birds at the density measured far away (blue)", ...
    FontSize=16, Color="k")
centre_dbz = (dbz_edges(1:end - 1) + dbz_edges(2:end)) / 2;
for i_c = 1:height(cases)
    for i_b = 1:numel(band_edges_km) - 1
        ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
        b = panels{i_c}.bands(i_b);
        s = summary(summary.radar == cases.radar(i_c) & summary.band == string(band_edges_km(i_b)) + "-" + string(band_edges_km(i_b + 1)) + " km", :);
        stairs(ax, dbz_edges(1:end - 1), b.sim * 100, Color=[0 0.45 0.74], LineWidth=2)
        stairs(ax, dbz_edges(1:end - 1), b.real * 100, Color=[0.85 0.33 0.1], LineWidth=2)
        set(ax, FontSize=11); xlim(ax, [-30 40])
        text(ax, 0.98, 0.95, sprintf("no echo: real %.0f%%, sim %.0f%%\nmean / profile: real %.2f, sim %.2f\n%.2f birds per gate\nCFP point %.0f%%, removed ≥3 dB %.0f%%", ...
            100 * b.real_no_echo, 100 * b.sim_no_echo, s.real_mean_over_profile, s.sim_mean_over_profile, s.birds_per_gate, ...
            100 * s.cfp_point, 100 * s.cfp_removed_3dB), Units="normalized", HorizontalAlignment="right", VerticalAlignment="top", ...
            FontSize=10, BackgroundColor="w")
        if i_c == 1, title(ax, string(band_edges_km(i_b)) + "-" + string(band_edges_km(i_b + 1)) + " km", FontSize=13), end
        if i_b == 1
            ylabel(ax, sprintf("%s %s, %.1f°\n%% of gates per 2 dB", cases.radar(i_c), cases.night(i_c), panels{i_c}.elev), ...
                FontWeight="bold", FontSize=11)
        end
        if i_c == height(cases), xlabel(ax, "Gate reflectivity (dBZ)"), end
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "gate_compare.png"), Resolution=150)
close(fig)
