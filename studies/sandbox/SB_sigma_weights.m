%% Sandbox: observation errors from the expected birds per gate
% The inversion is least squares on linear density, so an observation
% error enters as a weight 1 / sigma^2 per row and time
% (solve_advection_lean ROW_WEIGHTS [n_rows n_time]). On quiet nights the
% sweeps disagree near the radar by far more than the gate estimator's own
% sd (README, "Does the estimator's own uncertainty cover ..."): the error
% is a model error that grows as the gates hold fewer birds. Here:
% 1. Sweep pairs: the gate estimator on every sweep of one volume per
%    radar-night (the ten radar-nights of SB_disagreement_cases), median
%    over azimuth per 1 km range bin; for upper-sweep bins within 20 km the
%    offset ln(upper / lower sweep at the same height), against the
%    expected birds per gate lambda = lower-sweep density x counting
%    volume (250 m gates).
% 2. Error model sigma(lambda) = sqrt(sigma_0^2 + a^2 / lambda), the RMS
%    offset per lambda bin (sd in log density; no threshold, rows are never
%    removed).
% 3. Solves with weights from that model, lambda from the radar prior
%    (sandbox_run SIGMA_MODEL): only the lambda dependence (absolute
%    errors as now, down-weighted below a few birds per gate) and relative
%    errors (sigma x prediction). Scored by radar cross-validation on the
%    sandbox case (7 folds of 4 radars, as the removed two-level study)
%    and by the radar-distance bias of the full solve on both nights.
% Outputs: figures/sigma_pairs.csv, sigma_weights.png, sigma_weights.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
cache_dir = bmm4d.radar.volume_cache_dir(c);
cases = table(["KGRB" "KABR" "KVNX" "KHPX" "KIWX" "KGRR" "KDDC" "KLOT" "KSGF" "KDMX"]', ...
    ["2023-05-20" "2025-09-27" "2021-10-21" "2023-05-20" "2023-05-16" "2021-05-11" "2021-10-21" "2026-10-21" "2025-10-07" "2023-05-16"]', ...
    VariableNames=["radar" "night"]);

%% Sweep pairs ------------------------------------------------------------
table_s = bmm4d.radar.gate.bird_sum_table(bmm4d.radar.gate.beam(), 2000, 20000);
pairs = table();
for i_c = 1:height(cases)
    % Volume nearest 04 UTC of the night
    files = dir(fullfile(cache_dir, cases.radar(i_c) + "*.mat"));
    time = datetime(extractBetween(string({files.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
    in_night = find(string(dateshift(time - hours(12), "start", "day"), "yyyy-MM-dd") == cases.night(i_c));
    [~, k] = min(abs(time(in_night) - (datetime(cases.night(i_c)) + days(1) + hours(4))));
    raw_file = load(fullfile(files(in_night(k)).folder, files(in_night(k)).name)).volume.source_file;
    sweeps = bmm4d.radar.read_level2_moments(raw_file, "REF", 7);
    sweeps = sweeps(arrayfun(@(s) ~isempty(s.REF), sweeps));
    [~, keep] = unique(round([sweeps.elev], 1));
    sweeps = sweeps(keep);
    % Gate estimator as in estimate_density, median per 1 km range bin
    section = table();
    for sweep = sweeps
        range_m = sweep.REF_range(:);
        rows = range_m <= 40000;
        range_m = range_m(rows);
        dbz = sweep.REF(rows, :);
        dbz(sweep.REF_code(rows, :) == 1) = NaN;          % range folded: no information
        no_echo = sweep.REF_code(rows, :) == 0;
        beam = bmm4d.radar.gate.beam(gate_length_m=median(diff(range_m)));
        detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, 0.1);
        estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, 360 / size(dbz, 2), beam, table_s, radius_m=1000);
        [~, height_m] = slant2ground(range_m, sweep.elev);
        range_bin = discretize(range_m, 3000:1000:40000);
        for b = unique(range_bin(~isnan(range_bin)))'
            in_bin = range_bin == b;
            section = [section; table(round(sweep.elev, 1), mean(range_m(in_bin)) / 1000, mean(height_m(in_bin)), ...
                median(estimate.density(in_bin, :), "all", "omitnan"), median(estimate.log_sd(in_bin, :), "all", "omitnan"), ...
                VariableNames=["elev" "range_km" "height_m" "density" "log_sd"])]; %#ok<AGROW>
        end
    end
    % Upper-sweep bins within 20 km against the next lower sweep at the same height
    elevations = unique(section.elev);
    section.lower_density = nan(height(section), 1);
    for i_e = 2:numel(elevations)
        upper = find(section.elev == elevations(i_e) & section.range_km <= 20);
        lower = section(section.elev == elevations(i_e - 1), :);
        section.lower_density(upper) = exp(interp1(lower.height_m, log(lower.density), section.height_m(upper)));
    end
    section = section(~isnan(section.lower_density), :);
    section.radar = repmat(cases.radar(i_c), height(section), 1);
    pairs = [pairs; section]; %#ok<AGROW>
end
pairs.offset = log(pairs.density ./ pairs.lower_density);
beam = bmm4d.radar.gate.beam();
pairs.lambda = pairs.lower_density .* beam.volume_sample_km3(pairs.range_km * 1000);
writetable(pairs, fullfile(figure_dir, "sigma_pairs.csv"))

%% Error model ------------------------------------------------------------
% RMS offset per half-decade of lambda. sigma_0 and a are set by hand to
% the bins from 0.1 to 10 birds per gate (a least-squares fit is pulled to
% sigma_0 0.57 by the few, noisy busy bins, where the offset also holds
% real differences between the two sweeps' azimuths).
lambda_edges = 10 .^ (-2:0.5:3);
lambda_bin = discretize(pairs.lambda, lambda_edges);
ok = ~isnan(lambda_bin) & isfinite(pairs.offset);
rms_offset = accumarray(lambda_bin(ok), pairs.offset(ok), [numel(lambda_edges) - 1 1], @(x) sqrt(mean(x .^ 2)), NaN);
median_log_sd = accumarray(lambda_bin(ok), pairs.log_sd(ok), [numel(lambda_edges) - 1 1], @median, NaN);
n_pair = accumarray(lambda_bin(ok), 1, [numel(lambda_edges) - 1 1]);
lambda_centre = sqrt(lambda_edges(1:end - 1) .* lambda_edges(2:end))';
sigma_0 = 0.4;
sigma_a = 0.5;
sigma_model = @(lambda) sqrt(sigma_0 ^ 2 + sigma_a ^ 2 ./ max(lambda, 1e-6));
disp(table(lambda_centre, n_pair, round(rms_offset, 2), round(median_log_sd, 2), round(sigma_model(lambda_centre), 2), ...
    VariableNames=["lambda" "n" "rms_offset" "estimator_log_sd" "sigma_model"]))

%% Solves -------------------------------------------------------------------
methods = {
    struct("name", "unweighted")
    struct("name", "sigma(lambda)", "sigma_model", sigma_model)
    struct("name", "sigma(lambda), relative", "sigma_model", sigma_model, "sigma_relative", true)};
results = table();
for case_name = ["sandbox" "night_20210512"]
    S = sandbox_load(case_name);
    for i_m = 1:numel(methods)
        option = methods{i_m};
        option.name = option.name + " " + case_name;
        option.video = false;
        result = sandbox_run(S, option);
        results = [results; table(case_name, string(methods{i_m}.name), result.bias_200km, result.misfit, ...
            result.misfit_trusted, result.detail, result.seconds, VariableNames=["case" "method" "bias_200km" ...
            "misfit" "misfit_trusted" "detail" "seconds"])]; %#ok<AGROW>
        ratios.(matlab.lang.makeValidName(case_name + "_" + i_m)) = struct("ratio", mean(result.ratio, 2, "omitnan"), ...
            "edges_km", result.edges_km);
    end
end
disp(results)

%% Radar cross-validation (sandbox case) ----------------------------------------------
% 7 folds of 4 radars, one per longitude band (as the removed two-level study
% SB_validation). Scores on the withheld rows by range band: log bias
% mean log((obs + 1) / (pred + 1)) (> 0: predicts too little) and its RMS.
S = sandbox_load("sandbox");
n_radar = size(S.radar_scans, 1);
radar_lon = cellfun(@(s) s.lon, S.radar_scans(:, 1));
radar_lat = cellfun(@(s) s.lat, S.radar_scans(:, 1));
n_fold = 7;
[~, by_lon] = sort(radar_lon);
band = zeros(n_radar, 1);
band(by_lon) = ceil((1:n_radar)' / (n_radar / 4));
fold = zeros(n_radar, 1);
for i_b = 1:4
    in_band = find(band == i_b);
    [~, by_lat] = sort(radar_lat(in_band));
    fold(in_band(by_lat)) = mod(0:numel(in_band) - 1, n_fold)' + 1;
end
observed = full(S.observations.values);
F_active = S.observations.F(:, S.active_id);
range_km = double(S.observations.row_range_m) / 1000;
band_edges_km = [3 15 30 60 120 250];
validation = table();
for i_m = 1:numel(methods)
    predicted = zeros(size(observed));
    for i_f = 1:n_fold
        S_fold = sandbox_fold(S, find(fold == i_f));
        option = methods{i_m};
        option.name = "held out " + option.name + " fold " + i_f;
        option.video = false;
        result = sandbox_run(S_fold, option);
        withheld = S_fold.withheld_rows;
        predicted(withheld, :) = F_active(withheld, :) * double(result.field);
    end
    for i_b = 1:numel(band_edges_km) - 1
        use = observed ~= 0 & range_km >= band_edges_km(i_b) & range_km < band_edges_km(i_b + 1);
        log_ratio = log((observed(use) + 1) ./ (max(predicted(use), 0) + 1));
        validation = [validation; table(string(methods{i_m}.name), string(band_edges_km(i_b)) + "-" + ...
            string(band_edges_km(i_b + 1)) + " km", mean(log_ratio), sqrt(mean(log_ratio .^ 2)), ...
            norm(predicted(use) - observed(use)) / norm(observed(use)), ...
            VariableNames=["method" "band" "log_bias" "log_rms" "relative_l2"])]; %#ok<AGROW>
    end
end
disp(validation)
writetable(results, fullfile(figure_dir, "sigma_weights.csv"))
writetable(validation, fullfile(figure_dir, "sigma_weights_validation.csv"))

%% Plot results -----------------------------------------------------------
% Left: the error model over the sweep pairs. Middle: column density
% against distance to the nearest radar (relative to its 15-200 km mean),
% per method, both nights. Right: cross-validation log bias by range band.
fig = bmm4d.plot.create_export_figure(3600, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
scatter(ax, pairs.lambda, abs(pairs.offset), 12, [0.7 0.7 0.7], "filled", DisplayName="|offset|, one range bin")
plot(ax, lambda_centre, rms_offset, "o", MarkerSize=9, LineWidth=2, Color=[0.85 0.33 0.1], DisplayName="RMS offset")
plot(ax, lambda_centre, median_log_sd, "s", MarkerSize=9, LineWidth=2, Color=[0 0.45 0.74], DisplayName="estimator log sd")
lambda_line = logspace(-2, 3, 200);
plot(ax, lambda_line, sigma_model(lambda_line), "k-", LineWidth=2, DisplayName=sprintf("\\sigma = (%.1f^2 + %.1f^2/\\lambda)^{1/2}", sigma_0, sigma_a))
set(ax, XScale="log", FontSize=13); xlim(ax, [0.01 300]); ylim(ax, [0 3])
xlabel(ax, "Expected birds per gate \lambda"); ylabel(ax, "Error in log density")
legend(ax, Location="northeast", FontSize=11)
title(ax, "Same-height sweep pairs, 10 radar-nights", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
colors = [0.5 0.5 0.5; 0 0.45 0.74; 0.85 0.33 0.1];
styles = ["-" "--"];
case_list = ["sandbox" "night_20210512"];
for i_c = 1:2
    for i_m = 1:numel(methods)
        r = ratios.(matlab.lang.makeValidName(case_list(i_c) + "_" + i_m));
        centre = (r.edges_km(1:end - 1) + r.edges_km(2:end)) / 2;
        reference = mean(r.ratio(centre >= 15 & centre <= 200), "omitnan");
        plot(ax, centre, r.ratio / reference, styles(i_c), Color=colors(i_m, :), LineWidth=2, ...
            DisplayName=methods{i_m}.name + ", " + case_list(i_c))
    end
end
yline(ax, 1, "k:", HandleVisibility="off")
xlim(ax, [0 200]); ylim(ax, [0.6 1.1]); set(ax, FontSize=13)
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / 15-200 km mean")
legend(ax, Location="southeast", FontSize=10, Interpreter="none")
title(ax, "Radar-distance imprint of the solution", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
bands = unique(validation.band, "stable");
bar(ax, categorical(bands, bands), reshape(validation.log_bias, numel(bands), []))
colororder(ax, colors)
yline(ax, 0, "k-", HandleVisibility="off")
set(ax, FontSize=13)
ylabel(ax, "Held-out log bias, log((obs+1)/(pred+1))")
legend(ax, string(cellfun(@(m) m.name, methods, UniformOutput=false)), Location="northeast", FontSize=11)
title(ax, "Radar cross-validation, sandbox (7 folds)", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sigma_weights.png"), Resolution=150)
close(fig)
