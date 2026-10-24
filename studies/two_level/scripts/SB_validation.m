%% Sandbox: radar cross-validation (prototype step 5)
% Each radar's observations are predicted by a solve that did not see
% them. The 28 sandbox radars are split into 7 folds of 4 radars far
% apart (one per longitude band), so each fold withholds 4 radars; every
% radar is withheld once (leave-four-radars-out, 7 solves per method
% instead of 28). The withheld radars are also removed from the radar
% prior and the regional shape; the advection wind still uses all radars.
% Predictions use each model's own observation operator (the two-level
% near-range recovery a(r) of that fold; none for the current solver).
% Scores on the withheld rows, by range band:
% - relative L2 misfit (observed no-echo gates count as 0);
% - log bias: mean log((obs + 1) / (pred + 1)) on echo rows (> 0: the
%   model predicts too little), and its RMS (multiplicative error).
% Out-of-sample predictions are cached per method in the sandbox data
% folder, so the script can be rerun section by section.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "sandbox"), fullfile(project_root, "studies", "two_level", "functions"))
figure_dir = fullfile(project_root, "studies", "two_level", "figures", "sandbox");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end
cache_dir = fullfile(S.c.data_dir, "validation");
if ~exist(cache_dir, "dir")
    mkdir(cache_dir)
end

%% Folds: 4 longitude bands of 7 radars, fold j takes the j-th radar (by latitude) of each band
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
radar_names = S.radar_selected.name;
disp(table((1:n_fold)', arrayfun(@(f) strjoin(radar_names(fold == f), " "), (1:n_fold)'), VariableNames=["fold" "withheld"]))

%% Methods ----------------------------------------------------------------
logistic28 = @(r, e, z) (1 ./ (1 + exp(-(r - 28000) / 3000))) .* (~z + z ./ (1 + exp((r - 75000) / 10000)));
base = {"fine_prior", 0.2, "shape_prior", 0.03, "shape_length_m", 150000, "level_smoothness", 0.03};
before_far_fix = {"fine_amplitude", "relative", "advection_on", "density"};
methods = {
    struct("name", "current, no cuts, prior 0.3")
    struct("name", "current, logistic, prior 0.3", "row_weight", logistic28)
    struct("name", "current, logistic, prior 0.5", "row_weight", logistic28, "prior", 0.5)
    struct("name", "2L K3, logistic, no obs model", "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
        "row_weight", logistic28, "two_level", struct(base{:}, before_far_fix{:}, "near_range", "off"))
    struct("name", "2L K3, near elev (default)", "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct(base{:}, before_far_fix{:}))
    struct("name", "2L K6, near elev", "solver", "two_level", "n_pattern", 6, "regional_patterns", 5, ...
        "two_level", struct(base{:}, before_far_fix{:}))
    struct("name", "2L K3, near elev + error 0.4", "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct(base{:}, before_far_fix{:}, "error_cv", 0.4, "error_floor", 15))
    struct("name", "2L K3, far fix", "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
        "two_level", struct(base{:}))
    struct("name", "2L K3, far fix, no advection", "solver", "two_level", "n_pattern", 3, "regional_patterns", 2, ...
        "advection", 0, "two_level", struct(base{:}))};
% "far fix" (24 October 2026, now the solve_two_level defaults): unweighted
% amplitude of the fine deviation and advection of the deviation only.

%% Out-of-sample predictions (cached per method) ----------------------------------
observed = full(S.observations.values);
observed(observed > 0 & observed < 1e-6) = 0;            % no-echo gates
n_time = size(observed, 2);
F_active = S.observations.F(:, S.active_id);
range_m = double(S.observations.row_range_m);
elev_deg = double(S.observations.row_elev_deg);
predictions = cell(numel(methods), 1);
for i_m = 1:numel(methods)
    cache_file = fullfile(cache_dir, matlab.lang.makeValidName(methods{i_m}.name) + ".mat");
    if exist(cache_file, "file")
        loaded = load(cache_file, "predicted", "seconds");
        predictions{i_m} = loaded;
        continue
    end
    predicted = zeros(size(observed), "single");
    seconds = 0;
    for i_f = 1:n_fold
        S_fold = sandbox_fold(S, find(fold == i_f));
        option = methods{i_m};
        option.panel = false;
        r = two_level_run(S_fold, option);
        rows = find(S_fold.withheld_rows);
        a = ones(numel(rows), 1);
        if isfield(r.info, "near_curve") && r.info.controls.near_range ~= "off"
            curve = r.info.near_curve;
            group = discretize(elev_deg(rows), [curve.elev_from; curve.elev_to(end)]);
            group(isnan(group)) = height(curve);
            a = 1 - exp(-(range_m(rows) ./ curve.r0_m(group)) .^ curve.k(group));
        end
        predicted(rows, :) = single(a .* (F_active(rows, :) * double(r.field)));
        seconds = seconds + r.seconds;
        fprintf("[validation] %-30s fold %d (%s) %4.0f s\n", methods{i_m}.name, i_f, strjoin(radar_names(fold == i_f), " "), r.seconds)
    end
    save(cache_file, "predicted", "seconds")
    predictions{i_m} = struct("predicted", predicted, "seconds", seconds);
end

%% Scores by range band ---------------------------------------------------------
bands = struct("name", ["near 0-30 km" "trusted 30-120 km" "far > 120 km"], "from", [0 30000 120000], "to", [30000 120000 inf]);
is_observed = full(S.observations.values ~= 0);
is_echo = is_observed & observed > 0;
scores = table();
for i_m = 1:numel(methods)
    predicted = double(predictions{i_m}.predicted);
    row = table(string(methods{i_m}.name), predictions{i_m}.seconds / n_fold, VariableNames=["method" "seconds_per_fold"]);
    for i_b = 1:3
        in_band = repmat(range_m >= bands.from(i_b) & range_m < bands.to(i_b), 1, n_time);
        selected = is_observed & in_band;
        echo_selected = is_echo & in_band;
        log_ratio = log((observed(echo_selected) + 1) ./ (max(predicted(echo_selected), 0) + 1));
        tag = ["near" "trusted" "far"];
        row.(tag(i_b) + "_misfit") = norm(predicted(selected) - observed(selected)) / norm(observed(selected));
        row.(tag(i_b) + "_log_bias") = mean(log_ratio);
        row.(tag(i_b) + "_log_rms") = sqrt(mean(log_ratio .^ 2));
    end
    scores = [scores; row]; %#ok<AGROW>
end
disp(scores)
writetable(scores, fullfile(figure_dir, "validation_scores.csv"))

%% Scores per withheld radar (trusted band) ---------------------------------------
radar_misfit = nan(n_radar, numel(methods));
radar_log_bias = nan(n_radar, numel(methods));
trusted_rows = range_m >= 30000 & range_m < 120000;
for i_m = 1:numel(methods)
    predicted = double(predictions{i_m}.predicted);
    for i_r = 1:n_radar
        rows = S.observations.row_radar_id == i_r & trusted_rows;
        selected = full(is_observed(rows, :));
        obs_r = observed(rows, :);
        pred_r = predicted(rows, :);
        radar_misfit(i_r, i_m) = norm(pred_r(selected) - obs_r(selected)) / norm(obs_r(selected));
        echo_r = is_echo(rows, :);
        radar_log_bias(i_r, i_m) = mean(log((obs_r(echo_r) + 1) ./ (max(pred_r(echo_r), 0) + 1)));
    end
end
per_radar = array2table(round(radar_misfit, 3), VariableNames=cellfun(@(m) char(m.name), methods, UniformOutput=false), RowNames=cellstr(radar_names));
disp(per_radar)
writetable(per_radar, fullfile(figure_dir, "validation_per_radar_misfit.csv"), WriteRowNames=true)

%% Plot results -----------------------------------------------------------------
fig = create_export_figure(3600, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
method_names = categorical(scores.method, scores.method);
ax = nexttile(t);
barh(ax, method_names, [scores.near_misfit scores.trusted_misfit scores.far_misfit])
legend(ax, ["near 0-30 km" "trusted 30-120 km" "far > 120 km"], Location="southoutside", Orientation="horizontal")
xlabel(ax, "Relative misfit on withheld radars"); set(ax, FontSize=12); grid(ax, "on")
title(ax, "Out-of-sample misfit", FontWeight="normal")
ax = nexttile(t);
barh(ax, method_names, [scores.near_log_bias scores.trusted_log_bias scores.far_log_bias])
xline(ax, 0, "k-")
xlabel(ax, "Log bias (> 0: predicted too low)"); set(ax, FontSize=12, YTickLabel=[]); grid(ax, "on"); xlim(ax, [-1.2 0.3])
title(ax, "Out-of-sample log bias, log((obs+1)/(pred+1))", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
reference = 2;   % current solver, logistic, prior 0.3
for i_m = 1:numel(methods)
    plot(ax, radar_misfit(:, reference), radar_misfit(:, i_m), "o", LineWidth=1.5, DisplayName=methods{i_m}.name)
end
plot(ax, [0 1.2], [0 1.2], "k--", HandleVisibility="off")
xlabel(ax, "Trusted misfit, " + methods{reference}.name); ylabel(ax, "Trusted misfit, each method")
legend(ax, Location="southoutside", FontSize=10); set(ax, FontSize=12); axis(ax, "equal"); xlim(ax, [0 1.2]); ylim(ax, [0 1.2])
title(ax, "Per withheld radar", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "validation.png"), Resolution=150)
close(fig)
