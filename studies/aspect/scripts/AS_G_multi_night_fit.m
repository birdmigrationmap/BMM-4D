%% Aspect G: the aspect function over many nights (dev)
% Pools the overlap cells of the us_east night (6 frames, 01-06 UTC, from
% AS_A_bins.m) and of every night finished by AS_F_multi_night_run.m, and
% asks what one night could not answer:
% 1. Per night: the "height" model (b2 = b2_per_km h, b1) and the
%    harmonics by height band. How much do they change between nights?
% 2. Shape: are cos 2 psi and cos psi enough in spring and autumn (sin
%    terms and the 4th harmonic near zero, 24-bin curves alike)?
% 3. Aspect against fixed radar patterns: the same fit with a fixed
%    azimuthal pattern per radar in the geographic frame
%    (as_overlap_fit_radar). With tracks that differ between nights the
%    two separate; if b2 and b1 hold, the one-night values were aspect.
% 4. What explains the change between nights: per night, height band and
%    region harmonics (stage 1), then a weighted regression on height,
%    speed, season and longitude with a between-stratum variance tau^2
%    (stage 2). tau is what no covariate explains.
% 5. Leave one night out: models fitted on the other nights (one
%    climatological function, with speed, by season) and the current
%    us_east function, scored on the left-out night (overlap disagreement
%    and aspect left) against no correction and a fit to that night itself
%    (two-fold by even / odd frames).
% 6. Below the overlap (< 0.8 km above the antennas): single-radar
%    harmonics with a fixed effect per radar, sweep and range band, now
%    over many nights (as_time_fit): the aspect turns with the track
%    between nights, the radar's own pattern does not.
% Decision rule (PLAN_multi_night.md): a fixed function (with covariates
% if they matter) when it gets at least 90% of the per-night fit's
% reduction of the disagreement on left-out nights and leaves |c2|, |c1|
% < 0.03 on most nights; otherwise a fit per night in the pipeline.
% Data: the shared night set (studies/nights): 17 nights plus the
% reference night 21-22 Oct 2026 (night 1), binned by as_night_set_bins
% into <data>/aspect/night_set/<name>/frame_XX.mat. The fit uses the
% frames on the hour (every other frame), so each night counts by its
% length, not its frame rate, and memory stays moderate.
% Runs on the cluster (hpc/submit.sh); tables and the figure go to
% <data>/aspect/night_set/results (the job's code snapshot is read only).
% Also: A1 of PLAN_robust_aspect.md (section 3b: can a design separate the
% aspect from radar patterns? simulated log ratios refitted) and B1
% (section 3c: radar patterns from rain overlaps, then the bird fit with
% those patterns removed).

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"), fullfile(project_root, "examples"))
data_root = bmm4d.util.local_data_dir();
multi_dir = fullfile(data_root, "aspect", "night_set");
figure_dir = fullfile(multi_dir, "results");
if ~isfolder(figure_dir)
    mkdir(figure_dir)
end
us = case_settings("us_east");
n_radar = numel(us.radars);
rng(1)

%% Read the nights
night_list = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string");
night_list = night_list(startsWith(night_list.name, "night_us_"), :);
night_list.season(night_list.name == "night_us_20261021") = "autumn";
% Reference night first (night 1)
night_list = [night_list(night_list.name == "night_us_20261021", :); night_list(night_list.name ~= "night_us_20261021", :)];
night_list = night_list(arrayfun(@(name) isfile(fullfile(multi_dir, name, "done.txt")), night_list.name), :);
night_table = table(night_list.name, night_list.season, datetime(night_list.date), VariableNames=["name" "season" "date"]);
n_night = height(night_table);
pair_parts = {};
volume_parts = {};
for k = 1:n_night
    manifest = load(fullfile(data_root, "prod", night_table.name(k), "manifest.mat"), "radars", "times");
    assert(isequal(string(manifest.radars(:)), string(us.radars(:))), "The radar order must be the us_east one.")
    hourly = find(minute(manifest.times) == 0);
    files = fullfile(multi_dir, night_table.name(k), compose("frame_%02d.mat", hourly(:)));
    files = files(arrayfun(@isfile, files));
    for i = 1:numel(files)
        loaded = load(files(i), "bins", "volumes");
        if isempty(loaded.bins) || height(loaded.bins) == 0
            continue
        end
        % Frames numbered per night: 100 k + frame
        p = as_overlap_pairs(loaded.bins);
        p = removevars(p, p.Properties.VariableNames(startsWith(p.Properties.VariableNames, ["ell_" "fix_"])));
        p.frame = 100 * k + double(p.frame);
        p.night = repmat(k, height(p), 1);
        pair_parts{end + 1} = p; %#ok<SAGROW>
        vols = rmfield(loaded.volumes, intersect(fieldnames(loaded.volumes), ["fix_track_deg" "fix_amplitude"]));
        for j = 1:numel(vols)
            vols(j).frame = 100 * k + vols(j).frame;
            vols(j).night = k;
        end
        volume_parts{end + 1} = orderfields(vols); %#ok<SAGROW>
    end
end
pairs = vertcat(pair_parts{:});
volumes = [volume_parts{:}];
clear pair_parts volume_parts loaded
radar_height = nan(n_radar, 1);
radar_height([volumes.radar]) = [volumes.height];
radar_lon = nan(n_radar, 1);
radar_lon([volumes.radar]) = [volumes.lon];
D = as_overlap_design(pairs, "raw");
D.night = pairs.night(D.keep);
D.h_agl = D.height_asl_m - (radar_height(D.radar_a) + radar_height(D.radar_b)) / 2;
D.lon = (radar_lon(D.radar_a) + radar_lon(D.radar_b)) / 2;
clear pairs
row_season = night_table.season(D.night);
fprintf("%d nights, %d cell pairs, %d radar pairs, %d volumes\n", n_night, numel(D.y), max(D.pair), numel(volumes))

%% 1. Per night: height model and harmonics by height band
height_edges = [750 1000 1500 2000 2500 3500];
per_night = table();
per_night_height = table();
night_curves = cell(n_night, 1);
for k = 1:n_night
    rows = D.night == k;
    model = as_fit_aspect_model(D, rows, radar_height, "height");
    strata = {};
    labels = strings(0);
    for i = 1:numel(height_edges) - 1
        in_band = rows & D.h_agl >= height_edges(i) & D.h_agl < height_edges(i + 1);
        if nnz(in_band) >= 2000
            strata{end + 1} = in_band; %#ok<SAGROW>
            labels(end + 1) = string((height_edges(i) + height_edges(i + 1)) / 2); %#ok<SAGROW>
        end
    end
    H = as_fit_strata(D, strata, labels);
    H.height_m = double(H.stratum);
    H.night = repmat(k, height(H), 1);
    per_night_height = [per_night_height; H]; %#ok<AGROW>
    night_curves{k} = as_aspect_curve(D, rows, 24);
    track_vector = mean(exp(1i * deg2rad(D.track(rows))));
    per_night = [per_night; table(k, night_table.name(k), night_table.season(k), nnz(rows), ...
        model.b2_per_km, model.se(1), model.b1, model.se(2), mod(rad2deg(angle(track_vector)), 360), ...
        abs(mean(exp(2i * deg2rad(D.track(rows))))), median(D.speed(rows)), ...
        VariableNames=["night" "name" "season" "n" "b2_per_km" "b2_per_km_se" "b1" "b1_se" "track_deg" ...
        "axis_coherence" "speed_ms"])]; %#ok<AGROW>
end
disp(per_night)
writetable(per_night, fullfile(figure_dir, "AS_G_per_night.csv"))
writetable(per_night_height, fullfile(figure_dir, "AS_G_per_night_height.csv"))

%% 2. Shape: harmonics and 24-bin curves by season
shape_strata = {true(size(D.y)), row_season == "spring", row_season == "autumn"};
shape = as_fit_strata(D, shape_strata, ["all" "spring" "autumn"]);
disp(shape(:, ["stratum" "n" "c2" "c2_se" "s2" "s2_se" "c1" "c1_se" "s1" "s1_se" "c4" "c4_se"]))
writetable(shape, fullfile(figure_dir, "AS_G_shape.csv"))
season_curves = struct("spring", as_aspect_curve(D, row_season == "spring", 24), ...
    "autumn", as_aspect_curve(D, row_season == "autumn", 24));

%% 3. Aspect against fixed radar patterns
all_rows = true(size(D.y));
fit_plain = as_overlap_fit_radar(D, all_rows, radar_height, n_radar, false);
fit_radar = as_overlap_fit_radar(D, all_rows, radar_height, n_radar, true);
fit_radar_night0 = as_overlap_fit_radar(D, D.night == 1, radar_height, n_radar, true);
% With a prior of 0.1 log units on each radar harmonic (the free patterns
% of weakly constrained radars are not physical; see as_overlap_fit_radar)
fit_radar_prior = as_overlap_fit_radar(D, all_rows, radar_height, n_radar, true, 0.1);
fit_plain_night0 = as_overlap_fit_radar(D, D.night == 1, radar_height, n_radar, false);
fits = {fit_plain, fit_radar, fit_radar_prior, fit_plain_night0, fit_radar_night0};
radar_test = table(["all nights" "all nights, radar patterns" "all nights, radar patterns, prior 0.1" ...
    "reference night" "reference night, radar patterns"]', cellfun(@(f) f.beta(1), fits)', cellfun(@(f) f.se(1), fits)', ...
    cellfun(@(f) f.beta(2), fits)', cellfun(@(f) f.se(2), fits)', ...
    VariableNames=["fit" "b2_per_km" "b2_per_km_se" "b1" "b1_se"]);
disp(radar_test)
writetable(radar_test, fullfile(figure_dir, "AS_G_radar_patterns_test.csv"))
% Fixed pattern of each radar (geographic frame): 1st and 2nd harmonics
radar_patterns = table(string(us.radars(:)), fit_radar.radar(:, 1), fit_radar.radar(:, 2), fit_radar.radar(:, 3), ...
    fit_radar.radar(:, 4), hypot(fit_radar.radar(:, 1), fit_radar.radar(:, 2)), ...
    hypot(fit_radar.radar(:, 3), fit_radar.radar(:, 4)), max(fit_radar.radar_se, [], 2), ...
    VariableNames=["radar" "p1" "q1" "p2" "q2" "amplitude1" "amplitude2" "max_se"]);
writetable(radar_patterns, fullfile(figure_dir, "AS_G_radar_patterns.csv"))

%% 3b. Can the design separate aspect from radar patterns? (PLAN_robust_aspect A1)
% Simulated log ratios on the real design (same cells, pairs, frames and
% weights): (a) aspect only, the current shared function; (b) radar
% patterns only, those of the joint fit on the reference night; (c) both.
% Noise from the residual mismatch and the estimator variance. Each is
% refitted with and without radar terms, on the reference night and on all
% nights. If (a) keeps its amplitude with radar terms, the design separates
% the two; if it collapses, the data cannot decide.
current_model = struct("b2", 0, "b2_per_km", 0.092, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, "b1", -0.21, ...
    "b1_ramp_m", 1000);
[log_ta, log_tb] = as_pair_templates(D, volumes, current_model);
aspect_part = log_ta - log_tb;
pattern_of = @(P, az, radar) sum([cosd(az) sind(az) cosd(2 * az) sind(2 * az)] .* P(radar, :), 2);
P_night0 = fit_radar_night0.radar;
P_night0(~isfinite(P_night0)) = 0;
az_a = mod(D.psi_a + D.track, 360);
az_b = mod(D.psi_b + D.track, 360);
pattern_part = pattern_of(P_night0, az_a, D.radar_a) - pattern_of(P_night0, az_b, D.radar_b);
var_noise = D.var_est;
var_noise(~isfinite(var_noise)) = median(var_noise, "omitnan");
noise = randn(size(D.y)) .* sqrt(var_noise + fit_plain.sigma_mismatch ^ 2);
scenarios = struct("aspect", aspect_part, "patterns", pattern_part, "both", aspect_part + pattern_part);
identifiability = table();
for scope = ["reference night" "all nights"]
    rows = all_rows;
    if scope == "reference night"
        rows = D.night == 1;
    end
    for name = string(fieldnames(scenarios))'
        Ds = D;
        Ds.y = scenarios.(name) + noise;
        for radar_terms = [false true]
            f = as_overlap_fit_radar(Ds, rows, radar_height, n_radar, radar_terms);
            identifiability = [identifiability; table(scope, name, radar_terms, f.beta(1), f.se(1), f.beta(2), ...
                f.se(2), VariableNames=["scope" "truth" "radar_terms" "b2_per_km" "b2_per_km_se" "b1" "b1_se"])]; %#ok<AGROW>
        end
    end
end
clear Ds log_ta log_tb
disp(identifiability)
writetable(identifiability, fullfile(figure_dir, "AS_G_identifiability.csv"))

%% 3c. Radar patterns from rain (PLAN_robust_aspect B1)
% Rain has no aspect: the radar patterns of the rain overlaps (raw
% reflectivity, as_rain_bins; dummy track to the north, so the "aspect"
% columns are a pattern common to all radars, expected near zero) are
% compared with those of the joint bird fit, then removed from the bird
% log ratios before the plain fit.
rain_names = ["rain_us_20240411" "rain_us_20250405" "rain_us_20211028"];
rain_parts = {};
for k = 1:numel(rain_names)
    listing = dir(fullfile(multi_dir, rain_names(k), "frame_*.mat"));
    for i = 1:numel(listing)
        loaded = load(fullfile(listing(i).folder, listing(i).name), "bins");
        if isempty(loaded.bins) || height(loaded.bins) == 0
            continue
        end
        p = as_overlap_pairs(loaded.bins);
        p.frame = 100 * k + double(p.frame);
        rain_parts{end + 1} = p; %#ok<SAGROW>
    end
end
[rain_fit, rain_compare, rain_test] = deal([]);
if isempty(rain_parts)
    disp("No rain bins: section 3c skipped")
else
    R = as_overlap_design(vertcat(rain_parts{:}), "raw");
    rain_fit = as_overlap_fit_radar(R, true(size(R.y)), radar_height, n_radar, true);
    fprintf("Rain: %d cell pairs, %d radar pairs, common pattern c2 %.3f (%.3f) per km, c1 %.3f (%.3f)\n", ...
        rain_fit.n, rain_fit.n_pair, rain_fit.beta(1), rain_fit.se(1), rain_fit.beta(2), rain_fit.se(2))
    % Agreement of the patterns, per harmonic, over radars with both estimates
    harmonic_names = ["p1" "q1" "p2" "q2"];
    rain_compare = table();
    for source = ["joint fit, all nights" "joint fit, reference night"]
        P_bird = fit_radar.radar;
        if source == "joint fit, reference night"
            P_bird = fit_radar_night0.radar;
        end
        both = all(isfinite(P_bird), 2) & all(isfinite(rain_fit.radar), 2);
        for h = 1:4
            x = rain_fit.radar(both, h);
            y = P_bird(both, h);
            rain_compare = [rain_compare; table(source, harmonic_names(h), nnz(both), corr(x, y), ...
                (x' * y) / (x' * x), std(x), std(y), VariableNames=["bird_fit" "harmonic" "n_radar" "r" ...
                "slope" "sd_rain" "sd_bird"])]; %#ok<AGROW>
        end
    end
    disp(rain_compare)
    writetable(rain_compare, fullfile(figure_dir, "AS_G_rain_patterns_compare.csv"))
    rain_patterns = table(string(us.radars(:)), rain_fit.radar(:, 1), rain_fit.radar(:, 2), rain_fit.radar(:, 3), ...
        rain_fit.radar(:, 4), max(rain_fit.radar_se, [], 2), VariableNames=["radar" "p1" "q1" "p2" "q2" "max_se"]);
    writetable(rain_patterns, fullfile(figure_dir, "AS_G_rain_patterns.csv"))
    % Bird fit with the rain patterns removed
    P_rain = rain_fit.radar;
    P_rain(~isfinite(P_rain)) = 0;
    Dr = D;
    Dr.y = D.y - (pattern_of(P_rain, az_a, D.radar_a) - pattern_of(P_rain, az_b, D.radar_b));
    rain_test = table();
    for scope = ["reference night" "all nights"]
        rows = all_rows;
        if scope == "reference night"
            rows = D.night == 1;
        end
        for removed = [false true]
            Duse = D;
            if removed
                Duse = Dr;
            end
            f = as_overlap_fit_radar(Duse, rows, radar_height, n_radar, false);
            rain_test = [rain_test; table(scope, removed, f.beta(1), f.se(1), f.beta(2), f.se(2), f.sigma_mismatch, ...
                VariableNames=["scope" "rain_patterns_removed" "b2_per_km" "b2_per_km_se" "b1" "b1_se" ...
                "sigma_mismatch"])]; %#ok<AGROW>
        end
    end
    clear Dr Duse R
    disp(rain_test)
    writetable(rain_test, fullfile(figure_dir, "AS_G_rain_test.csv"))
end

%% 4. What explains the change between nights (two stages)
region_edges = [-inf -92 -84 inf];
stage1 = table();
for k = 1:n_night
    for i_h = 1:numel(height_edges) - 1
        for i_r = 1:numel(region_edges) - 1
            rows = D.night == k & D.h_agl >= height_edges(i_h) & D.h_agl < height_edges(i_h + 1) & ...
                D.lon >= region_edges(i_r) & D.lon < region_edges(i_r + 1);
            if nnz(rows) < 5000 || numel(unique(D.pair(rows))) < 5
                continue
            end
            fit = as_overlap_fit(D, ["c2" "c1" D.names(D.beam_columns)], rows);
            stage1 = [stage1; table(k, night_table.season(k), median(D.h_agl(rows)), median(D.speed(rows)), ...
                mean(D.lon(rows)), fit.beta(1), fit.se(1), fit.beta(2), fit.se(2), nnz(rows), ...
                VariableNames=["night" "season" "h_m" "speed_ms" "lon" "c2" "c2_se" "c1" "c1_se" "n"])]; %#ok<AGROW>
        end
    end
end
covariates = [ones(height(stage1), 1), (stage1.h_m - 2000) / 1000, stage1.speed_ms - 9, ...
    double(stage1.season == "spring"), (stage1.lon + 88) / 10];
covariate_names = ["intercept" "per km above 2 km" "per m/s above 9" "spring" "per 10 deg lon"];
meta = struct("c2", meta_regression(stage1.c2, stage1.c2_se, covariates, covariate_names), ...
    "c1", meta_regression(stage1.c1, stage1.c1_se, covariates, covariate_names), ...
    "c2_height_only", meta_regression(stage1.c2, stage1.c2_se, covariates(:, 1:2), covariate_names(1:2)), ...
    "c1_constant", meta_regression(stage1.c1, stage1.c1_se, covariates(:, 1), covariate_names(1)));
for name = string(fieldnames(meta))'
    fprintf("%s: tau %.3f\n", name, meta.(name).tau)
    disp(meta.(name).coefficients)
end
writetable(stage1, fullfile(figure_dir, "AS_G_stage1.csv"))

%% 5. Leave one night out
current = struct("b2", 0, "b2_per_km", 0.092, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, "b1", -0.21, ...
    "b1_ramp_m", 1000);
lono = table();
for k = 1:n_night
    test = D.night == k;
    train = ~test;
    Dk = subset_design(D, test);
    models = struct("current", current, ...
        "climatology", as_fit_aspect_model(D, train, radar_height, "height"), ...
        "climatology_speed", as_fit_aspect_model(D, train, radar_height, "height_speed"));
    % Same season only when another night of that season is there
    same_season = train & row_season == night_table.season(k);
    if nnz(same_season) >= 1e5
        models.season = as_fit_aspect_model(D, same_season, radar_height, "height");
    else
        models.season = models.climatology;
    end
    for name = ["climatology" "climatology_speed" "season"]
        models.(name).b1_ramp_m = 1000;
    end
    raw_fit = as_overlap_fit(Dk, ["c2" "c1" Dk.names(Dk.beam_columns)]);
    weights = 1 ./ (Dk.var_est + raw_fit.sigma_mismatch ^ 2);
    y_methods = struct("none", Dk.y);
    for name = string(fieldnames(models))'
        [log_ta, log_tb] = as_pair_templates(Dk, volumes, models.(name));
        y_methods.(name) = Dk.y - log_ta + log_tb;
    end
    % The night's own fit, two-fold by even and odd frames
    y_methods.own_night = nan(size(Dk.y));
    [~, ~, frame_rank] = unique(Dk.frame);
    even = mod(frame_rank, 2) == 0;
    for fold = 1:2
        fit_rows = even == (fold == 1);
        own = as_fit_aspect_model(Dk, fit_rows, radar_height, "height");
        own.b1_ramp_m = 1000;
        [log_ta, log_tb] = as_pair_templates(Dk, volumes, own);
        y_methods.own_night(~fit_rows) = Dk.y(~fit_rows) - log_ta(~fit_rows) + log_tb(~fit_rows);
    end
    common = isfinite(y_methods.own_night);
    for name = string(fieldnames(y_methods))'
        s = as_overlap_score(Dk, y_methods.(name), common, weights);
        lono = [lono; table(k, night_table.name(k), name, s.disagreement, s.c2, s.c2_se, s.c1, s.c1_se, ...
            VariableNames=["night" "name" "method" "disagreement" "c2_left" "c2_se" "c1_left" "c1_se"])]; %#ok<AGROW>
    end
end
% Reduction of the disagreement relative to no correction, and as a share of the night's own fit
none = lono(lono.method == "none", ["night" "disagreement"]);
own = lono(lono.method == "own_night", ["night" "disagreement"]);
[~, i_none] = ismember(lono.night, none.night);
[~, i_own] = ismember(lono.night, own.night);
lono.reduction = 1 - lono.disagreement ./ none.disagreement(i_none);
lono.share_of_own = lono.reduction ./ (1 - own.disagreement(i_own) ./ none.disagreement(i_none));
lono_summary = groupsummary(lono, "method", "median", ["reduction" "share_of_own"]);
lono_summary.share_aspect_left_small = splitapply(@(c2, c1) mean(abs(c2) < 0.03 & abs(c1) < 0.03), ...
    lono.c2_left, lono.c1_left, findgroups(lono.method));
disp(lono_summary)
writetable(lono, fullfile(figure_dir, "AS_G_leave_one_night_out.csv"))
writetable(lono_summary, fullfile(figure_dir, "AS_G_leave_one_night_out_summary.csv"))
% The us_east night is in-sample for "current" (its parameters come from it)

%% 6. Below the overlap: single-radar harmonics over many nights
cells = as_cell_harmonics(volumes, "raw");
cells.night = floor(cells.frame / 100);
single_edges = [0 250 500 750 1000 1500 2000 2500 3500];
single = table();
for i = 1:numel(single_edges) - 1
    in_band = cells.band >= 2 & cells.height_m >= single_edges(i) & cells.height_m < single_edges(i + 1);
    for scope = ["all nights" "us_east night"]
        rows = in_band & (scope == "all nights" | cells.night == 1);
        if nnz(rows) < 200
            continue
        end
        tf = as_time_fit(cells, rows, 100);
        single = [single; table(scope, (single_edges(i) + single_edges(i + 1)) / 2, tf.n, tf.b2, tf.b2_se, ...
            tf.b2_turned, tf.b1, tf.b1_se, tf.b1_turned, tf.track_spread_deg, VariableNames=["scope" "height_m" ...
            "n" "b2" "b2_se" "b2_turned" "b1" "b1_se" "b1_turned" "track_spread_deg"])]; %#ok<AGROW>
    end
end
disp(single)
writetable(single, fullfile(figure_dir, "AS_G_single_radar.csv"))

%% Decision summary
climatology = lono_summary(lono_summary.method == "climatology", :);
fprintf("\nPooled (all nights): b2_per_km %.3f (%.3f), b1 %.3f (%.3f); with radar patterns: %.3f (%.3f), %.3f (%.3f)\n", ...
    fit_plain.beta(1), fit_plain.se(1), fit_plain.beta(2), fit_plain.se(2), fit_radar.beta(1), fit_radar.se(1), ...
    fit_radar.beta(2), fit_radar.se(2))
fprintf("Between-night sd not explained by covariates: c2 %.3f, c1 %.3f\n", meta.c2.tau, meta.c1.tau)
fprintf("Climatology on left-out nights: %.0f%% of the own-night reduction, aspect left small on %.0f%% of nights\n", ...
    100 * climatology.median_share_of_own, 100 * climatology.share_aspect_left_small)
if climatology.median_share_of_own >= 0.9 && climatology.share_aspect_left_small >= 0.75
    disp("-> a fixed function is enough (check the covariates of section 4 before choosing its form)")
else
    disp("-> fit the aspect function per night in the pipeline (overlap fit after cleaning)")
end

%% Save
save(fullfile(multi_dir, "multi_night_results.mat"), "night_table", "per_night", "per_night_height", ...
    "night_curves", "shape", "season_curves", "radar_test", "radar_patterns", "stage1", "meta", "lono", ...
    "lono_summary", "single", "fit_plain", "fit_radar", "fit_radar_night0", "identifiability", "rain_fit", ...
    "rain_compare", "rain_test")

%% Plot results
% (a, b) side-on and front/back harmonics by height, one line per night
% (spring green, autumn brown), the pooled fit in black; (c) reduction of
% the overlap disagreement on left-out nights by method; (d) below and in
% the overlap: single-radar (many nights) against overlap harmonics.
fig = bmm4d.plot.create_export_figure(3000, 2200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 2, 2, TileSpacing="tight", Padding="tight");
season_color = struct("spring", [0.30 0.60 0.30], "autumn", [0.65 0.40 0.15]);
h_line = (0:3500)' / 1000;
for i_panel = 1:2
    ax = nexttile(layout);
    hold(ax, "on");
    box(ax, "on")
    column = ["c2" "c1"];
    for k = 1:n_night
        rows = per_night_height.night == k;
        plot(ax, per_night_height.(column(i_panel))(rows), per_night_height.height_m(rows) / 1000, "-o", ...
            Color=season_color.(night_table.season(k)), MarkerSize=4)
    end
    if i_panel == 1
        plot(ax, fit_plain.beta(1) * h_line, h_line, "k-", LineWidth=2)
        xlabel(ax, "Side-on harmonic c_2")
        title(ax, "(a) Side-on, per night (spring green, autumn brown; black: all)", FontSize=13)
    else
        plot(ax, fit_plain.beta(2) * min(1, h_line), h_line, "k-", LineWidth=2)
        xlabel(ax, "Front/back harmonic c_1")
        title(ax, "(b) Front/back, per night", FontSize=13)
    end
    xline(ax, 0, ":")
    ylabel(ax, "Height above the antennas (km)")
    set(ax, FontSize=13)
end
ax = nexttile(layout);
hold(ax, "on");
box(ax, "on")
methods = ["current" "climatology" "climatology_speed" "season" "own_night"];
for i_m = 1:numel(methods)
    rows = lono.method == methods(i_m);
    scatter(ax, i_m + 0.1 * randn(nnz(rows), 1), 100 * lono.reduction(rows), 40, [0.5 0.5 0.5], "filled")
    plot(ax, i_m + [-0.3 0.3], 100 * median(lono.reduction(rows)) * [1 1], "k-", LineWidth=2)
end
xticks(ax, 1:numel(methods))
xticklabels(ax, ["us\_east fit" "all other nights" "+ speed" "same season" "own night"])
ylabel(ax, "Disagreement reduction (%)")
title(ax, "(c) Left-out nights: fitted on", FontSize=13)
set(ax, FontSize=13)
ax = nexttile(layout);
hold(ax, "on");
box(ax, "on")
rows = single.scope == "all nights";
errorbar(ax, single.b2(rows), single.height_m(rows) / 1000, single.b2_se(rows), "horizontal", "o-", ...
    Color=[0.85 0.33 0.10], DisplayName="c_2, single radars, all nights")
errorbar(ax, single.b1(rows), single.height_m(rows) / 1000, single.b1_se(rows), "horizontal", "s-", ...
    Color=[0.00 0.45 0.74], DisplayName="c_1, single radars, all nights")
pooled = groupsummary(per_night_height, "height_m", "mean", ["c2" "c1"]);
plot(ax, pooled.mean_c2, pooled.height_m / 1000, "o--", Color=[0.85 0.33 0.10], DisplayName="c_2, overlap (mean of nights)")
plot(ax, pooled.mean_c1, pooled.height_m / 1000, "s--", Color=[0.00 0.45 0.74], DisplayName="c_1, overlap (mean of nights)")
xline(ax, 0, ":", HandleVisibility="off")
xlim(ax, "padded")
legend(ax, Location="best")
title(ax, "(d) Single radars over many nights against the overlap", FontSize=13)
xlabel(ax, "Harmonic amplitude")
ylabel(ax, "Height above the antennas (km)")
set(ax, FontSize=13)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AS_G_multi_night.png"), Resolution=300)
close(fig)

%% Local functions
function E = subset_design(D, rows)
% The pairs ROWS of an as_overlap_design struct (fields with one value per pair)
E = D;
n = numel(D.y);
for name = string(fieldnames(D))'
    value = D.(name);
    if size(value, 1) == n && name ~= "keep"
        E.(name) = value(rows, :);
    end
end
end

function out = meta_regression(y, se, X, names)
% Weighted regression of stratum estimates Y (standard errors SE) on X,
% with a between-stratum variance tau^2 by the method of moments
% (DerSimonian-Laird with covariates): what the covariates do not explain.
w = 1 ./ se .^ 2;
A = X' * (w .* X);
beta = A \ (X' * (w .* y));
Q = sum(w .* (y - X * beta) .^ 2);
df = numel(y) - size(X, 2);
trace_term = sum(w) - trace(A \ (X' * (w .^ 2 .* X)));
tau2 = max(0, (Q - df) / trace_term);
w = 1 ./ (se .^ 2 + tau2);
A = X' * (w .* X);
beta = A \ (X' * (w .* y));
out = struct("tau", sqrt(tau2), "coefficients", table(names(:), beta, sqrt(diag(inv(A))), ...
    VariableNames=["term" "estimate" "se"]));
end
