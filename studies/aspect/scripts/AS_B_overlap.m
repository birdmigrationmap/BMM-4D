%% Aspect B: the aspect function from overlapping radars, and the corrections compared (dev)
% Where two radars see the same 4 km x 4 km x 100 m cell at the same time
% (as_frame_bins, AS_A_bins.m), the log ratio of their densities does not
% depend on the density field, only on the two viewing directions:
%   log eta_a - log eta_b = log T(psi_a) - log T(psi_b) + pair offset + noise,
% psi = azimuth from the radar - track of the birds (mean of the two
% radars' velocity profiles at the cell height). Fixed effects per pair,
% frame and 500 m layer take calibration, per-volume scale and the
% different beam smoothing of the profile; a slope of log(range_a /
% range_b) per layer takes the beam broadening (as_overlap_design,
% as_overlap_fit). This gives the aspect function free of real gradients.
% Then:
% - the harmonics in subsets (height, speed, region, time, beam matching);
% - the aspect function without a shape (24 bins of psi);
% - single-radar harmonics (population means: real gradients average out
%   only if they have no preferred direction relative to the track), the
%   variance of the side-on harmonic against the 45 deg one (how much the
%   amplitude really varies between volumes), and the time fit (the track
%   turning during the night at each radar, as_time_fit);
% - the corrections scored on the overlap (two-fold by even / odd frames:
%   shared parameters fitted on one half, scored on the other):
%   disagreement and aspect harmonics left; and on the single-radar
%   harmonics at 10-60 km (first harmonic = gradient across the disc).
% Outputs: <data>/aspect/overlap_results.mat, figures/AS_*.csv.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"))
data_root = bmm4d.util.local_data_dir();
bin_dir = fullfile(data_root, "aspect", "bins");
figure_dir = fullfile(project_root, "studies", "aspect", "figures");
manifest = load(fullfile(data_root, "prod", "us_east", "manifest.mat"));

%% Read the overlap pairs and the volume summaries
frame_files = dir(fullfile(bin_dir, "frame_*.mat"));
pair_parts = cell(numel(frame_files), 1);
volume_parts = cell(numel(frame_files), 1);
for i_frame = 1:numel(frame_files)
    loaded = load(fullfile(bin_dir, frame_files(i_frame).name));
    pair_parts{i_frame} = as_overlap_pairs(loaded.bins);
    volume_parts{i_frame} = loaded.volumes;
end
pairs = vertcat(pair_parts{:});
volumes = [volume_parts{:}];
radar_height = nan(numel(manifest.radars), 1);
radar_height([volumes.radar]) = [volumes.height];
radar_lon = nan(numel(manifest.radars), 1);
radar_lon([volumes.radar]) = [volumes.lon];
D = as_overlap_design(pairs, "raw");
D.h_agl = D.height_asl_m - (radar_height(D.radar_a) + radar_height(D.radar_b)) / 2;
fprintf("%d pairs of cells, %d radar pairs, %d volumes\n", numel(D.y), max(D.pair), numel(volumes))

%% Aspect harmonics: all pairs and subsets
all_rows = true(size(D.y));
strata = {all_rows};
labels = "all";
height_edges = [750 1000 1250 1500 2000 2500 3000 3500];
for i = 1:numel(height_edges) - 1
    strata{end + 1} = D.h_agl >= height_edges(i) & D.h_agl < height_edges(i + 1); %#ok<SAGROW>
    labels(end + 1) = sprintf("height %d-%d m", height_edges(i), height_edges(i + 1)); %#ok<SAGROW>
end
speed_edges = [0 6 10 14 30];
for i = 1:numel(speed_edges) - 1
    strata{end + 1} = D.speed >= speed_edges(i) & D.speed < speed_edges(i + 1); %#ok<SAGROW>
    labels(end + 1) = sprintf("speed %d-%d m/s", speed_edges(i), speed_edges(i + 1)); %#ok<SAGROW>
end
pair_lon = (radar_lon(D.radar_a) + radar_lon(D.radar_b)) / 2;
strata = [strata {pair_lon < -92, pair_lon >= -92 & pair_lon < -84, pair_lon >= -84, ...
    D.frame <= 8, D.frame > 8 & D.frame <= 16, D.frame > 16, ...
    abs(log(D.range_a ./ D.range_b)) < 0.05}];
labels = [labels "lon < -92" "lon -92 to -84" "lon >= -84" "23:30-03:00 UTC" "03:30-07:00 UTC" ...
    "07:30-12:00 UTC" "equal ranges (5%)"];
harmonics = as_fit_strata(D, strata, labels);
disp(harmonics(:, ["stratum" "n" "n_pair" "c2" "c2_se" "s2" "c1" "c1_se" "s1" "c4" "c4_se"]))
writetable(harmonics, fullfile(figure_dir, "AS_overlap_harmonics.csv"))

% Track or fixed geographic axis (the night's median track, 190 deg) together
az_a = mod(D.psi_a + D.track, 360);
az_b = mod(D.psi_b + D.track, 360);
E = D;
E.X = [D.X(:, [1 3]), -(cosd(2 * (az_a - 190)) - cosd(2 * (az_b - 190))), cosd(az_a - 190) - cosd(az_b - 190), ...
    D.X(:, D.beam_columns)];
E.names = ["c2" "c1" "c2_geographic" "c1_geographic" D.names(D.beam_columns)];
axis_fit = as_overlap_fit(E, E.names);
axis_test = rmfield(axis_fit, ["resid" "y_tilde" "weights"]);
disp([axis_test.names(1:4); compose("%.3f (%.3f)", [axis_test.beta(1:4); axis_test.se(1:4)]')'])

%% Aspect function without a shape (24 bins of 15 deg)
curve_raw = as_aspect_curve(D, all_rows, 24);

%% Single-radar harmonics: population means, variance, track turning in time
cells = as_cell_harmonics(volumes, "raw");
height_bins = [0 250 500 750 1000 1250 1500 2000 2500 3000 3500];
single_by_height = table();
for i = 1:numel(height_bins) - 1
    rows = cells.band >= 2 & cells.height_m >= height_bins(i) & cells.height_m < height_bins(i + 1);
    time_fit = as_time_fit(cells, rows, 100);
    single_by_height = [single_by_height; table((height_bins(i) + height_bins(i + 1)) / 2, nnz(rows), ...
        mean(cells.c2(rows)), std(cells.c2(rows)) / sqrt(nnz(rows)), mean(cells.c1(rows)), ...
        std(cells.c1(rows)) / sqrt(nnz(rows)), std(cells.c2(rows)), std(cells.s2(rows)), ...
        time_fit.b2, time_fit.b2_se, time_fit.b1, time_fit.b1_se, time_fit.track_spread_deg, ...
        VariableNames=["height_m" "n" "c2" "c2_se" "c1" "c1_se" "c2_sd" "s2_sd" "time_b2" "time_b2_se" ...
        "time_b1" "time_b1_se" "track_spread_deg"])]; %#ok<AGROW>
end
disp(single_by_height)
writetable(single_by_height, fullfile(figure_dir, "AS_single_radar_by_height.csv"))
overlap_by_height = harmonics(startsWith(harmonics.stratum, "height"), :);
overlap_by_height.height_m = (height_edges(1:end - 1) + height_edges(2:end))' / 2;

%% Shared models: two-fold by even and odd frames, then all frames
forms = ["symmetric" "constant" "height" "height_speed"];
even = mod(D.frame, 2) == 0;
fold_models = struct();
for fold = 1:2
    fit_rows = even == (fold == 1);
    for form = forms
        fold_models(fold).(form) = as_fit_aspect_model(D, fit_rows, radar_height, form);
    end
    % Front/back term ramped from zero at the antenna to full at 1 km (below
    % the lowest overlap; the single-radar front/back harmonic there is half)
    fold_models(fold).height_ramp = fold_models(fold).height;
    fold_models(fold).height_ramp.b1_ramp_m = 1000;
end
models = struct();
for form = forms
    models.(form) = as_fit_aspect_model(D, all_rows, radar_height, form);
end
models.height_ramp = models.height;
models.height_ramp.b1_ramp_m = 1000;
disp(models.height_ramp)

%% Score the corrections on the overlap (out of fold)
keep_id = find(D.keep);
y_methods = struct("none", D.y, ...
    "ellipse", log(double(pairs.ell_a(keep_id))) - log(double(pairs.ell_b(keep_id))), ...
    "per_volume", log(double(pairs.fix_a(keep_id))) - log(double(pairs.fix_b(keep_id))));
model_names = [forms "height_ramp"];
for name = model_names
    y_methods.(name) = nan(size(D.y));
    for fold = 1:2
        score_rows = even ~= (fold == 1);
        [log_ta, log_tb] = as_pair_templates(D, volumes, fold_models(fold).(name));
        y_methods.(name)(score_rows) = D.y(score_rows) - log_ta(score_rows) + log_tb(score_rows);
    end
end
common = isfinite(y_methods.ellipse) & isfinite(y_methods.per_volume);
raw_fit = as_overlap_fit(D, ["c2" "c1" D.names(D.beam_columns)]);
weights = 1 ./ (D.var_est + raw_fit.sigma_mismatch ^ 2);   % the same for every method
method_names = string(fieldnames(y_methods))';
overlap_scores = table();
for name = method_names
    s = as_overlap_score(D, y_methods.(name), common, weights);
    overlap_scores = [overlap_scores; table(name, s.disagreement, s.c2, s.c2_se, s.c1, s.c1_se, ...
        VariableNames=["method" "disagreement" "c2_left" "c2_se" "c1_left" "c1_se"])]; %#ok<AGROW>
end
disp(overlap_scores)
writetable(overlap_scores, fullfile(figure_dir, "AS_overlap_scores.csv"))

% Aspect function left after the recommended correction
D_corrected = D;
D_corrected.y = y_methods.height_ramp;
curve_corrected = as_aspect_curve(D_corrected, isfinite(D_corrected.y), 24);

%% Score the corrections on the single-radar harmonics at 10-60 km
versions = {"raw", "ell", "fix"};
for name = model_names
    versions{end + 1} = @(volume) as_sector_corrected(volume, models.(name)); %#ok<SAGROW>
end
kdmx = find(manifest.radars == "KDMX");
single_scores = table();
for i = 1:numel(versions)
    H = as_volume_harmonics(volumes, versions{i}, [2 3], "local");
    ok = isfinite(H.c2);
    kdmx_night = H.radar == kdmx & ismember(H.frame, 11:17);   % 04:30-07:30 UTC
    single_scores = [single_scores; table(method_names(i), median(H.first(ok)), mean(H.c1(ok)), ...
        mean(H.c2(ok)), mean(H.s2(ok)), median(H.wedge(ok)), mean(H.wedge(ok) < 0.8), min(H.wedge(kdmx_night)), ...
        VariableNames=["method" "first_median" "c1_mean" "c2_mean" "s2_mean" "wedge_median" ...
        "wedge_share_below_08" "kdmx_deepest"])]; %#ok<AGROW>
    if i == 1
        harmonics_raw = H;
    end
end
disp(single_scores)
writetable(single_scores, fullfile(figure_dir, "AS_single_radar_scores.csv"))

%% Per-radar amplitudes (one night): reproducible, but not the aspect alone
radar_folds = {as_fit_radar_offsets(D, even, radar_height, fold_models(1).height), ...
    as_fit_radar_offsets(D, ~even, radar_height, fold_models(2).height)};
both = isfinite(radar_folds{1}.estimate) & isfinite(radar_folds{2}.estimate);
radar_offsets = struct("fold_correlation", corr(radar_folds{1}.estimate(both), radar_folds{2}.estimate(both)), ...
    "sd", std(radar_folds{1}.estimate(both)), "median_se", median(radar_folds{1}.se(both)), ...
    "tau", radar_folds{1}.tau, "estimate", [radar_folds{1}.estimate radar_folds{2}.estimate]);
fprintf("per-radar b2 offsets: sd %.3f (se %.3f), fold correlation %.2f\n", radar_offsets.sd, ...
    radar_offsets.median_se, radar_offsets.fold_correlation)

%% Save
volume_harmonics = harmonics_raw;
save(fullfile(data_root, "aspect", "overlap_results.mat"), "harmonics", "axis_test", "curve_raw", ...
    "curve_corrected", "single_by_height", "overlap_by_height", "fold_models", "models", "overlap_scores", ...
    "single_scores", "radar_offsets", "volume_harmonics")
