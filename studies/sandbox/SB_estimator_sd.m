%% Sandbox: is the estimator's own uncertainty large enough to matter?
% The gate estimator turns the gates of a 2 km window into a likelihood of
% the window's density (studies/gate_density report, Eq. 5). Its spread,
% sd (birds km^-3, estimate_sweep), is the measurement error of each
% observation. The inversion leaves its own misfit: the residual
% y - F rho of each observation, which holds the measurement error and
% everything the reconstruction cannot represent (variation within a grid
% cell, the forward model, smoothing).
% If the residual is far larger than sd almost everywhere, weighting by sd
% changes little; if they are comparable, propagating sd is worth it.
% Steps:
% 1. Clean the sandbox volumes again with the production settings,
%    keeping sd (clean_volume), into their own folder.
% 2. Build the observation rows a second time with sd in place of the
%    density (same gates, same rows), so each row gets its sd.
% 3. Solve the sandbox (unweighted) and compare, by expected birds per
%    gate lambda = prediction x counting volume: the RMS residual, the RMS
%    sd, and the share of the residual variance that sd explains.
% 4. Solve unweighted and with weights 1 / (sd^2 + 1.5^2), sd the spread
%    the estimator has at a predicted density (from the prior: one solve;
%    from the unweighted solution: two), scaled to mean variance.
% 5. Animate the three solutions: map and three cross-sections per time.
% Outputs: figures/estimator_sd.png, estimator_sd.csv,
% estimator_sd_solves.csv, estimator_sd_frames.mp4/.gif/.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
S = sandbox_load("sandbox");
c = case_settings("sandbox");
manifest = load(fullfile(c.data_dir, "manifest.mat"));
sd_dir = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes", "sandbox_sd");
density_dir = fullfile(sd_dir, "density");
spread_dir = fullfile(sd_dir, "sd");

%% 1. Clean again, keeping sd ---------------------------------------------------------
[~, stem] = fileparts(manifest.volume_files);
has = manifest.volume_files ~= "";
density_files = strings(size(manifest.volume_files));
density_files(has) = fullfile(density_dir, stem(has) + ".mat");
spread_files = strings(size(manifest.volume_files));
spread_files(has) = fullfile(spread_dir, stem(has) + ".mat");
if ~isfolder(density_dir), mkdir(density_dir); end
if ~isfolder(spread_dir), mkdir(spread_dir); end
todo = find(has & ~isfile(spread_files));
raw_files = manifest.raw_files(todo);
out_files = density_files(todo);
read_options = c.read_options; vp = c.vp; filter = c.filter;
status = strings(numel(todo), 1);
if ~isempty(todo)
    parfor k = 1:numel(todo)
        status(k) = bmm4d.radar.clean_volume(raw_files(k), read_options, vp, filter, out_files(k));
    end
end
disp(unique(status))
% Copy with sd as the data (the rows of step 2), and check the density
% against the production cache
max_difference = 0;
for k = todo(:)'
    volume = load(density_files(k)).volume;
    cached = load(manifest.volume_files(k)).volume;
    for i_e = 1:numel(volume.dens.sweeps)
        max_difference = max(max_difference, max(abs(volume.dens.sweeps(i_e).data - cached.dens.sweeps(i_e).data), [], "all", "omitnan"));
        volume.dens.sweeps(i_e).data = volume.dens.sweeps(i_e).sd;
    end
    volume.dens.sweeps = rmfield(volume.dens.sweeps, "sd");
    save(spread_files(k), "volume")
end
fprintf("Largest difference to the cached density: %.3g birds km^-3\n", max_difference)

%% 2. Observation rows of sd ---------------------------------------------------------------
c_sd = c;
c_sd.name = "sandbox_sd";
c_sd.data_dir = fullfile(bmm4d.util.local_data_dir(), "prod", c_sd.name);
if ~isfile(fullfile(c_sd.data_dir, "case.mat"))
    if ~isfolder(c_sd.data_dir), mkdir(c_sd.data_dir); end
    volume_files = spread_files; raw_files = manifest.raw_files; radars = manifest.radars; times = manifest.times;
    save(fullfile(c_sd.data_dir, "manifest.mat"), "raw_files", "radars", "times", "volume_files")
    copyfile(fullfile(c.data_dir, "static_clutter.mat"), c_sd.data_dir)
    bmm4d.build_observations(c_sd);
end
loaded = load(fullfile(c_sd.data_dir, "case.mat"), "observations");
sd_values = loaded.observations.values;
assert(isequal(size(sd_values), size(S.observations.values)), "The sd rows differ from the sandbox rows.")
% A few windows have an sd but no density row (46 of 8.5e6 on 21 Oct): the density's rows are used
fprintf("Observations: %d with a density, %d of them without sd\n", nnz(S.observations.values), ...
    nnz(S.observations.values ~= 0 & sd_values == 0))

%% 3. Residual of the solution against sd ------------------------------------------------
result = sandbox_run(S, struct("name", "estimator sd check (unweighted)", "video", false));
observed_mask = S.observations.values ~= 0;
observed = full(S.observations.values(observed_mask));
predicted_all = S.observations.F(:, S.active_id) * double(result.field);
predicted = predicted_all(observed_mask);
sd = full(sd_values(observed_mask));
range_m = repmat(double(S.observations.row_range_m), 1, size(observed_mask, 2));
range_m = range_m(observed_mask);
beam = bmm4d.radar.gate.beam();
lambda = max(predicted, 0) .* beam.volume_sample_km3(range_m);
residual = observed - predicted;
lambda_edges = 10 .^ (-2:0.5:3);
lambda_bin = discretize(lambda, lambda_edges);
ok = ~isnan(lambda_bin);
n_bin = numel(lambda_edges) - 1;
summary = table(sqrt(lambda_edges(1:end - 1) .* lambda_edges(2:end))', accumarray(lambda_bin(ok), 1, [n_bin 1]), ...
    accumarray(lambda_bin(ok), observed(ok), [n_bin 1], @median, NaN), ...
    accumarray(lambda_bin(ok), residual(ok), [n_bin 1], @(x) sqrt(mean(x .^ 2)), NaN), ...
    accumarray(lambda_bin(ok), sd(ok), [n_bin 1], @(x) sqrt(mean(x .^ 2)), NaN), ...
    VariableNames=["lambda" "n" "median_observed" "rms_residual" "rms_sd"]);
summary.sd_share_of_variance = summary.rms_sd .^ 2 ./ summary.rms_residual .^ 2;
summary.share_of_rows = summary.n / sum(summary.n);
% Same for the relative size, sd / observed and |residual| / observed (medians)
summary.median_relative_sd = accumarray(lambda_bin(ok), sd(ok) ./ observed(ok), [n_bin 1], @median, NaN);
summary.median_relative_residual = accumarray(lambda_bin(ok), abs(residual(ok)) ./ observed(ok), [n_bin 1], @median, NaN);
disp(summary)
writetable(summary, fullfile(figure_dir, "estimator_sd.csv"))
fprintf("All rows: RMS residual %.1f, RMS sd %.1f birds km^-3, sd share of the residual variance %.2f\n", ...
    sqrt(mean(residual .^ 2)), sqrt(mean(sd .^ 2)), mean(sd .^ 2) / mean(residual .^ 2))

%% Plot results -----------------------------------------------------------
% Left: |residual| and sd of each observation (a sample) against lambda,
% with the RMS per half decade. Right: share of the residual variance
% that sd explains, raw and x2 (the likelihood is too narrow, as radials
% and windows overlap), with the share of rows per bin.
fig = bmm4d.plot.create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Sandbox night (21 Oct 2026): estimator sd against the misfit of the solution, per observation", FontSize=16)
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
sample = find(ok);
sample = sample(randperm(numel(sample), min(20000, numel(sample))));
scatter(ax, lambda(sample), abs(residual(sample)), 4, [0.85 0.33 0.1], "filled", MarkerFaceAlpha=0.15, HandleVisibility="off")
scatter(ax, lambda(sample), sd(sample), 4, [0 0.45 0.74], "filled", MarkerFaceAlpha=0.15, HandleVisibility="off")
plot(ax, summary.lambda, summary.rms_residual, "-o", Color=[0.85 0.33 0.1], LineWidth=2.5, MarkerFaceColor=[0.85 0.33 0.1], ...
    DisplayName="residual y - F\rho (RMS)")
plot(ax, summary.lambda, summary.rms_sd, "-s", Color=[0 0.45 0.74], LineWidth=2.5, MarkerFaceColor=[0 0.45 0.74], ...
    DisplayName="estimator sd (RMS)")
set(ax, XScale="log", YScale="log", FontSize=13); xlim(ax, [0.01 1000]); ylim(ax, [0.01 3000])
xlabel(ax, "Expected birds per gate \lambda"); ylabel(ax, "Birds km^{-3}")
legend(ax, Location="northwest", FontSize=12)
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
yyaxis(ax, "right")
bar(ax, log10(summary.lambda), summary.share_of_rows, 0.9, FaceColor=[0.75 0.75 0.75], EdgeColor="none", DisplayName="share of observations")
ylabel(ax, "Share of observations"); ax.YAxis(2).Color = [0.4 0.4 0.4]; ylim(ax, [0 1])
yyaxis(ax, "left")
plot(ax, log10(summary.lambda), summary.sd_share_of_variance, "-s", Color=[0 0.45 0.74], LineWidth=2.5, ...
    MarkerFaceColor=[0 0.45 0.74], DisplayName="sd^2 / residual^2")
plot(ax, log10(summary.lambda), min(4 * summary.sd_share_of_variance, 1.5), "--s", Color=[0 0.45 0.74], LineWidth=2, ...
    DisplayName="(2 sd)^2 / residual^2")
yline(ax, 1, "k:", HandleVisibility="off")
ylabel(ax, "Share of the residual variance"); ylim(ax, [0 1.5]); ax.YAxis(1).Color = [0 0.45 0.74];
set(ax, FontSize=13); xlim(ax, [-2 3]); xticks(ax, -2:3); xticklabels(ax, ["0.01" "0.1" "1" "10" "100" "1000"])
xlabel(ax, "Expected birds per gate \lambda")
legend(ax, Location="northeast", FontSize=12)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "estimator_sd.png"), Resolution=150)
close(fig)

%% 4. Solves weighted by sd ----------------------------------------------------------------
% sigma of each observation: the spread the estimator has, at the
% window's range, for a predicted density (not for the window's own
% reading, which would give low readings more weight because they are
% low). From all observations of the night: median sd / density by
% expected birds per gate lambda and range, read at lambda = prediction x
% counting volume and times the prediction. The prediction is F x the
% radar prior (one solve), or F x the unweighted solution (two solves).
% An absolute floor, sigma^2 = sd^2 + floor^2, stands for what the
% estimator does not see (forward model, density varying within a cell);
% it also caps the largest weights, which set the number of PCG
% iterations (no floor: weights up to 230 and 187 iterations).
% Scaling: the prior, smoothness and advection strengths were tuned with
% every weight 1. Scaled to mean weight 1, a weight-1 observation stands
% for sigma = 1 / sqrt(mean(1 / sigma^2)), set by the most precise
% (near-empty) windows: 0.5 birds km^-3 on 21 Oct, so typical observations
% lose weight against the regularization and the field is smoothed.
% Scaled to mean variance, w = mean(sigma^2) / sigma^2, a weight-1
% observation stands for the RMS sigma.
sd_full = full(sd_values);
sd_full(~observed_mask) = 0;
observed_full = full(S.observations.values);
row_range = repmat(double(S.observations.row_range_m), 1, size(observed_mask, 2));
volume_sample = beam.volume_sample_km3(row_range(observed_mask));
% Relative spread of the estimator by lambda and range
lambda_grid_edges = 10 .^ (-2:0.25:3.5);
range_grid_edges_km = [0 5 10 15 20 30 40 60 80 120 160 250];
relative_sd = lookup_relative_sd(observed_full(observed_mask), sd_full(observed_mask), volume_sample, ...
    row_range(observed_mask) / 1000, lambda_grid_edges, range_grid_edges_km);
prior_keys = S.prior_cache.keys;
variants = table(["unweighted" "sd + floor 1.5, prior" "sd + floor 1.5, solution"]', [NaN 1.5 1.5]', ...
    ["" "prior" "solution"]', VariableNames=["name" "floor" "prediction"]);
solves = table();
solved_fields = cell(height(variants), 1);
for i_v = 1:height(variants)
    option = struct("name", "estimator sd weights: " + variants.name(i_v), "video", false);
    if i_v > 1
        if variants.prediction(i_v) == "prior"
            predicted_field = S.prior_cache(prior_keys{1});
        else
            predicted_field = solved_fields{1};
        end
        prediction = max(S.observations.F(:, S.active_id) * double(predicted_field), 0.1);
        prediction = prediction(observed_mask);
        sigma2 = (prediction .* read_relative_sd(relative_sd, prediction .* volume_sample, row_range(observed_mask) / 1000, ...
            lambda_grid_edges, range_grid_edges_km)) .^ 2 + variants.floor(i_v) ^ 2;
        W = zeros(size(observed_mask));
        W(observed_mask) = mean(sigma2) ./ sigma2;
        option.row_weight = @(range_m, elev_deg, is_zero) W;
    end
    r = sandbox_run(S, option);
    solved_fields{i_v} = r.field;
    solves = [solves; table(variants.name(i_v), r.bias_200km, r.misfit, r.misfit_trusted, r.detail, ...
        sum(double(r.field), "all") / sum(double(solved_fields{1}), "all"), r.seconds, r.iterations, ...
        VariableNames=["weights" "bias_200km" "misfit" "misfit_trusted" "detail" "total_birds" "seconds" "iterations"])]; %#ok<AGROW>
end
disp(solves)
writetable(solves, fullfile(figure_dir, "estimator_sd_solves.csv"))

%% 5. Animation: map and cross-sections of each solve ------------------------------------
% Rows: the three solves. Columns: column density from above (section
% lines in white), west-east sections through KDMX and KLSX, and a
% south-north section through KDVN. One frame per time (mp4 and gif), the
% middle frame also as PNG.
[ix, iy, ~] = ind2sub(S.g.sz, S.active_id);
column_id = sub2ind(S.g.sz(1:2), ix, iy);
air = bmm4d.grid.air_fraction(S.g, S.active_id);
to_column = @(x) accumarray(column_id, double(x) .* air, [prod(S.g.sz(1:2)) 1]) * S.g.dz / 1000;   % birds km^-2
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~S.g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* S.g.mask_water, [2 1 3]);
height_m = squeeze(S.g.z);
sections = table(["KDMX" "KLSX" "KDVN"]', ["we" "we" "sn"]', VariableNames=["radar" "direction"]);
for i_k = 1:height(sections)
    k = find(S.radar_selected.name == sections.radar(i_k), 1);
    sections.lon(i_k) = S.radar_selected.lon(k);
    sections.lat(i_k) = S.radar_selected.lat(k);
    [~, sections.x_id(i_k)] = min(abs(S.g.lon - sections.lon(i_k)));
    [~, sections.y_id(i_k)] = min(abs(S.g.lat - sections.lat(i_k)));
end
active_column = to_column(ones(numel(S.active_id), 1)) > 0;
n_time = size(solved_fields{1}, 2);
file_stem = fullfile(figure_dir, "estimator_sd_frames");
video = VideoWriter(file_stem + ".mp4", "MPEG-4");
video.FrameRate = 1;
open(video)
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, height(variants), 1 + height(sections), TileSpacing="tight", Padding="tight");
    title(t, "Sandbox night, " + string(S.c.times(i_t), "dd MMM yyyy HH:mm") + " UTC", FontSize=18, FontWeight="bold")
    for i_v = 1:height(variants)
        field = nan(S.g.sz, "single");
        field(S.active_id) = solved_fields{i_v}(:, i_t);
        % From above
        ax = nexttile(t);
        column = nan(S.g.sz(1:2));
        values = to_column(solved_fields{i_v}(:, i_t));
        column(active_column) = log10(max(values(active_column), 1));
        bmm4d.plot.plot_night_map(ax, S.g, background, column, log10([5 300]), states, S.radar_selected, S.c.figure_bounds, ones(S.g.sz(1:2)));
        colormap(ax, parula)
        hold(ax, "on")
        for i_k = 1:height(sections)
            if sections.direction(i_k) == "we"
                plot(ax, sections.lon(i_k) + [-3 3], S.g.lat(sections.y_id(i_k)) * [1 1], "w-", LineWidth=2)
            else
                plot(ax, S.g.lon(sections.x_id(i_k)) * [1 1], sections.lat(i_k) + [-3 3], "w-", LineWidth=2)
            end
        end
        ylabel(ax, variants.name(i_v), FontWeight="bold", FontSize=15)
        if i_v == 1, title(ax, "Column, 5-300 birds km^{-2} (log)", FontSize=14), end
        % Sections
        for i_k = 1:height(sections)
            ax = nexttile(t);
            if sections.direction(i_k) == "we"
                values = squeeze(field(:, sections.y_id(i_k), :))';
                imagesc(ax, S.g.lon, height_m, values, AlphaData=~isnan(values))
                xlim(ax, sections.lon(i_k) + [-3 3]); hold(ax, "on"); xline(ax, sections.lon(i_k), "k:", LineWidth=1.5)
                header = sprintf("W-E through %s (%.2f°N)", sections.radar(i_k), S.g.lat(sections.y_id(i_k)));
            else
                values = squeeze(field(sections.x_id(i_k), :, :))';
                imagesc(ax, S.g.lat, height_m, values, AlphaData=~isnan(values))
                xlim(ax, sections.lat(i_k) + [-3 3]); hold(ax, "on"); xline(ax, sections.lat(i_k), "k:", LineWidth=1.5)
                header = sprintf("S-N through %s (%.2f°W)", sections.radar(i_k), -S.g.lon(sections.x_id(i_k)));
            end
            set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12); ylim(ax, [0 3000]); clim(ax, [0 150]); colormap(ax, parula)
            if i_v == 1, title(ax, header, FontSize=14), end
            if i_k == 1, ylabel(ax, "Height (m a.s.l.)"), end
            if i_v == height(variants), xlabel(ax, ternary_label(sections.direction(i_k))), end
        end
    end
    cb = colorbar(ax); cb.Layout.Tile = "east"; cb.Label.String = "Density in the sections (birds km^{-3})"; cb.FontSize = 13;
    drawnow
    frame_image = getframe(fig);
    writeVideo(video, frame_image)
    [indexed, color_map] = rgb2ind(frame_image.cdata, 256);
    if i_t == 1
        imwrite(indexed, color_map, file_stem + ".gif", LoopCount=Inf, DelayTime=1.2)
    else
        imwrite(indexed, color_map, file_stem + ".gif", WriteMode="append", DelayTime=1.2)
    end
    if i_t == ceil(n_time / 2)
        bmm4d.plot.export_fixed_figure(fig, file_stem + ".png", Resolution=150)
    end
end
close(video)
close(fig)

%% Local functions ----------------------------------------------------------
function relative_sd = lookup_relative_sd(observed, sd, volume_sample, range_km, lambda_edges, range_edges_km)
% Median sd / density of the windows by expected birds per gate (from the
% window's own density) and range; empty cells take the nearest lambda
lambda_bin = discretize(observed .* volume_sample, lambda_edges);
range_bin = discretize(range_km, range_edges_km);
ok = ~isnan(lambda_bin) & ~isnan(range_bin) & observed > 0;
relative_sd = accumarray([lambda_bin(ok) range_bin(ok)], sd(ok) ./ observed(ok), ...
    [numel(lambda_edges) - 1, numel(range_edges_km) - 1], @median, NaN);
relative_sd = fillmissing(relative_sd, "nearest", 1);
relative_sd = fillmissing(relative_sd, "nearest", 2);
end

function relative = read_relative_sd(relative_sd, lambda, range_km, lambda_edges, range_edges_km)
% Relative sd at a given lambda and range (cells of the lookup, ends clamped)
lambda_bin = min(max(discretize(min(max(lambda, lambda_edges(1)), lambda_edges(end)), lambda_edges), 1), numel(lambda_edges) - 1);
range_bin = min(max(discretize(min(range_km, range_edges_km(end)), range_edges_km), 1), numel(range_edges_km) - 1);
relative = relative_sd(sub2ind(size(relative_sd), lambda_bin, range_bin));
end

function label = ternary_label(direction)
% Axis label of a section
if direction == "we"
    label = "Longitude (°)";
else
    label = "Latitude (°)";
end
end
