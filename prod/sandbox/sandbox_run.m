function result = sandbox_run(S, option)
%SANDBOX_RUN Solve the sandbox case with one option and score it.
%   OPTION is a struct with a name and any of these fields (defaults from
%   prod_case): prior, smoothness, advection, coeff_horver, shape,
%   prior_correlation_m, prior_correlation_vertical_m, prior_mean
%   ("radar": the radar-profile prior, default; "zero"; "profile": one
%   domain-mean vertical profile per time, no horizontal structure),
%   lambda_profile, prior_interpolation,
%   prior_length_scale_m, and row_weight, a function handle
%   w = f(range_m, elev_deg, is_zero) giving each observation row a weight.
%   NONNEGATIVE (default false) adds the projected FISTA step after PCG;
%   without it the solve is plain PCG (fast, negative cells allowed).
%   RESULT has the field and metrics: radar-distance bias (all bins and
%   within 200 km), data misfit on all original rows and on trusted rows
%   (30-120 km from the radar), detail (RMS of the column relative to its
%   regional mean, i.e. small-scale structure kept), negative share (negative
%   mass / total absolute mass), solve time and iterations. Each run appends
%   a row to prod/sandbox/figures/sandbox_log.csv and, unless option.video
%   is false, writes a video of all frames (column map, column / regional
%   mean, KDMX section, radar-distance curve) and its middle frame to
%   prod/sandbox/figures/runs/<name>.mp4 and .png.

inv = S.c.inversion;
names = ["prior" "smoothness" "advection" "coeff_horver" "lambda_profile" "prior_interpolation" "prior_length_scale_m"];
for name = names
    if ~isfield(option, name)
        option.(name) = inv.(name);
    end
end
defaults = struct("shape", 0, "prior_correlation_m", 0, "prior_correlation_vertical_m", 0, ...
    "row_weight", [], "prior_mean", "radar", "nonnegative", false, "video", true);
for name = string(fieldnames(defaults))'
    if ~isfield(option, name)
        option.(name) = defaults.(name);
    end
end
% Density prior (cached by its settings)
key = sprintf("%g|%s|%g", option.lambda_profile, option.prior_interpolation, option.prior_length_scale_m);
if isKey(S.prior_cache, key)
    prior_active = S.prior_cache(key);
else
    evalc("prior_active = build_density_prior_fast(S.g, S.radar_scans, S.observations, " + ...
        """lambda_profile"", option.lambda_profile, ""lambda_time_profile"", inv.lambda_time_profile, " + ...
        """maxit"", inv.maxit, ""active_id"", S.active_id, ""interpolation"", option.prior_interpolation, " + ...
        """length_scale_m"", option.prior_length_scale_m, ""include_water"", true);");
    S.prior_cache(key) = prior_active;   % containers.Map is a handle: cached for later calls
end

% Prior mean: the radar prior, zero, or its horizontal mean per layer
radar_prior = prior_active;
switch string(option.prior_mean)
    case "zero"
        prior_active = zeros(size(radar_prior), "like", radar_prior);
    case "profile"
        [~, ~, layer] = ind2sub(S.g.sz, S.active_id);
        layer_mean = zeros(S.g.sz(3), size(radar_prior, 2));
        for i_t = 1:size(radar_prior, 2)
            layer_mean(:, i_t) = accumarray(layer, double(radar_prior(:, i_t)), [S.g.sz(3) 1], @mean, 0);
        end
        prior_active = cast(layer_mean(layer, :), "like", radar_prior);
end

% Row weights
row_weights = [];
if ~isempty(option.row_weight)
    row_weights = option.row_weight(double(S.observations.row_range_m), double(S.observations.row_elev_deg), S.row_is_zero);
end

solve_tic = tic;
[field, info] = solve_advection_lean(S.observations, S.g, S.active_mask_3d, S.steps, prior_active, ...
    prior=option.prior, smoothness=option.smoothness, advection=option.advection, coeff_horver=option.coeff_horver, ...
    tol=inv.tol, maxit=inv.maxit, store_class="single", advection_relax_exponent=inv.advection_relax_exponent, ...
    nonnegative=option.nonnegative, projected_maxit=inv.projected_maxit, row_weights=row_weights, ...
    shape=option.shape, prior_correlation_m=option.prior_correlation_m, ...
    prior_correlation_vertical_m=option.prior_correlation_vertical_m, advection_reference=radar_prior);
seconds_used = toc(solve_tic);

% Metrics: radar-distance bias and misfit on all original rows
[ratio, edges_km, bias_all] = radar_distance_bias(field, S.active_id, S.g, S.radar_selected);
within = edges_km(2:end) <= 200;
ratio_within = ratio(within, :) ./ mean(ratio(within, :), 1, "omitnan");
bias_200 = sqrt(mean(log(ratio_within) .^ 2, 1, "omitnan"));
predicted = S.observations.F(:, S.active_id) * double(field);
observed = full(S.observations.values);
mask = observed ~= 0;
misfit = norm(predicted(mask) - observed(mask)) / norm(observed(mask));
trusted = mask & S.observations.row_range_m >= 30000 & S.observations.row_range_m <= 120000;
misfit_trusted = norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted));
[column, relative, is_active] = column_maps(S, field);
% Detail: RMS of log(column / regional mean) per frame, mean over frames
detail = mean(sqrt(mean(log(relative) .^ 2, [1 2], "omitnan")), "omitnan");
negative_share = sum(max(-double(field), 0), "all") / max(sum(abs(double(field)), "all"), eps);
result = struct("name", string(option.name), "field", field, "ratio", ratio, "edges_km", edges_km, ...
    "bias_all", mean(bias_all, "omitnan"), "bias_200km", mean(bias_200, "omitnan"), "bias_frame", bias_200, ...
    "misfit", misfit, "misfit_trusted", misfit_trusted, "detail", detail, "negative_share", negative_share, ...
    "seconds", seconds_used, "iterations", info.solver.iterations, "option", option, "solver", info.solver);

% Log every run, so progress can be followed across sessions
figure_dir = fullfile(setup_paths(), "prod", "sandbox", "figures");
log_row = table(datetime("now"), result.name, option.nonnegative, result.bias_200km, result.misfit, ...
    result.misfit_trusted, result.detail, negative_share, seconds_used, result.iterations, ...
    VariableNames=["time" "name" "nonnegative" "bias_200km" "misfit" "misfit_trusted" "detail" ...
    "negative_share" "seconds" "iterations"]);
log_file = fullfile(figure_dir, "sandbox_log.csv");
writetable(log_row, log_file, WriteMode="append", WriteVariableNames=~isfile(log_file))
if option.video
    sandbox_video(S, result, column, relative, is_active, fullfile(figure_dir, "runs"));
end
end

function sandbox_video(S, result, column, relative, is_active, out_dir)
% One frame per time: column map, column / regional mean (the radar
% imprint), west-east section through KDMX, radar-distance curve of the
% frame. The middle frame is also saved as a PNG.
g = S.g;
n_time = size(result.field, 2);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
k = find(S.radar_selected.name == "KDMX", 1);
[~, section_y] = min(abs(g.lat - S.radar_selected.lat(k)));
center = (result.edges_km(1:end-1) + result.edges_km(2:end)) / 2;
column(~is_active) = nan;
if ~isfolder(out_dir)
    mkdir(out_dir);
end
file_stem = fullfile(out_dir, matlab.lang.makeValidName(result.name));
video = VideoWriter(file_stem + ".mp4", "MPEG-4");
video.FrameRate = 2;
open(video)
fig = create_export_figure(2400, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, 2, 2, TileSpacing="tight", Padding="tight");
    title(t, [sprintf("%s | %s UTC | bias %.3f", result.name, string(S.c.times(i_t), "HH:mm"), result.bias_frame(i_t)), ...
        sprintf("Run: bias %.3f, trusted misfit %.3f, detail %.2f, negative %.1f%%, %.0f s", result.bias_200km, ...
        result.misfit_trusted, result.detail, 100 * result.negative_share, result.seconds)], ...
        FontSize=13, Color="k", Interpreter="none")
    ax = nexttile(t);
    plot_night_map(ax, g, background, column(:, :, i_t), [0 200], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    yline(ax, g.lat(section_y), "w--")
    cb = colorbar(ax, Color="k"); cb.Label.String = "Column density (birds km^{-2})";
    ax = nexttile(t);
    plot_night_map(ax, g, background, log2(relative(:, :, i_t)), [-1 1], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, signed_map); cb = colorbar(ax, Color="k"); cb.Label.String = "log_2(column / 150 km regional mean)";
    ax = nexttile(t);
    field = nan(g.sz, "single");
    field(S.active_id) = result.field(:, i_t);
    section = squeeze(field(:, section_y, :))';
    imagesc(ax, g.lon, squeeze(g.z), section, AlphaData=~isnan(section))
    set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12); clim(ax, [0 150])
    cb = colorbar(ax, Color="k"); cb.Label.String = "Density (birds km^{-3})";
    xlim(ax, S.radar_selected.lon(k) + [-3 3]); ylim(ax, [0 3000])
    title(ax, sprintf("Section through KDMX (%.2f°N)", g.lat(section_y)), FontSize=12, FontWeight="normal")
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, center, result.ratio, Color=[0.75 0.75 0.75])
    plot(ax, center, result.ratio(:, i_t), "k-o", LineWidth=2)
    yline(ax, 1, "k--"); ylim(ax, [0.6 1.4]); set(ax, FontSize=12)
    xlabel(ax, "Distance to nearest radar (km)"); ylabel(ax, "Column / regional mean")
    title(ax, "Radar-distance curve (grey: other frames)", FontSize=12, FontWeight="normal")
    drawnow
    writeVideo(video, getframe(fig))
    if i_t == ceil(n_time / 2)
        export_fixed_figure(fig, file_stem + ".png", Resolution=150)
    end
end
close(video)
close(fig)
end

function [column, relative, is_active] = column_maps(S, field)
% Column density [x y t] (birds km^-2, negatives clipped) and its ratio to
% the 150 km regional mean (NaN where the regional mean is below 5 birds
% km^-2 or outside the active columns).
g = S.g;
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
is_active = false(g.sz(1:2));
is_active(column_id) = true;
sigma_cells = 150000 / mean([g.dx g.dy]);
weight = imgaussfilt(double(is_active), sigma_cells);
n_time = size(field, 2);
column = zeros([g.sz(1:2) n_time]);
relative = nan([g.sz(1:2) n_time]);
for i_t = 1:n_time
    map = reshape(accumarray(column_id, double(field(:, i_t)), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    map = max(map, 0) .* is_active;
    regional = imgaussfilt(map, sigma_cells) ./ weight;
    ok = is_active & regional > 5 & map > 0;
    ratio_map = nan(g.sz(1:2));
    ratio_map(ok) = map(ok) ./ regional(ok);
    column(:, :, i_t) = map;
    relative(:, :, i_t) = ratio_map;
end
end
