function result = two_level_run(S, option)
%TWO_LEVEL_RUN sandbox_run with the separable and two-level solvers (archived).
%   OPTION is a struct with a name and any of these fields (defaults from
%   prod_case): prior, smoothness, advection, coeff_horver, shape,
%   prior_correlation_m, prior_correlation_vertical_m, prior_mean
%   ("radar": the radar-profile prior, default; "zero"; "profile": one
%   domain-mean vertical profile per time, no horizontal structure),
%   lambda_profile, prior_interpolation,
%   prior_length_scale_m, and row_weight, a function handle
%   w = f(range_m, elev_deg, is_zero) giving each observation row a weight.
%   SOLVER "separable" uses solve_separable (column amounts x pattern
%   shapes above the smoothed ground, shapes from shape_patterns.mat of the
%   us_east case) with its options n_pattern, shape_prior, shape_length_m,
%   column_prior, level_floor, negative_penalty, ground_smoothing_m, and
%   regional_patterns J (> 0: the mean shape is bent per column by the
%   radar profiles' coefficients on patterns 2..J+1, interpolated over
%   regional_length_m, profiles fitted with regional_lambda_profile, default 1;
%   regional_shape.m); the level and the advection
%   relaxation come from the radar prior. SOLVER "two_level" uses
%   solve_two_level (coarse level V_L and coarse shape coefficients, fine
%   deviation); its options go in the struct field two_level, and
%   regional_patterns J > 0 starts its shape from the radar profiles.
%   RESULT has the field and metrics: radar-distance bias (all bins and
%   within 200 km), data misfit on all original rows and on trusted rows
%   (30-120 km from the radar), detail (RMS of the column relative to its
%   regional mean, i.e. small-scale structure kept), time and iterations.
%   A panel (map, KDMX section, bias curve) goes to studies/two_level/figures/sandbox
%   unless option.panel is false.

inv = S.c.inversion;
names = ["prior" "smoothness" "advection" "coeff_horver" "lambda_profile" "prior_interpolation" "prior_length_scale_m"];
for name = names
    if ~isfield(option, name)
        option.(name) = inv.(name);
    end
end
defaults = struct("shape", 0, "prior_correlation_m", 0, "prior_correlation_vertical_m", 0, ...
    "row_weight", [], "prior_mean", "radar", "solver", "lean", "n_pattern", 3, "shape_prior", 0.3, ...
    "shape_length_m", 100000, "column_prior", 0, "level_floor", 10, "negative_penalty", 3, ...
    "ground_smoothing_m", 30000, "panel", true, "regional_patterns", 0, "regional_length_m", 150000, "regional_lambda_profile", 1, ...
    "two_level", struct());
for name = string(fieldnames(defaults))'
    if ~isfield(option, name)
        option.(name) = defaults.(name);
    end
end
run_tic = tic;

% Density prior (cached by its settings)
key = sprintf("%g|%s|%g", option.lambda_profile, option.prior_interpolation, option.prior_length_scale_m);
if isKey(S.prior_cache, key)
    prior_active = S.prior_cache(key);
else
    evalc("prior_active = bmm4d.prior.build_density_prior_fast(S.g, S.radar_scans, S.observations, " + ...
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

if option.solver == "separable" || option.solver == "two_level"
    patterns = load(fullfile(prod_case("us_east").data_dir, "shape_patterns.mat"));
    % Regional mean shape from the radar profiles (patterns 2..J+1), or the night's mean shape
    shape_coefficients = [];
    if option.regional_patterns > 0
        % Radar profiles with weak vertical smoothing (the prior's lambda_profile 10 flattens them)
        shape_key = sprintf("%g|profiles only", option.regional_lambda_profile);
        if ~isKey(S.prior_cache, shape_key)
            evalc("[~, shape_info] = bmm4d.prior.build_density_prior_fast(S.g, S.radar_scans, S.observations, " + ...
                """lambda_profile"", option.regional_lambda_profile, ""lambda_time_profile"", inv.lambda_time_profile, " + ...
                """maxit"", inv.maxit, ""active_id"", 1, ""include_water"", true);");
            S.prior_cache(shape_key) = shape_info.radar_profiles;
        end
        shape_coefficients = regional_shape(S.g, S.radar_scans, S.prior_cache(shape_key), patterns.shape_patterns, ...
            patterns.z_agl, option.regional_patterns + 1, option.regional_length_m, option.ground_smoothing_m);
    end
end
if option.solver == "two_level"
    % Two-level solver: OPTION.two_level holds its own options (solve_two_level)
    arguments_two_level = [{"n_pattern", option.n_pattern, "shape_sd", patterns.shape_sd, ...
        "shape_initial", shape_coefficients, "ground_smoothing_m", option.ground_smoothing_m, ...
        "level_floor", option.level_floor, "advection", option.advection, ...
        "advection_relax_exponent", inv.advection_relax_exponent, "negative_penalty", option.negative_penalty, ...
        "row_weights", row_weights, "tol", inv.tol, "maxit", inv.maxit, "projected_maxit", inv.projected_maxit}, ...
        namedargs2cell(option.two_level)];
    [field, info] = solve_two_level(S.observations, S.g, S.active_mask_3d, S.steps, radar_prior, ...
        patterns.shape_patterns, patterns.z_agl, arguments_two_level{:});
elseif option.solver == "separable"
    [field, info] = solve_separable(S.observations, S.g, S.active_mask_3d, S.steps, radar_prior, ...
        patterns.shape_patterns, patterns.z_agl, n_pattern=option.n_pattern, shape_sd=patterns.shape_sd, ...
        regional_shape=shape_coefficients, ...
        smoothness=option.smoothness, shape_prior=option.shape_prior, shape_length_m=option.shape_length_m, ...
        column_prior=option.column_prior, level_floor=option.level_floor, ground_smoothing_m=option.ground_smoothing_m, ...
        advection=option.advection, advection_relax_exponent=inv.advection_relax_exponent, ...
        negative_penalty=option.negative_penalty, row_weights=row_weights, tol=inv.tol, maxit=inv.maxit, ...
        projected_maxit=inv.projected_maxit);
else
    [field, info] = bmm4d.solve.solve_advection_lean(S.observations, S.g, S.active_mask_3d, S.steps, prior_active, ...
        prior=option.prior, smoothness=option.smoothness, advection=option.advection, coeff_horver=option.coeff_horver, ...
        tol=inv.tol, maxit=inv.maxit, store_class="single", advection_relax_exponent=inv.advection_relax_exponent, ...
        nonnegative=true, projected_maxit=inv.projected_maxit, row_weights=row_weights, ...
        shape=option.shape, prior_correlation_m=option.prior_correlation_m, ...
        prior_correlation_vertical_m=option.prior_correlation_vertical_m, advection_reference=radar_prior);
end
seconds_used = toc(run_tic);

% Metrics: radar-distance bias and misfit on all original rows
[ratio, edges_km, bias_all] = bmm4d.plot.radar_distance_bias(field, S.active_id, S.g, S.radar_selected);
within = edges_km(2:end) <= 200;
ratio_within = ratio(within, :) ./ mean(ratio(within, :), 1, "omitnan");
bias_200 = sqrt(mean(log(ratio_within) .^ 2, 1, "omitnan"));
predicted = S.observations.F(:, S.active_id) * double(field);
observed = full(S.observations.values);
mask = observed ~= 0;
misfit = norm(predicted(mask) - observed(mask)) / norm(observed(mask));
trusted = mask & S.observations.row_range_m >= 30000 & S.observations.row_range_m <= 120000;
misfit_trusted = norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted));
detail = column_detail(S, field);
result = struct("name", string(option.name), "field", field, "ratio", ratio, "edges_km", edges_km, ...
    "bias_all", mean(bias_all, "omitnan"), "bias_200km", mean(bias_200, "omitnan"), "misfit", misfit, ...
    "misfit_trusted", misfit_trusted, "detail", detail, ...
    "seconds", seconds_used, "iterations", info.solver.iterations, "option", option, "solver", info.solver, "info", info);
if option.panel
    sandbox_panel(S, result);
end
end

function sandbox_panel(S, result)
% Map at the middle frame, KDMX west-east section, radar-distance curve.
g = S.g;
t_mid = ceil(size(result.field, 2) / 2);
field = nan(g.sz, "single");
field(S.active_id) = result.field(:, t_mid);
column = sum(field, 3, "omitnan") .* g.dz ./ 1000;
column(~any(S.active_mask_3d, 3)) = nan;
k = find(S.radar_selected.name == "KDMX", 1);
[~, section_y] = min(abs(g.lat - S.radar_selected.lat(k)));
fig = bmm4d.plot.create_export_figure(3600, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
ax = nexttile(t);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
bmm4d.plot.plot_night_map(ax, g, background, column, [0 200], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
yline(ax, g.lat(section_y), "w--")
colorbar(ax, Color="k");
title(ax, result.name + " | " + string(S.c.times(t_mid), "HH:mm") + " UTC", FontSize=13, Interpreter="none")
ax = nexttile(t);
section = squeeze(field(:, section_y, :))';
imagesc(ax, g.lon, squeeze(g.z), section, AlphaData=~isnan(section))
set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12); clim(ax, [0 150]); colorbar(ax)
xlim(ax, S.radar_selected.lon(k) + [-3 3]); ylim(ax, [0 3000])
title(ax, sprintf("Section through KDMX (%.2f°N)", g.lat(section_y)), FontSize=13)
ax = nexttile(t); hold(ax, "on"); box(ax, "on")
center = (result.edges_km(1:end-1) + result.edges_km(2:end)) / 2;
plot(ax, center, result.ratio, Color=[0.7 0.7 0.7])
plot(ax, center, mean(result.ratio, 2, "omitnan"), "k-o", LineWidth=2)
yline(ax, 1, "k--"); ylim(ax, [0.6 1.4]); set(ax, FontSize=12)
xlabel(ax, "Distance to nearest radar (km)"); ylabel(ax, "Column / regional mean")
title(ax, sprintf("bias %.3f, misfit %.3f (trusted %.3f), detail %.2f", result.bias_200km, result.misfit, result.misfit_trusted, result.detail), FontSize=12)
out_dir = fullfile(setup_paths(), "studies", "two_level", "figures", "sandbox");
bmm4d.plot.export_fixed_figure(fig, fullfile(out_dir, matlab.lang.makeValidName(result.name) + ".png"), Resolution=150)
close(fig)
end

function detail = column_detail(S, field)
% RMS of log(column / regional mean) over active columns with birds, mean
% over time: the small-scale structure the map keeps (0 = fully smooth).
g = S.g;
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
is_active = false(g.sz(1:2));
is_active(column_id) = true;
sigma_cells = 150000 / mean([g.dx g.dy]);
weight = imgaussfilt(double(is_active), sigma_cells);
values = nan(1, size(field, 2));
for i_t = 1:size(field, 2)
    column = reshape(accumarray(column_id, double(field(:, i_t)), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    column = max(column, 0) .* is_active;
    regional = imgaussfilt(column, sigma_cells) ./ weight;
    ok = is_active & regional > 5 & column > 0;
    values(i_t) = sqrt(mean(log(column(ok) ./ regional(ok)) .^ 2));
end
detail = mean(values, "omitnan");
end
