%% Gate density M: a level per cut type and distance, bounded by the clutter filter
% Dev study. Within one volume the same height is seen by several cuts at
% different ranges. Birds vary smoothly over a few km; what the clutter
% filter removes belongs to one cut's geometry. So a cut that reads low only
% where its filter removed power, against a cut at the same height that was
% not filtered, lost birds there; where cuts agree after filtering, the
% removed power was clutter. A single sweep cannot tell.
% Per scan, one vertical profile p (100 m layers, lightly smoothed) is fitted to
% the estimator's mean density per sweep and range bin (GD_L, production
% estimator) through the beam's vertical weighting:
%   density = c(type, band) * beta(lambda) * sum_z W(z) p(z)
% - c, one level per cut type (split: surveillance pass of a cut scanned twice;
%   batch: one pass) and distance band (2-10, 10-20, 20-30, 30-45, 45-80 km),
%   kept within [P_out / P_all, 1.1] (GD_CFP_POWER: the share of the echo the
%   clutter filter left; the removed power is birds or clutter, so the birds
%   kept are at least that share). Split cuts at 45-80 km are the anchor (1).
%   In the refined fit each c is also pulled toward 1 by one sweep's worth of
%   pseudo-observations that read as predicted: a level backed by many sweeps
%   moves freely, one backed by a single sweep only halfway.
% - beta(lambda), the sampling deficit the estimator leaves where gates hold
%   few birds, one curve for all scans and cuts, by expected birds per gate
%   lambda = (W p) x counting volume (bmm4d.radar.gate.beam); 1 from 10 birds
%   per gate. Fitted in rounds with the scans.
% c is linear in distance between the band centres and beta linear in log lambda
% between its knots, so that neither steps at a band or bin edge.
% Three fits: before (p only), levels (p and bounded c, as on 4 Nov), refined
% (p, bounded and pulled c, beta). Elevation agreement: gate-weighted RMS and
% mean of log((observed + 1) / (fitted + 1)) per sweep and range bin within 45
% km (1 bird km^-3 added so that near-empty bins do not dominate).
% Outputs: studies/gate_density/figures/GD_M_cut_levels.csv, GD_M_summary.png,
% GD_M_examples_1-4.png and GD_M_cut_levels.mat.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");
L = load(fullfile(figure_dir, "GD_L_cfp_bounds.mat"));
scans = L.scans;
n_scan = height(scans);
beam = bmm4d.radar.gate.beam();

%% Settings ---------------------------------------------------------------
band_edges_km = [2 10 20 30 45 80];
knots_km = (band_edges_km(1:end - 1) + band_edges_km(2:end)) / 2;   % levels are linear between band centres
n_band = numel(band_edges_km) - 1;
group_name = ["split" "batch"];
upper_level = 1.1;
agreement_km = 45;
estimator_method = find(L.method_name == "estimator");
lambda_knots = [0.01 0.03 0.1 0.3 1 3 10];         % expected birds per gate; beta = 1 from the last knot
n_round = 4;                                        % rounds of beta with the scans
fit_name = ["before" "levels" "refined"];
n_fit = numel(fit_name);

%% 1. Reported and removed power per scan, cut and range bin ------------------------
power_file = fullfile(figure_dir, "GD_M_cfp_power.mat");
if isfile(power_file)
    load(power_file, "power")
else
    power = cell(n_scan, 1);
    pool = gcp("nocreate");
    if isempty(pool) || pool.NumWorkers > 4
        delete(pool);
        parpool("Processes", 4);
    end
    raw_files = scans.raw_file;
    bin_edges_km = L.bin_edges_km;
    parfor i_s = 1:n_scan
        power{i_s} = gd_cfp_power(raw_files(i_s), bin_edges_km, 10);
    end
    save(power_file, "power")
end

%% 2. Observations and bounds per scan -------------------------------------------------
n_bin = numel(L.bin_km);
band_of_bin = discretize(L.bin_km, band_edges_km);
obs = cell(n_scan, 1);
lower_bound = nan(n_band, 2, n_scan);
for i_s = 1:n_scan
    elevation = L.sweep_elevation(:, i_s);
    n_sweep = nnz(~isnan(elevation));
    % The power table has the same cut selection as GD_L
    assert(numel(power{i_s}.elev) >= n_sweep && all(abs(power{i_s}.elev(1:n_sweep) - elevation(1:n_sweep)) < 0.05), ...
        "Sweeps of scan %d do not match", i_s)
    is_split = power{i_s}.is_split(1:n_sweep);
    [i_bin, i_sweep] = ndgrid(1:n_bin, 1:n_sweep);
    density = reshape(L.density_sum(:, 1:n_sweep, i_s, estimator_method) ./ L.gate_count(:, 1:n_sweep, i_s, estimator_method), [], 1);
    count = reshape(L.gate_count(:, 1:n_sweep, i_s, estimator_method), [], 1);
    keep = count > 0 & ~isnan(density);
    obs{i_s} = struct("density", density(keep), "count", count(keep), "range_m", L.bin_km(i_bin(keep))' * 1000, ...
        "range_lo", L.bin_edges_km(i_bin(keep))' * 1000, "range_hi", L.bin_edges_km(i_bin(keep) + 1)' * 1000, ...
        "elevation", elevation(i_sweep(keep)), "sweep", i_sweep(keep), "band", band_of_bin(i_bin(keep))', ...
        "group", 2 - is_split(i_sweep(keep)), "antenna", scans.antenna_agl(i_s), "n_sweep", n_sweep, "is_split", is_split);
    % Bounds: share of the echo the filter left, per type and band (no data: fixed at 1)
    bounds = ones(n_band, 2, 2);
    for g = 1:2
        of_type = is_split == (g == 1);
        for b = 1:n_band
            p_out = sum(power{i_s}.p_out(band_of_bin == b, of_type), "all");
            p_all = sum(power{i_s}.p_all(band_of_bin == b, of_type), "all");
            if p_all > 0
                bounds(b, g, :) = [p_out / p_all, upper_level];
            end
        end
    end
    bounds(n_band, 1, :) = 1;                        % anchor: split cuts at 45-80 km
    obs{i_s}.bounds = bounds;
    lower_bound(:, :, i_s) = bounds(:, :, 1);
end

%% 3. Fits: before, levels, refined (with the shared sampling curve) -----------------------
fits = cell(n_scan, n_fit);
beta = ones(numel(lambda_knots), 1);
for i_s = 1:n_scan
    fits{i_s, 1} = fit_cut_levels(obs{i_s}, beam, knots_km, lambda_knots, beta, ones(n_band, 2, 2), false);
    fits{i_s, 2} = fit_cut_levels(obs{i_s}, beam, knots_km, lambda_knots, beta, obs{i_s}.bounds, false);
end
beta_rounds = nan(numel(beta), n_round);
for i_round = 1:n_round
    sum_observed = zeros(size(beta));
    sum_fitted = zeros(size(beta));
    for i_s = 1:n_scan
        fit = fit_cut_levels(obs{i_s}, beam, knots_km, lambda_knots, beta, obs{i_s}.bounds, true);
        fits{i_s, 3} = fit;
        sum_observed = sum_observed + fit.beta_observed;
        sum_fitted = sum_fitted + fit.beta_fitted;
    end
    % New curve: observed over fitted without beta, per lambda bin; 1 for many birds per gate
    beta = min(max(sum_observed ./ max(sum_fitted, eps), 0.1), 1.2);
    beta(end) = 1;
    beta(isnan(beta)) = 1;
    beta_rounds(:, i_round) = beta;
end
% Last fit with the final curve
for i_s = 1:n_scan
    fits{i_s, 3} = fit_cut_levels(obs{i_s}, beam, knots_km, lambda_knots, beta, obs{i_s}.bounds, true);
end
fprintf("beta at expected birds per gate %s:\n", mat2str(lambda_knots))
disp(round(beta_rounds, 2))

%% 4. Agreement between elevations ---------------------------------------------------
agreement = nan(n_scan, n_fit);
residual = nan(n_scan, 4, 2, n_fit);                 % mean log(observed / fitted) by band and type
level = nan(n_band, 2, n_scan, n_fit);
for i_s = 1:n_scan
    for k = 1:n_fit
        fit = fits{i_s, k};
        % 1 bird km^-3 added to both, so near-empty bins do not dominate
        ok = ~isnan(fit.density);
        log_ratio = log((fit.density + 1) ./ (fit.predicted + 1));
        within = ok & fit.range_m <= agreement_km * 1000;
        agreement(i_s, k) = sqrt(sum(fit.weight(within) .* log_ratio(within) .^ 2) / sum(fit.weight(within)));
        for b = 1:4
            for g = 1:2
                rows = ok & fit.band == b & fit.group == g;
                if any(rows)
                    residual(i_s, b, g, k) = sum(fit.weight(rows) .* log_ratio(rows)) / sum(fit.weight(rows));
                end
            end
        end
        level(:, :, i_s, k) = fit.level;
    end
end
filter_loss = 1 - squeeze(mean(lower_bound(2:3, 2, :), 1));
density_class = discretize(scans.bird_density, [0 20 100 Inf]);
result = table(scans.station, string(scans.night), scans.bird_density, filter_loss, agreement(:, 1), agreement(:, 2), agreement(:, 3), ...
    squeeze(level(1, 2, :, 2)), squeeze(level(1, 2, :, 3)), squeeze(level(2, 2, :, 2)), squeeze(level(2, 2, :, 3)), ...
    VariableNames=["station" "night" "bird_density" "batch_filter_loss" "rms_before" "rms_levels" "rms_refined" ...
    "batch_c_2_10_levels" "batch_c_2_10_refined" "batch_c_10_20_levels" "batch_c_10_20_refined"]);
writetable(result, fullfile(figure_dir, "GD_M_cut_levels.csv"))
save(fullfile(figure_dir, "GD_M_cut_levels.mat"), "level", "lower_bound", "agreement", "residual", "filter_loss", ...
    "band_edges_km", "knots_km", "lambda_knots", "beta", "beta_rounds")
fprintf("RMS log within %d km, median over scans: before %.3f, levels %.3f, refined %.3f\n", agreement_km, median(agreement, 1))
fprintf("Scans where the refined fit agrees better than before: %d of %d (levels: %d)\n", ...
    nnz(agreement(:, 3) < agreement(:, 1)), n_scan, nnz(agreement(:, 2) < agreement(:, 1)))
for k = 1:n_fit
    for g = 1:2
        fprintf("%-8s %-6s mean log(obs/fit) 2-10/10-20/20-30/30-45 km: %s\n", fit_name(k), group_name(g), sprintf("%6.2f", median(residual(:, :, g, k), 1, "omitnan")))
    end
end

%% 5. Examples: many scans, distance-height and the profile of every sweep ------------------
% Two scans per cell of filter loss (thirds) x density (below 20, 20-100,
% above 100 birds km^-3): the one whose elevations agree best and the one
% that agrees worst before the levels, plus the scans shown on 4 Nov.
loss_class = discretize(filter_loss, [-Inf prctile(filter_loss, [100/3 200/3]) Inf]);
example = strings(0);
for i_loss = 1:3
    for i_density = 1:3
        in_cell = find(loss_class == i_loss & density_class == i_density);
        if isempty(in_cell), continue, end
        [~, order] = sort(agreement(in_cell, 1));
        example = [example; scans.station(in_cell(order([1 end]))) + " " + string(scans.night(in_cell(order([1 end]))))]; %#ok<AGROW>
    end
end
example = unique([example; "KDLH 2021-10-16"; "KLNX 2022-03-27"; "KDGX 2024-08-21"; "KGRB 2023-10-14"; "KAMA 2021-04-03"; "KDVN 2022-05-07"], "stable");
example_id = arrayfun(@(e) find(scans.station + " " + string(scans.night) == e, 1), example);
[~, order] = sort(filter_loss(example_id));
example_id = example_id(order);
loss_name = ["low" "medium" "high"];
per_page = 6;
for i_page = 1:ceil(numel(example_id) / per_page)
    page = example_id((i_page - 1) * per_page + 1:min(i_page * per_page, numel(example_id)));
    fig = bmm4d.plot.create_export_figure(5000, 5400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
    t = tiledlayout(fig, per_page, 5, TileSpacing="tight", Padding="tight");
    for i_row = 1:numel(page)
        i_s = page(i_row);
        header = sprintf("%s %s", scans.station(i_s), string(scans.night(i_s)));
        note = sprintf("%.0f birds km^{-3}, %s filter loss (%.2f)", scans.bird_density(i_s), loss_name(loss_class(i_s)), filter_loss(i_s));
        plot_distance_height(nexttile(t, (i_row - 1) * 5 + 1), fits{i_s, 1}, obs{i_s}, false);
        ylabel(nexttile(t, (i_row - 1) * 5 + 1), ["\bf" + header "\rm" + note "Height above radar (km)"], FontSize=11)
        plot_distance_height(nexttile(t, (i_row - 1) * 5 + 2), fits{i_s, 3}, obs{i_s}, true);
        cb = colorbar(nexttile(t, (i_row - 1) * 5 + 2));
        cb.Label.String = "birds km^{-3}";
        for k = 1:n_fit
            plot_profiles(nexttile(t, (i_row - 1) * 5 + 2 + k), fits{i_s, k}, obs{i_s}, k > 1, agreement_km, agreement(i_s, k));
        end
        if i_row == 1
            titles = ["before" "refined: ÷ c β(λ)" "profiles, before" "profiles, levels only" "profiles, refined"];
            for j = 1:5
                title(nexttile(t, j), titles(j), FontSize=13)
            end
        end
    end
    bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, sprintf("GD_M_examples_%d.png", i_page)), Resolution=150)
    close(fig)
end

%% 6. Summary over all scans ---------------------------------------------------------
band_mid = (band_edges_km(1:4) + band_edges_km(2:5)) / 2;
group_colors = [0.25 0.25 0.25; 0.85 0.37 0.05];
fit_style = [":" "--" "-"];
fig = bmm4d.plot.create_export_figure(4400, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 4, TileSpacing="tight", Padding="tight");
ax = nexttile(t);
hold(ax, "on")
for k = 1:n_fit
    for g = 1:2
        plot(ax, band_mid, median(residual(:, :, g, k), 1, "omitnan"), fit_style(k) + "o", Color=group_colors(g, :), LineWidth=2, ...
            MarkerFaceColor=group_colors(g, :), DisplayName=group_name(g) + ", " + fit_name(k))
    end
end
yline(ax, 0, "k:", HandleVisibility="off")
set(ax, FontSize=12, XLim=[0 45], YLim=[-0.8 0.2])
grid(ax, "on")
xlabel(ax, "Distance (km)")
ylabel(ax, "Median log((observed + 1) / (fitted + 1))")
title(ax, "Do the elevations agree?", FontSize=14)
legend(ax, Location="southeast", FontSize=10, NumColumns=2)
ax = nexttile(t);
semilogx(ax, lambda_knots, beta, "-o", Color=[0.2 0.45 0.75], LineWidth=2, MarkerFaceColor=[0.2 0.45 0.75])
yline(ax, 1, "k:")
set(ax, FontSize=12, YLim=[0 1.25])
grid(ax, "on")
xlabel(ax, "Expected birds per gate \lambda")
ylabel(ax, "\beta(\lambda)")
title(ax, "Sampling curve, all scans", FontSize=14)
ax = nexttile(t);
hold(ax, "on")
band_colors = parula(5);
for b = 1:4
    scatter(ax, squeeze(lower_bound(b, 2, :)), squeeze(level(b, 2, :, 3)), 30, band_colors(b, :), "filled", MarkerFaceAlpha=0.7, ...
        DisplayName=sprintf("%g-%g km", band_edges_km(b), band_edges_km(b + 1)))
end
plot(ax, [0 1.1], [0 1.1], "k:", HandleVisibility="off")
set(ax, FontSize=12, XLim=[0 1.05], YLim=[0 1.15])
grid(ax, "on")
xlabel(ax, "Lower bound P_{out} / P_{all}")
ylabel(ax, "Batch level c, refined")
title(ax, "Batch levels within their bounds", FontSize=14)
legend(ax, Location="southeast", FontSize=10)
ax = nexttile(t);
hold(ax, "on")
plot(ax, [0 1], [0 1], "k:", HandleVisibility="off")
scatter(ax, agreement(:, 1), agreement(:, 2), 25, [0.6 0.6 0.6], "filled", DisplayName="levels only")
scatter(ax, agreement(:, 1), agreement(:, 3), 25, [0.85 0.37 0.05], "filled", DisplayName="refined")
set(ax, FontSize=12, XLim=[0 1], YLim=[0 1])
grid(ax, "on")
xlabel(ax, "Disagreement within 45 km, before")
ylabel(ax, "after")
title(ax, "Per scan (below the line: better)", FontSize=14)
legend(ax, Location="northwest", FontSize=10)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_M_summary.png"), Resolution=300)
close(fig)

%% Local functions ----------------------------------------------------------
function fit = fit_cut_levels(obs, beam, knots_km, lambda_knots, beta, bounds, pulled)
% Vertical profile (100 m layers to 4 km, smoothed over three layers at each
% step), levels c per cut type at KNOTS_KM (band centres) within BOUNDS
% (equal bounds fix them; PULLED adds one sweep's worth of
% pseudo-observations at the prediction) and the fixed sampling curve BETA
% at LAMBDA_KNOTS (expected birds per gate). Both are piecewise linear (c in
% distance, beta in log lambda, constant beyond the end knots), so that no
% step appears at a band or bin edge. Multiplicative updates as
% GD_FIT_RANGE_BIAS, gates weighted by the scan's mean density.
layer_edges = 0:100:4000;
n_layer = numel(layer_edges) - 1;
n_knot = numel(knots_km);
effective_radius = 4 / 3 * 6371000;
height_m = obs.antenna + obs.range_m .* sind(obs.elevation) + obs.range_m .^ 2 / (2 * effective_radius);
height_sd = obs.range_m * beam.width_rad / (4 * sqrt(log(2)));
W = normcdf((layer_edges(2:end) - height_m) ./ height_sd) - normcdf((layer_edges(1:end - 1) - height_m) ./ height_sd);
volume = beam.volume_sample_km3(obs.range_m);
scan_mean = sum(obs.density .* obs.count) / sum(obs.count);
w = obs.count / max(scan_mean, 1);
% Hat weights of each observation on the distance knots, and its cut type
H = interp1(knots_km(:), eye(n_knot), min(max(obs.range_m / 1000, knots_km(1)), knots_km(end)), "linear");
is_type = [obs.group == 1, obs.group == 2];
sweeps_at_knot = zeros(n_knot, 2);
for g = 1:2
    for k = 1:n_knot
        sweeps_at_knot(k, g) = numel(unique(obs.sweep(is_type(:, g) & H(:, k) > 0)));
    end
end
hat_lambda = @(lambda) interp1(log(lambda_knots(:)), eye(numel(lambda_knots)), ...
    min(max(log(lambda), log(lambda_knots(1))), log(lambda_knots(end))), "linear");
profile = scan_mean * ones(n_layer, 1);
level = ones(n_knot, 2);
for it = 1:300
    clean = W * profile;
    G = hat_lambda(max(clean .* volume, 1e-6));
    sampling = G * beta;
    level_obs = sum(H .* level(:, obs.group)', 2);
    factor = level_obs .* sampling;
    ratio = obs.density ./ max(factor .* clean, 1e-9);
    profile = profile .* ((w .* factor .* ratio)' * W)' ./ max(((w .* factor)' * W)', 1e-12);
    profile = conv([profile(1); profile; profile(end)], [1 2 1]' / 4, "valid");
    % Levels: multiplicative update per knot and type, optionally pulled toward 1
    clean = W * profile;
    base = sampling .* clean;
    level_obs = sum(H .* level(:, obs.group)', 2);
    ratio = obs.density ./ max(level_obs .* base, 1e-9);
    for g = 1:2
        rows = is_type(:, g);
        update_sum = H(rows, :)' * (w(rows) .* base(rows) .* ratio(rows));
        weight_sum = H(rows, :)' * (w(rows) .* base(rows));
        if pulled
            pseudo = weight_sum ./ max(sweeps_at_knot(:, g), 1);
            new_level = (level(:, g) .* update_sum + pseudo) ./ max(weight_sum + pseudo, eps);
        else
            new_level = level(:, g) .* update_sum ./ max(weight_sum, eps);
        end
        new_level(weight_sum == 0) = 1;
        level(:, g) = min(max(new_level, bounds(:, g, 1)), bounds(:, g, 2));
    end
end
clean = W * profile;
G = hat_lambda(max(clean .* volume, 1e-6));
sampling = G * beta;
level_obs = sum(H .* level(:, obs.group)', 2);
fit = struct("profile", profile, "layer_mid_m", (layer_edges(1:end - 1) + layer_edges(2:end))' / 2, "level", level, ...
    "factor", level_obs .* sampling, "predicted", level_obs .* sampling .* clean, "density", obs.density, ...
    "weight", w, "range_m", obs.range_m, "height_m", height_m, "band", obs.band, "group", obs.group, "sweep", obs.sweep, ...
    "beta_observed", G' * (w .* obs.density), "beta_fitted", G' * (w .* level_obs .* clean));
end

function plot_distance_height(ax, fit, obs, corrected)
% Each sweep's range bins as bands of the beam (+- half width), lowest drawn last
effective_radius = 4 / 3 * 6371000;
half_width_deg = 0.95 / 2;
values = fit.density;
if corrected
    values = values ./ fit.factor;
end
hold(ax, "on")
color_limit = [0 max(prctile(fit.density(fit.range_m <= 45000), 98), 1)];
for i_e = obs.n_sweep:-1:1
    rows = find(fit.sweep == i_e);
    if isempty(rows), continue, end
    elevation = obs.elevation(rows(1));
    r0 = obs.range_lo(rows);
    r1 = obs.range_hi(rows);
    height_at = @(r, e) (obs.antenna + r .* sind(e) + r .^ 2 / (2 * effective_radius)) / 1000;
    x = [r0 r1 r1 r0]' / 1000;
    y = [height_at(r0, elevation - half_width_deg) height_at(r1, elevation - half_width_deg) ...
        height_at(r1, elevation + half_width_deg) height_at(r0, elevation + half_width_deg)]';
    patch(ax, x, y, values(rows)', EdgeColor="none")
end
set(ax, XLim=[0 80], YLim=[0 3], FontSize=11, Color=[0.93 0.93 0.93], Box="on")
clim(ax, color_limit)
colormap(ax, parula)
end

function plot_profiles(ax, fit, obs, corrected, agreement_km, rms)
% Every sweep within AGREEMENT_KM as a line of density against beam height
% (grey: split cuts, orange: batch, darker = lower cut), and the fitted profile
type_colors = {[0.15 0.15 0.15; 0.6 0.6 0.6], [0.8 0.3 0.0; 0.98 0.75 0.4]};
values = fit.density;
if corrected
    values = values ./ fit.factor;
end
hold(ax, "on")
near = fit.range_m <= agreement_km * 1000;
for i_e = 1:obs.n_sweep
    rows = near & fit.sweep == i_e;
    if ~any(rows), continue, end
    shades = type_colors{2 - obs.is_split(i_e)};
    shade = shades(1, :) + (shades(2, :) - shades(1, :)) * (i_e - 1) / max(obs.n_sweep - 1, 1);
    plot(ax, values(rows), fit.height_m(rows) / 1000, "-", Color=shade, LineWidth=1.2)
end
plot(ax, fit.profile, fit.layer_mid_m / 1000, "k-", LineWidth=2.5)
x_max = max([prctile(values(near), 99) max(fit.profile(fit.layer_mid_m <= 3000)) 1]) * 1.1;
set(ax, YLim=[0 3], XLim=[0 x_max], FontSize=11, Box="on")
grid(ax, "on")
text(ax, 0.97, 0.95, sprintf("RMS log %.2f", rms), Units="normalized", HorizontalAlignment="right", FontSize=11)
if corrected
    text(ax, 0.97, 0.84, "batch c " + strjoin(compose("%.2f", fit.level(1:4, 2)'), " "), Units="normalized", ...
        HorizontalAlignment="right", FontSize=9, Color=[0.8 0.3 0.0])
    text(ax, 0.97, 0.75, "split c " + strjoin(compose("%.2f", fit.level(1:4, 1)'), " "), Units="normalized", ...
        HorizontalAlignment="right", FontSize=9, Color=[0.3 0.3 0.3])
end
end
