% SB_DOUBLE_LAYER How many vertical shape patterns do the radar profiles need?
%   On the us_east night (101 radars x 26 frames), the radar profiles are
%   shifted to height above the smoothed ground, normalized to unit column,
%   and decomposed into shape patterns (EOFs). For K = 1..6 patterns it
%   reports the shape variance explained and the misplaced share (the part
%   of the column put at the wrong height by the K-pattern fit), and it
%   counts double-layer profiles (two peaks >= 500 m apart, the secondary
%   >= 50% of the primary). Output: figures/double_layer_*.png and the
%   table figures/double_layer_profiles.csv.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "sandbox"), fullfile(project_root, "studies", "two_level", "functions"))
out_dir = fullfile(project_root, "studies", "two_level", "figures", "sandbox");

%% Read the us_east case and fit the radar profiles
c = prod_case("us_east");
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
% active_id = 1: only the per-radar profiles are needed, not the gridded prior
[~, prior_info] = bmm4d.prior.build_density_prior_fast(g, radar_scans, observations, "lambda_profile", 1, ...
    "active_id", 1, "interpolation", "gaussian", "include_water", true);
profiles_asl = prior_info.radar_profiles;   % [nz n_radar n_time], birds km^-3, layers above sea level
clear observations
[nz, n_radar, n_time] = size(profiles_asl);
times = c.times(:);

%% Shift the profiles to height above the smoothed ground
% Shapes are defined above the DEM low-pass filtered over ~30 km (decision 1
% in MODEL_REVIEW): a radar profile averages beams over ~100 km, so the
% smoothed ground is its natural reference. Profiles are resampled on
% 50:100:2450 m AGL; radars whose grid top (2950 m asl) leaves less than
% 2.5 km above ground are dropped (high plains radars).
ground_smooth = imgaussfilt(max(g.dem, 0), 30000 / mean([g.dx g.dy]));
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
radar_ground = interp2(lon_grid', lat_grid', ground_smooth', radar_selected.lon, radar_selected.lat);
z_asl = squeeze(g.z);
z_agl = (50:100:2450)';
n_agl = numel(z_agl);
profiles_agl = nan(n_agl, n_radar, n_time);
for i_r = 1:n_radar
    if z_asl(end) - radar_ground(i_r) < z_agl(end)
        continue
    end
    profiles_agl(:, i_r, :) = max(interp1(z_asl - radar_ground(i_r), reshape(profiles_asl(:, i_r, :), nz, n_time), z_agl, "linear", 0), 0);
end
vid = squeeze(sum(profiles_agl, 1) * g.dz / 1000);   % birds km^-2, [n_radar n_time]

%% Normalized shapes of the profiles with birds
% Keep radar-times with VID >= 5 birds km^-2 (below that the shape is noise)
[radar_id, time_id] = ndgrid(1:n_radar, 1:n_time);
keep = ~isnan(vid) & vid >= 5;
P = reshape(profiles_agl, n_agl, []);
P = P(:, keep(:));
shapes = P ./ sum(P, 1);                     % unit column (sum over layers = 1)
weights = vid(keep);                         % for density-weighted summaries
n_profile = size(shapes, 2);
fprintf("%d radars kept (ground <= %.0f m), %d radar-times with VID >= 5\n", ...
    nnz(any(~isnan(vid), 2)), z_asl(end) - z_agl(end), n_profile)

%% Shape patterns (EOFs) and how well K patterns represent the profiles
% Uncentered SVD: the solver fits rho = sum_k x_k phi_k with free x_k, so
% the K-pattern fit of a profile is its projection on the first K patterns.
% phi_1 is then close to the mean shape.
[U, sv, ~] = svd(shapes, "econ");
U = U .* sign(sum(U, 1));                    % phi_1 positive
eof = U(:, 1:6);
k_list = 1:6;
variance_explained = cumsum(diag(sv) .^ 2) / sum(diag(sv) .^ 2);
% Centered variance (anomalies from the mean shape), for comparison with
% the 76/91/96% quoted earlier
[~, sv_centered, ~] = svd(shapes - mean(shapes, 2), "econ");
variance_centered = cumsum(diag(sv_centered) .^ 2) / sum(diag(sv_centered) .^ 2);
misplaced = nan(numel(k_list), n_profile);
for k = k_list
    fit = eof(:, 1:k) * (eof(:, 1:k)' * shapes);
    fit = max(fit, 0);
    fit = fit ./ sum(fit, 1);
    % Misplaced share: half the L1 distance between the two unit-column
    % shapes, i.e. the fraction of birds the K-pattern fit puts at the wrong height
    misplaced(k, :) = 0.5 * sum(abs(fit - shapes), 1);
end
summary = table(k_list', variance_explained(k_list), [nan; variance_centered(1:5)], ...
    median(misplaced, 2), prctile(misplaced, 90, 2), prctile(misplaced, 99, 2), ...
    sum(misplaced .* weights', 2) / sum(weights), mean(misplaced > 0.2, 2), ...
    VariableNames=["K" "variance" "variance_centered_Kminus1" "misplaced_median" "misplaced_p90" ...
    "misplaced_p99" "misplaced_vid_weighted" "share_misplaced_gt20pct"]);
disp(summary)
writetable(summary, fullfile(out_dir, "double_layer_eof_summary.csv"))

%% Spread of the shape coefficients and saved patterns
% c_k = x_k / x_1 for each profile: its spread sets the prior sd of pattern
% k in the separable solver (penalty growing with pattern order).
coefficients = eof' * shapes;
shape_sd = std(coefficients(2:6, :) ./ coefficients(1, :), 0, 2);
disp(table((2:6)', shape_sd, VariableNames=["k" "sd_of_c_k"]))
% Saved for the sandbox solver, scaled so that sum(phi_1) * dz = 1 km:
% x_1 is then the column amount (birds km^-2) and x_k / x_1 = c_k
shape_patterns = eof / (sum(eof(:, 1)) * 0.1);
save(fullfile(c.data_dir, "shape_patterns.mat"), "shape_patterns", "z_agl", "shape_sd")

%% Double-layer profiles
% Peaks of the shape with a prominence of at least 10% of the maximum;
% double layer = a secondary peak >= 50% of the primary, >= 500 m away.
is_double = false(1, n_profile);
peak_low = nan(1, n_profile);
peak_high = nan(1, n_profile);
secondary_ratio = nan(1, n_profile);
for i_p = 1:n_profile
    s = shapes(:, i_p);
    % Zero padding so a peak in the lowest or highest layer counts
    peak_layer = find(islocalmax([0; s; 0], MinProminence=0.1 * max(s), FlatSelection="center")) - 1;
    peak_value = s(peak_layer);
    if numel(peak_value) < 2
        continue
    end
    [~, order] = sort(peak_value, "descend");
    primary = order(1);
    for j = order(2:end)'
        if peak_value(j) >= 0.5 * peak_value(primary) && abs(z_agl(peak_layer(j)) - z_agl(peak_layer(primary))) >= 500
            is_double(i_p) = true;
            pair = sort(z_agl(peak_layer([primary j])));
            peak_low(i_p) = pair(1);
            peak_high(i_p) = pair(2);
            secondary_ratio(i_p) = peak_value(j) / peak_value(primary);
            break
        end
    end
end
kept_radar = radar_id(keep);
kept_time = time_id(keep);
double_table = table(radar_selected.name(kept_radar(is_double)), times(kept_time(is_double)), ...
    vid(sub2ind([n_radar n_time], kept_radar(is_double), kept_time(is_double))), ...
    peak_low(is_double)', peak_high(is_double)', secondary_ratio(is_double)', ...
    misplaced(3, is_double)', misplaced(4, is_double)', misplaced(6, is_double)', ...
    VariableNames=["radar" "time" "vid" "peak_low_m" "peak_high_m" "secondary_ratio" ...
    "misplaced_K3" "misplaced_K4" "misplaced_K6"]);
double_table = sortrows(double_table, ["radar" "time"]);
writetable(double_table, fullfile(out_dir, "double_layer_profiles.csv"))
fprintf("Double layers: %d of %d radar-times (%.1f%%), %.1f%% of the birds (VID-weighted), on %d radars\n", ...
    nnz(is_double), n_profile, 100 * mean(is_double), 100 * sum(weights(is_double)) / sum(weights), ...
    numel(unique(double_table.radar)))
disp(groupsummary(double_table, "radar", ["min" "max"], "time"))
fprintf("Misplaced share on double layers, median K=3/4/6: %.3f / %.3f / %.3f; single layers K=3: %.3f\n", ...
    median(misplaced(3, is_double)), median(misplaced(4, is_double)), median(misplaced(6, is_double)), ...
    median(misplaced(3, ~is_double)))

%% Plot results
fig = bmm4d.plot.create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 4, TileSpacing="tight", Padding="tight");
% Shape patterns
ax = nexttile(t); hold(ax, "on"); box(ax, "on")
colors = lines(6);
for k = 1:6
    plot(ax, eof(:, k), z_agl / 1000, LineWidth=2, Color=colors(k, :), DisplayName="\phi_" + k)
end
xline(ax, 0, "k:", HandleVisibility="off")
legend(ax, Location="northeast"); xlabel(ax, "Pattern value"); ylabel(ax, "Height above smoothed ground (km)")
title(ax, "Shape patterns (uncentered EOFs)"); set(ax, FontSize=13)
% Variance and misplaced share against K
ax = nexttile(t); hold(ax, "on"); box(ax, "on")
plot(ax, k_list, 1 - variance_explained(k_list), "-o", Color=[0.5 0.5 0.5], LineWidth=2, DisplayName="unexplained variance")
plot(ax, k_list, summary.misplaced_median, "-o", Color=[0.85 0.33 0.1], LineWidth=2, DisplayName="misplaced, median")
plot(ax, k_list, summary.misplaced_p90, "--o", Color=[0.85 0.33 0.1], LineWidth=2, DisplayName="misplaced, 90%")
plot(ax, k_list, median(misplaced(:, is_double), 2), "-s", Color=[0.2 0.2 0.7], LineWidth=2, DisplayName="misplaced, double layers")
set(ax, YScale="log", FontSize=13); xlabel(ax, "Number of patterns K"); ylabel(ax, "Share")
legend(ax, Location="southwest"); title(ax, "Representation error")
% Double-layer examples with their K = 3 and K = 4 fits
[~, order] = sort(misplaced(3, :) .* is_double, "descend");
for i_example = 1:2
    ax = nexttile(t); hold(ax, "on"); box(ax, "on")
    i_p = order(i_example);
    s = shapes(:, i_p);
    for k = [1 3 4]
        fit = max(eof(:, 1:k) * (eof(:, 1:k)' * s), 0);
        plot(ax, fit / sum(fit), z_agl / 1000, LineWidth=1.5, DisplayName="K = " + k)
    end
    plot(ax, s, z_agl / 1000, "k-", LineWidth=2.5, DisplayName="profile")
    legend(ax, Location="northeast"); xlabel(ax, "Shape (unit column)"); set(ax, FontSize=13)
    title(ax, sprintf("%s %s UTC, VID %.0f", radar_selected.name(kept_radar(i_p)), ...
        string(times(kept_time(i_p)), "HH:mm"), weights(i_p)))
end
bmm4d.plot.export_fixed_figure(fig, fullfile(out_dir, "double_layer_eofs.png"), Resolution=150)
close(fig)
