%% Sandbox: a radar-profile prior that keeps the low-level peak
% The prior mean fits one profile per radar and time to all its rows, with
% a first-difference smoothness in height (lambda_profile, 10 in
% production). The profile was too low near the ground and too high
% aloft, and between radars the solution falls back to it (README, "Where
% the column goes with distance"). Weighting the rows by range did not
% change it; the vertical smoothness does. Here lambda_profile = 10 (current),
% 1 and 0.1.
% For each value: the prior and the solution as column contribution by height
% band and distance to the nearest radar, the radar-distance curve (1 km
% bins), bias and misfits.
% Outputs: figures/profile_prior.png, profile_prior.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var") || ~isstruct(S)
    S = sandbox_load("sandbox");
end
g = S.g;
lambda_values = [10 1 0.1];

%% Solve ------------------------------------------------------------------
runs = cell(numel(lambda_values), 1);
priors = cell(numel(lambda_values), 1);
inv = S.c.inversion;
for i_q = 1:numel(lambda_values)
    runs{i_q} = sandbox_run(S, struct("name", sprintf("lambda_profile %g", lambda_values(i_q)), ...
        "lambda_profile", lambda_values(i_q), "video", false));
    priors{i_q} = S.prior_cache(sprintf("%g|%s|%g", lambda_values(i_q), inv.prior_interpolation, inv.prior_length_scale_m));
end

%% Column contribution by height band and distance ------------------------------
[ix, iy, iz] = ind2sub(g.sz, S.active_id);
z_layer = squeeze(g.z);
z_cell = z_layer(iz);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for i_r = 1:height(S.radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r)))));
end
column_id = sub2ind(g.sz(1:2), ix, iy);
distance_edges = [0 25 50 100 150 250];
height_edges = [0 1000 2000 3000];
distance_bin = discretize(distance_km(column_id), distance_edges);
height_bin = discretize(z_cell(:), height_edges);
n_columns = arrayfun(@(b) numel(unique(column_id(distance_bin == b))), 1:numel(distance_edges) - 1)';
ok = ~isnan(distance_bin) & ~isnan(height_bin);
band = @(field) accumarray([distance_bin(ok) height_bin(ok)], field(ok), [numel(distance_edges) - 1 numel(height_edges) - 1]) ...
    * g.dz / 1000 ./ n_columns;
summary = table();
for i_q = 1:numel(lambda_values)
    for source = ["prior" "solution"]
        if source == "prior"
            field = mean(double(priors{i_q}), 2);
        else
            field = mean(double(runs{i_q}.field), 2);
        end
        contribution = band(field);
        summary = [summary; table(repmat(lambda_values(i_q), numel(n_columns), 1), repmat(source, numel(n_columns), 1), ...
            string(distance_edges(1:end - 1)') + "-" + string(distance_edges(2:end)'), contribution(:, 1), contribution(:, 2), ...
            contribution(:, 3), sum(contribution, 2), VariableNames=["lambda_profile" "source" "distance_km" "below_1km" ...
            "from_1_to_2km" "above_2km" "column"])]; %#ok<AGROW>
    end
end
summary{:, 4:7} = round(summary{:, 4:7}, 1);
disp(summary)
writetable(summary, fullfile(figure_dir, "profile_prior.csv"))
scores = table(lambda_values', cellfun(@(r) r.bias_200km, runs), cellfun(@(r) r.misfit, runs), cellfun(@(r) r.misfit_trusted, runs), ...
    cellfun(@(r) r.detail, runs), VariableNames=["lambda_profile" "bias" "misfit" "trusted_misfit" "detail"]);
disp(scores)

%% Plot results -----------------------------------------------------------
% Left: radar-distance curve in 1 km bins. Right: share of the column
% below 1 km by distance, prior (dashed) and solution (solid).
is_active = false(g.sz(1:2));
is_active(column_id) = true;
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
fine_edges = 0:1:250;
fine_center = fine_edges(1:end - 1)' + 0.5;
fine_bin = discretize(distance_km, fine_edges);
low = z_cell(:) < 1000;
colors = lines(numel(lambda_values));
fig = bmm4d.plot.create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Vertical smoothness of the radar-profile prior (estimator, no cuts; mean over frames)", FontSize=14)
ax_curve = nexttile(t); hold(ax_curve, "on"); box(ax_curve, "on"); grid(ax_curve, "on")
ax_low = nexttile(t); hold(ax_low, "on"); box(ax_low, "on"); grid(ax_low, "on")
for i_q = 1:numel(lambda_values)
    field = mean(double(runs{i_q}.field), 2);
    column = reshape(accumarray(column_id, field, [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    column = max(column, 0) .* is_active;
    regional = imgaussfilt(column, sigma_cells) ./ active_weight;
    use = is_active & ~isnan(fine_bin);
    curve = accumarray(fine_bin(use), column(use) ./ max(regional(use), 1), [numel(fine_center) 1], @mean, nan);
    plot(ax_curve, fine_center, curve / mean(curve, "omitnan"), Color=colors(i_q, :), LineWidth=1.5, ...
        DisplayName=sprintf("lambda_{profile} = %g: bias %.3f, trusted misfit %.3f", lambda_values(i_q), runs{i_q}.bias_200km, runs{i_q}.misfit_trusted))
    for source = 1:2
        if source == 1
            values = field;
        else
            values = mean(double(priors{i_q}), 2);
        end
        low_column = accumarray(column_id, values .* low, [prod(g.sz(1:2)) 1]);
        all_column = accumarray(column_id, values, [prod(g.sz(1:2)) 1]);
        share = accumarray(fine_bin(use), low_column(use) ./ max(all_column(use), eps), [numel(fine_center) 1], @mean, nan);
        styles = ["-" "--"];
        plot(ax_low, fine_center, movmean(share, 5, "omitnan"), styles(source), Color=colors(i_q, :), LineWidth=1.5, ...
            HandleVisibility=string(matlab.lang.OnOffSwitchState(source == 1)), DisplayName=sprintf("lambda_{profile} = %g", lambda_values(i_q)))
    end
end
yline(ax_curve, 1, "k--", HandleVisibility="off"); xlim(ax_curve, [0 250]); ylim(ax_curve, [0.6 1.4]); set(ax_curve, FontSize=12)
xlabel(ax_curve, "Distance to the nearest radar (km)"); ylabel(ax_curve, "Column / regional mean")
legend(ax_curve, Location="northeast", FontSize=10); title(ax_curve, "Radar-distance curve, 1 km bins", FontWeight="normal")
xlim(ax_low, [0 250]); set(ax_low, FontSize=12)
xlabel(ax_low, "Distance to the nearest radar (km)"); ylabel(ax_low, "Share of the column below 1 km")
legend(ax_low, Location="northeast", FontSize=10)
title(ax_low, "Low-layer share: solution (solid), prior (dashed)", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "profile_prior.png"), Resolution=150)
close(fig)
