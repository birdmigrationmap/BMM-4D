%% Sandbox: the near-radar artefact at fine range, per radar
% Three views of what happens close to each radar, with fine distance bins
% (1 km to 20 km, then 2.5, 5 and 10 km) and one line per radar:
% 1. Columns: column density / 150 km regional mean against the distance
%    to the column's own (nearest) radar, for the solution without cuts
%    and with hard cuts.
% 2. Observations: observed / predicted per elevation, where the
%    prediction comes from the hard-cut solution, which never saw the rows
%    within 25 km. Below 25 km this is the held-out near-range deficit.
% 3. Vertical: the reconstructed density by distance to the radar and
%    height above ground (mean over radars and frames), for both solutions
%    and their ratio, with the beam centre heights.
% Needs `results` from SB_fast_solver.m (or reruns the two solves).
% Outputs: studies/sandbox/figures/near_range_*.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end
run_names = ["no cuts, PCG" "hard cuts, PCG"];
if ~exist("results", "var") || ~all(ismember(run_names, cellfun(@(r) r.name, results)))
    hard_cut = @(r, e, z) double(r >= 25000 & ~(z & r > 75000));
    results = {sandbox_run(S, struct("name", run_names(1), "video", false))
        sandbox_run(S, struct("name", run_names(2), "row_weight", hard_cut, "video", false))};
end
runs = cellfun(@(name) results{find(cellfun(@(r) r.name == name, results), 1)}, cellstr(run_names), UniformOutput=false);
g = S.g;
n_radar = height(S.radar_selected);
n_time = size(S.observations.values, 2);
edges_km = [0:1:20, 22.5:2.5:50, 55:5:100, 110:10:150];
center_km = (edges_km(1:end - 1) + edges_km(2:end)) / 2;
n_bin = numel(center_km);

%% Distance of every column to its own (nearest) radar ----------------------
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
own_radar = zeros(g.sz(1:2));
for i_r = 1:n_radar
    d = 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r))));
    closer = d < distance_km;
    distance_km(closer) = d(closer);
    own_radar(closer) = i_r;
end
[ix, iy, iz] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
is_active = false(g.sz(1:2));
is_active(column_id) = true;
column_bin = discretize(distance_km, edges_km);

%% 1. Column / regional mean per radar --------------------------------------
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
column_ratio = nan(n_bin, n_radar, numel(runs));
for i_run = 1:numel(runs)
    sum_ratio = zeros(n_bin, n_radar);
    count = zeros(n_bin, n_radar);
    for i_t = 1:n_time
        column = reshape(accumarray(column_id, double(runs{i_run}.field(:, i_t)), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
        column = max(column, 0) .* is_active;
        regional = imgaussfilt(column, sigma_cells) ./ active_weight;
        ok = is_active & regional > 5 & ~isnan(column_bin);
        sum_ratio = sum_ratio + accumarray([column_bin(ok) own_radar(ok)], column(ok) ./ regional(ok), [n_bin n_radar]);
        count = count + accumarray([column_bin(ok) own_radar(ok)], 1, [n_bin n_radar]);
    end
    column_ratio(:, :, i_run) = sum_ratio ./ count;
end

%% 2. Held-out observed / predicted per elevation and radar --------------------
% Ratio of sums per (range bin, radar), all frames, rows with a prediction
% above 1 bird km^-3; no-echo rows (stored as eps) count as 0.
% Gates are 2 km apart from 4 km, so the rows get 2 km bins below 20 km
row_edges_km = [3:2:21, 22.5:2.5:50, 55:5:100, 110:10:150];
row_center_km = (row_edges_km(1:end - 1) + row_edges_km(2:end)) / 2;
n_row_bin = numel(row_center_km);
elevation_edges = [0 0.7 1.1 1.5 2.1 2.8 90];
elevation_label = ["0.5°" "0.9°" "1.3°" "1.8°" "2.4°" "≥ 3.1°"];
o = S.observations;
predicted = o.F(:, S.active_id) * double(runs{2}.field);
observed = full(o.values);
observed(observed < 1e-6) = 0;
has_value = full(o.values ~= 0);
row_bin = discretize(double(o.row_range_m) / 1000, row_edges_km);
row_elevation = discretize(double(o.row_elev_deg), elevation_edges);
n_elevation = numel(elevation_label);
held_out_ratio = nan(n_row_bin, n_radar, n_elevation);
for i_e = 1:n_elevation
    sum_observed = zeros(n_row_bin, n_radar);
    sum_predicted = zeros(n_row_bin, n_radar);
    for i_t = 1:n_time
        ok = has_value(:, i_t) & predicted(:, i_t) > 1 & row_elevation == i_e & ~isnan(row_bin);
        sum_observed = sum_observed + accumarray([row_bin(ok) double(o.row_radar_id(ok))], observed(ok, i_t), [n_row_bin n_radar]);
        sum_predicted = sum_predicted + accumarray([row_bin(ok) double(o.row_radar_id(ok))], predicted(ok, i_t), [n_row_bin n_radar]);
    end
    ratio = sum_observed ./ sum_predicted;
    ratio(sum_predicted < 50) = nan;
    held_out_ratio(:, :, i_e) = ratio;
end
clear predicted observed has_value

%% 3. Density by distance and height above ground -----------------------------
height_edges = 0:100:3000;
cell_agl = squeeze(g.z(iz)) - g.dem(column_id);
cell_agl = cell_agl(:);
cell_height_bin = discretize(cell_agl, height_edges);
cell_range_bin = column_bin(column_id);
ok = ~isnan(cell_height_bin) & ~isnan(cell_range_bin);
range_height = zeros(numel(height_edges) - 1, n_bin, numel(runs));
for i_run = 1:numel(runs)
    mean_density = mean(double(runs{i_run}.field), 2);
    range_height(:, :, i_run) = accumarray([cell_height_bin(ok) cell_range_bin(ok)], mean_density(ok), ...
        [numel(height_edges) - 1 n_bin], @mean, nan);
end
% Beam centre height above the radar's ground (4/3 earth radius, median
% antenna height above ground; the ground slope away from the radar is ignored)
beam_range_m = (0:1:150) * 1000;
antenna_agl = median(S.radar_selected.height - S.radar_selected.dem);
beam_agl = antenna_agl + beam_range_m(:) .* sind([0.5 0.9 1.3 1.8 2.4 3.1]) + beam_range_m(:) .^ 2 / (2 * 4 / 3 * 6371000);

%% Plot -------------------------------------------------------------------------
radar_colors = parula(n_radar + 4);
fig = create_export_figure(3600, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Column density / 150 km regional mean, per radar (mean over 6 frames)", FontSize=14, Color="k")
for i_run = 1:numel(runs)
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_r = 1:n_radar
        plot(ax, center_km, column_ratio(:, i_r, i_run), Color=[radar_colors(i_r, :) 0.5], LineWidth=1)
    end
    plot(ax, center_km, median(column_ratio(:, :, i_run), 2, "omitnan"), "k-", LineWidth=3)
    yline(ax, 1, "k--"); xline(ax, 25, ":", "25 km cut")
    set(ax, XScale="log", FontSize=12); xlim(ax, [0.5 150]); ylim(ax, [0 2])
    xticks(ax, [1 2 5 10 15 25 50 100 150])
    xlabel(ax, "Distance to own radar (km)"); ylabel(ax, "Column / regional mean")
    title(ax, runs{i_run}.name + " (black: median over radars)", FontWeight="normal", Interpreter="none")
end
export_fixed_figure(fig, fullfile(figure_dir, "near_range_columns.png"), Resolution=150)
close(fig)

fig = create_export_figure(3600, 2000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
title(t, "Observed / predicted per radar, prediction from the hard-cut solution (rows < 25 km held out)", FontSize=14, Color="k")
for i_e = 1:n_elevation
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_r = 1:n_radar
        plot(ax, row_center_km, held_out_ratio(:, i_r, i_e), Color=[radar_colors(i_r, :) 0.5], LineWidth=1)
    end
    plot(ax, row_center_km, median(held_out_ratio(:, :, i_e), 2, "omitnan"), "k-", LineWidth=3)
    yline(ax, 1, "k--"); xline(ax, 25, ":")
    set(ax, XScale="log", YScale="log", FontSize=12); xlim(ax, [3 150]); ylim(ax, [0.01 4])
    xticks(ax, [4 6 10 15 25 50 100 150]); yticks(ax, [0.01 0.03 0.1 0.3 1 3])
    xlabel(ax, "Range (km)"); ylabel(ax, "Observed / predicted")
    title(ax, "Elevation " + elevation_label(i_e), FontWeight="normal")
end
export_fixed_figure(fig, fullfile(figure_dir, "near_range_observations.png"), Resolution=150)
close(fig)

fig = create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "Reconstructed density by distance to own radar and height above ground (mean over radars and frames)", FontSize=14, Color="k")
height_center = (height_edges(1:end - 1) + height_edges(2:end)) / 2;
panels = {range_height(:, :, 1), range_height(:, :, 2), log2(range_height(:, :, 1) ./ range_height(:, :, 2))};
panel_titles = [run_names(1) run_names(2) "log_2(no cuts / hard cuts)"];
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
for i_p = 1:3
    ax = nexttile(t); hold(ax, "on"); box(ax, "on")
    pcolor(ax, center_km, height_center, panels{i_p}); shading(ax, "flat")
    plot(ax, beam_range_m / 1000, beam_agl, "w-", LineWidth=1)
    set(ax, XScale="log", FontSize=12, Layer="top"); xlim(ax, [0.5 150]); ylim(ax, [0 3000])
    xticks(ax, [1 2 5 10 25 50 150]); xtickangle(ax, 0)
    xlabel(ax, "Distance to own radar (km)"); ylabel(ax, "Height above ground (m)")
    cb = colorbar(ax, Color="k");
    if i_p < 3
        clim(ax, [0 quantile(range_height(:), 0.99)]); cb.Label.String = "Density (birds km^{-3})";
    else
        colormap(ax, signed_map); clim(ax, [-1 1]); cb.Label.String = "log_2 ratio";
    end
    title(ax, panel_titles(i_p) + " (white: beam centres 0.5-3.1°)", FontWeight="normal", Interpreter="none")
end
export_fixed_figure(fig, fullfile(figure_dir, "near_range_vertical.png"), Resolution=150)
close(fig)
