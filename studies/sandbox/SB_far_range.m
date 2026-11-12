%% Sandbox: the far-range dip of the radar-distance curve
% With the gate estimator and no cuts, the column density is 0.92 of the
% regional mean at 100-150 km from the radar. Three checks:
% 1. Solves: a flat prior (no radar imprint in the prior) and rows to
%    100 km only. If the dip goes with the far rows, the data cause it.
% 2. Held-out fit: rows beyond 100 km, predicted by the solve without
%    them (so from nearer radars), by frame and by the height the
%    prediction comes from.
% 3. Near-radar profile: the solved density by height within 50 km of a
%    radar, per frame, to compare with the far-range deficit.
% 4. Figure of the base solve: column map, column / 150 km regional mean,
%    and the radar-distance curve and vertical imprint in 1 km bins.
% Outputs: tables printed; figures/far_range_fine.png, .csv; results in
% studies/sandbox/README.md.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end
g = S.g;
obs = S.observations;
n_time = size(obs.values, 2);

%% 1. Solves --------------------------------------------------------------
runs = {sandbox_run(S, struct("name", "base", "video", false))
    sandbox_run(S, struct("name", "flat prior", "prior_mean", "profile", "video", false))
    sandbox_run(S, struct("name", "rows to 100 km", "row_weight", @(r, e, z) double(r <= 100000), "video", false))};
center_km = (runs{1}.edges_km(1:end - 1) + runs{1}.edges_km(2:end)) / 2;
curves = cell2mat(cellfun(@(r) mean(r.ratio, 2, "omitnan"), runs', UniformOutput=false));
disp(array2table(round(curves, 3), RowNames=string(center_km), VariableNames=cellfun(@(r) char(r.name), runs, UniformOutput=false)'))
disp(cellfun(@(r) r.bias_200km, runs)')

%% 2. Held-out fit of the rows beyond 100 km ------------------------------
% Prediction-weighted height of a row: where its predicted signal comes from
[~, ~, iz] = ind2sub(g.sz, S.active_id);
z_layer = squeeze(g.z);
z_cell = z_layer(iz);
F_active = obs.F(:, S.active_id);
far = obs.row_range_m > 100000 & obs.row_elev_deg < 1.1;
observed = full(obs.values(far, :));
held = F_active(far, :) * double(runs{3}.field);
row_height = (F_active(far, :) * (double(runs{3}.field) .* z_cell(:))) ./ max(held, eps);
use = observed ~= 0 & held > 1;
height_edges = 1000:250:2750;
height_bin = discretize(row_height, height_edges);
ratio_height = arrayfun(@(k) sum(observed(use & height_bin == k)) / sum(held(use & height_bin == k)), 1:numel(height_edges) - 1);
ratio_frame = arrayfun(@(t) sum(observed(use(:, t), t)) / sum(held(use(:, t), t)), 1:n_time);
disp(table(height_edges(1:end - 1)' + 125, ratio_height', VariableNames=["height_m" "observed_over_held_out"]))
disp(table(string(S.c.times(:), "HH:mm"), ratio_frame', VariableNames=["time_utc" "observed_over_held_out"]))

%% 3. Near-radar density profile per frame --------------------------------
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for i_r = 1:height(S.radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r)))));
end
near = distance_km(sub2ind(g.sz(1:2), ix, iy)) < 50;
profile = zeros(g.sz(3), n_time);
for i_t = 1:n_time
    profile(:, i_t) = accumarray(iz(near), double(runs{1}.field(near, i_t)), [g.sz(3) 1], @mean);
end
disp(array2table(round(profile(5:5:25, :)), RowNames=string(z_layer(5:5:25)) + " m", VariableNames=cellstr(string(S.c.times(:), "HH:mm"))))

%% 4. Reconstruction map and 1 km radar-distance curve ---------------------
% Distance of each grid column to the nearest radar, in 1 km bins. Grid
% columns are 4 km wide, so the curve cannot show anything finer than that.
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
column_id = sub2ind(g.sz(1:2), ix, iy);
is_active = false(g.sz(1:2));
is_active(column_id) = true;
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
edges_km = 0:1:250;
center_km = edges_km(1:end - 1)' + 0.5;
distance_bin = discretize(distance_km, edges_km);
height_edges = 0:200:3000;
cell_bin = distance_bin(column_id);
height_bin = discretize(z_cell(:), height_edges);
fine_runs = runs([1 3]);
curve = nan(numel(center_km), n_time, numel(fine_runs));
vertical = nan(numel(center_km), numel(height_edges) - 1, numel(fine_runs));
column_mean = zeros(g.sz(1:2));
relative_mean = zeros(g.sz(1:2));
for i_run = 1:numel(fine_runs)
    field = double(fine_runs{i_run}.field);
    for i_t = 1:n_time
        column = reshape(accumarray(column_id, field(:, i_t), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
        column = max(column, 0) .* is_active;
        regional = imgaussfilt(column, sigma_cells) ./ active_weight;
        relative = column ./ max(regional, 1);
        ok = is_active & ~isnan(distance_bin);
        curve(:, i_t, i_run) = accumarray(distance_bin(ok), relative(ok), [numel(center_km) 1], @mean, nan);
        if i_run == 1
            column_mean = column_mean + column / n_time;
            relative_mean = relative_mean + relative / n_time;
        end
    end
    % Vertical imprint: mean density per (distance, height) / mean at that height
    ok = ~isnan(cell_bin) & ~isnan(height_bin);
    cell_mean = mean(field, 2);
    by_cell = accumarray([cell_bin(ok) height_bin(ok)], cell_mean(ok), [numel(center_km) numel(height_edges) - 1], @mean, nan);
    layer_mean = accumarray(height_bin(ok), cell_mean(ok), [numel(height_edges) - 1 1], @mean, nan);
    vertical(:, :, i_run) = by_cell ./ layer_mean';
end
writetable(table(center_km, mean(curve(:, :, 1), 2, "omitnan"), mean(curve(:, :, 2), 2, "omitnan"), ...
    VariableNames=["distance_km" "base" "rows_to_100km"]), fullfile(figure_dir, "far_range_fine.csv"))

% Near range in 0.5 km bins: mean of the frame-mean ratio per radar, then
% the mean over radars and a 95% bootstrap interval over radars
own_radar = zeros(g.sz(1:2));
best_km = inf(g.sz(1:2));
for i_r = 1:height(S.radar_selected)
    d = 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r))));
    closer = d < best_km;
    best_km(closer) = d(closer);
    own_radar(closer) = i_r;
end
near_edges_km = 0:0.5:25;
near_center_km = near_edges_km(1:end - 1)' + 0.25;
near_bin = discretize(distance_km, near_edges_km);
ok = is_active & ~isnan(near_bin);
per_radar = accumarray([near_bin(ok) own_radar(ok)], relative_mean(ok), [numel(near_center_km) height(S.radar_selected)], @mean, nan);
near_mean = mean(per_radar, 2, "omitnan");
rng(1)
n_boot = 2000;
boot = nan(numel(near_center_km), n_boot);
for i_b = 1:n_boot
    boot(:, i_b) = mean(per_radar(:, randi(size(per_radar, 2), 1, size(per_radar, 2))), 2, "omitnan");
end
near_ci = prctile(boot, [2.5 97.5], 2);
n_near = sum(~isnan(per_radar), 2);

% Beam centre height (m a.s.l.) of each sweep at the mean radar height, 4/3 Earth
sweep_elev = [0.48 0.88 1.27 1.80 2.42 3.08 4.00 5.10 6.42];
beam_range_m = (0:1:250)' * 1000;
beam_height = mean(S.radar_selected.height) + beam_range_m .* sind(sweep_elev) + beam_range_m .^ 2 / (2 * 4 / 3 * 6371000);
beam_edges = mean(S.radar_selected.height) + beam_range_m .* sind(sweep_elev(1) + [-1 1] * 0.925 / 2) + beam_range_m .^ 2 / (2 * 4 / 3 * 6371000);

fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
title(t, "Sandbox, gate estimator, no cuts (mean over the 6 frames)", FontSize=14)
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
column_mean(~is_active) = nan;
relative_mean(~is_active) = nan;
ax = nexttile(t);
bmm4d.plot.plot_night_map(ax, g, background, column_mean, [0 200], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
cb = colorbar(ax); cb.Label.String = "Column density (birds km^{-2})";
title(ax, "Reconstruction", FontWeight="bold")
ax = nexttile(t);
bmm4d.plot.plot_night_map(ax, g, background, log2(relative_mean), [-0.5 0.5], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
colormap(ax, signed_map); cb = colorbar(ax); cb.Label.String = "log_2(column / 150 km regional mean)";
title(ax, "Radar imprint", FontWeight="bold")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
plot(ax, center_km, curve(:, :, 1), Color=[0.8 0.8 0.8])
plot(ax, center_km, mean(curve(:, :, 1), 2, "omitnan"), "k-", LineWidth=2)
plot(ax, center_km, mean(curve(:, :, 2), 2, "omitnan"), "-", Color=[0.85 0.33 0.1], LineWidth=1.5)
yline(ax, 1, "k--"); xlim(ax, [0 250]); ylim(ax, [0.6 1.4]); set(ax, FontSize=12)
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column / regional mean")
legend(ax, ["frames" repmat("", 1, n_time - 1) "mean" "rows to 100 km only"], Location="northeast")
title(ax, "Radar-distance curve, 1 km bins", FontWeight="bold")
for i_run = 1:2
    ax = nexttile(t, 3 + i_run);
    imagesc(ax, center_km, height_edges(1:end - 1) + 100, vertical(:, :, i_run)', AlphaData=~isnan(vertical(:, :, i_run)'))
    set(ax, YDir="normal", Color="w", FontSize=12); clim(ax, [0.5 1.5]); colormap(ax, signed_map)
    cb = colorbar(ax); cb.Label.String = "Density / mean at that height";
    xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Height (m a.s.l.)")
    title(ax, "Vertical imprint: " + fine_runs{i_run}.name, FontWeight="bold")
    % Beam centres of the sweeps (solid) and the half-power edges of the lowest (dashed)
    hold(ax, "on")
    plot(ax, beam_range_m / 1000, beam_height, "k-", LineWidth=0.8)
    plot(ax, beam_range_m / 1000, beam_edges, "k--", LineWidth=0.8)
    for i_s = 1:2
        i_label = find(beam_height(:, i_s) <= 2800 - 250 * i_s, 1, "last");
        text(ax, beam_range_m(i_label) / 1000 + 2, beam_height(i_label, i_s), sprintf("%.1f°", sweep_elev(i_s)), ...
            FontSize=11, HorizontalAlignment="left", VerticalAlignment="middle", BackgroundColor="w", Margin=0.5)
    end
    xlim(ax, [0 250]); ylim(ax, [0 3000])
end
ax = nexttile(t, 6); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
has_ci = ~isnan(near_ci(:, 1));
fill(ax, [near_center_km(has_ci); flipud(near_center_km(has_ci))], [near_ci(has_ci, 1); flipud(near_ci(has_ci, 2))], [0.8 0.8 0.8], ...
    EdgeColor="none", DisplayName="95% interval over radars")
plot(ax, near_center_km, near_mean, "k-o", LineWidth=1.5, MarkerSize=4, MarkerFaceColor="k", DisplayName="mean over radars")
yline(ax, 1, "k--", HandleVisibility="off"); xlim(ax, [0 25]); ylim(ax, [0.4 1.6]); set(ax, FontSize=12)
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column / regional mean")
legend(ax, Location="northwest")
title(ax, "First 25 km, 0.5 km bins (4 km grid)", FontWeight="bold")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "far_range_fine.png"), Resolution=150)
close(fig)
