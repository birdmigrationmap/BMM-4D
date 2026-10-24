%% Production test 5: maps and animation
% Vertically integrated density: a nine-time overview PNG and an MP4 of
% every time step, on a light land/water background.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"))
% Case: "keywest", "midwest" or "us_east"; suffix _rho or _none for other rain
% removal. A case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
if ~isfolder(c.figure_dir)
    mkdir(c.figure_dir)
end
load(fullfile(c.data_dir, "case.mat"), "g", "radar_selected")
load(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id")
n_time = size(active_field, 2);

%% Vertically integrated density ------------------------------------------
% Density is in birds km^-3 and g.dz in m, so the column is birds km^-2.
% Columns without any active cell are NaN (background shows through).
column_density = zeros([g.sz(1:2) n_time], "single");
for i_t = 1:n_time
    field = zeros(g.sz, "single");
    field(active_id) = active_field(:, i_t);
    column_density(:, :, i_t) = sum(field, 3) .* g.dz ./ 1000;
end
active_column = false(g.sz(1:2));
[ix, iy, ~] = ind2sub(g.sz, active_id);
active_column(sub2ind(g.sz(1:2), ix, iy)) = true;
column_density(repmat(~active_column, 1, 1, n_time)) = nan;
column_density(column_density < 0) = 0;
c_lim = [0 quantile(column_density(:), 0.995)];

% Fade the map to 40% opacity from the end of the radar data (r_lim) to the
% edge of the solved domain, where the field comes only from the prior and
% advection.
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
radar_distance = inf(g.sz(1:2));
for i_r = 1:height(radar_selected)
    d_km = 6371 * acos(min(1, sind(lat_grid) .* sind(radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(radar_selected.lat(i_r)) .* cosd(lon_grid - radar_selected.lon(i_r))));
    radar_distance = min(radar_distance, d_km * 1000);
end
alpha = 0.4 + 0.6 * min(1, max(0, (c.active_range_m - radar_distance) ./ (c.active_range_m - c.r_lim)));
if ~isempty(c.domain)
    % A geographic domain is shown uniformly, without the radar-range fade
    alpha = ones(g.sz(1:2));
end

% Hide cells where the sun is up (night window of prod_case), fading over
% +-1 deg of sun elevation around the limit.
night_alpha = zeros([g.sz(1:2) n_time], "single");
for i_t = 1:n_time
    elevation = sun_elevation(lat_grid, lon_grid, c.times(i_t));
    night_alpha(:, :, i_t) = min(1, max(0, (c.sun_limit_deg + 1 - elevation) / 2));
end
frame_labels = string(c.times, "dd MMM HH:mm") + " UTC";

%% Map background ---------------------------------------------------------
% Pale water and light land keep the density colors dominant.
land_color = [0.96 0.96 0.94];
water_color = [0.85 0.92 0.97];
background = reshape(land_color, 1, 1, 3) .* ~g.mask_water + reshape(water_color, 1, 1, 3) .* g.mask_water;
background = permute(background, [2 1 3]);  % image() expects [y x]
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
map_bounds = c.figure_bounds;

%% Nine-time overview -----------------------------------------------------
frame_ids = unique(round(linspace(1, n_time, 9)));
fig = create_export_figure(3600, 2600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, 3, TileSpacing="tight", Padding="tight");
for i_frame = 1:numel(frame_ids)
    ax = nexttile(t);
    plot_night_map(ax, g, background, column_density(:, :, frame_ids(i_frame)), c_lim, states, radar_selected, map_bounds, ...
        alpha .* night_alpha(:, :, frame_ids(i_frame)));
    title(ax, frame_labels(frame_ids(i_frame)), FontSize=14, FontWeight="bold")
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Label.String = "Vertically integrated density (birds km^{-2})";
set(cb, FontName="Helvetica", FontSize=13, Color="k")
export_fixed_figure(fig, fullfile(c.figure_dir, "column_density.png"), Resolution=300)
close(fig)

%% Animation of all times -------------------------------------------------
video = VideoWriter(fullfile(c.figure_dir, "column_density.mp4"), "MPEG-4");
video.FrameRate = 4;
open(video)
fig = create_export_figure(2400, 1800, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    ax = axes(fig);
    plot_night_map(ax, g, background, column_density(:, :, i_t), c_lim, states, radar_selected, map_bounds, ...
        alpha .* night_alpha(:, :, i_t));
    title(ax, frame_labels(i_t), FontSize=14, FontWeight="bold")
    cb = colorbar(ax, Color="k");
    cb.Label.String = "Vertically integrated density (birds km^{-2})";
    writeVideo(video, getframe(fig))
end
close(video)
close(fig)
fprintf("[prod/%s] figures in %s\n", c.name, c.figure_dir)
