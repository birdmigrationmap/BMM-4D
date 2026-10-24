%% US night reconstruction: maps and animation
% Map the vertically integrated bird density for the reconstructed night
% and export a nine-time overview and an MP4 of every time step.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "us_night", "functions"))

night_date = datetime("2026-10-21");
night_tag = string(night_date, "yyyyMMdd");
night_dir = fullfile(bmm4d.util.local_data_dir(), "us_night", "data", night_tag);
figure_dir = fullfile(project_root, "us_night", "figures", night_tag);
if ~exist(figure_dir, "dir")
    mkdir(figure_dir)
end

load(fullfile(night_dir, "case.mat"), "g", "radar_selected")
load(fullfile(night_dir, "inversion.mat"), "active_field", "active_id")
n_time = size(active_field, 2);

%% Vertically integrated density ------------------------------------------
% Density is in birds km^-3 and g.dz in m, so the column is birds km^-2.
column_density = nan([g.sz(1:2) n_time], "single");
active_xy = false(g.sz(1:2));
[x_id, y_id, ~] = ind2sub(g.sz, active_id);
active_xy(sub2ind(g.sz(1:2), x_id, y_id)) = true;
for i_t = 1:n_time
    field = zeros(g.sz, "single");
    field(active_id) = active_field(:, i_t);
    column_t = sum(field, 3) .* g.dz ./ 1000;
    column_t(~active_xy) = nan;
    column_density(:, :, i_t) = column_t;
end
c_lim = [0 quantile(column_density(:), 0.99)];
frame_labels = string(g.time, "dd MMM HH:mm") + " UTC";

%% Map background ---------------------------------------------------------
% Pale water and light land keep the density colors dominant.
land_color = [0.96 0.96 0.94];
water_color = [0.85 0.92 0.97];
background = reshape(land_color, 1, 1, 3) .* ~g.mask_water + reshape(water_color, 1, 1, 3) .* g.mask_water;
background = permute(background, [2 1 3]);  % image() expects [y x]
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
map_bounds = [g.lon(1) g.lon(end) g.lat(1) g.lat(end)];

%% Nine-time overview -----------------------------------------------------
frame_ids = unique(round(linspace(1, n_time, 9)));
fig = bmm4d.plot.create_export_figure(3600, 2600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, 3, TileSpacing="tight", Padding="tight");
for i_frame = 1:numel(frame_ids)
    ax = nexttile(t);
    bmm4d.plot.plot_night_map(ax, g, background, column_density(:, :, frame_ids(i_frame)), c_lim, states, radar_selected, map_bounds);
    title(ax, frame_labels(frame_ids(i_frame)), FontSize=14, FontWeight="bold")
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Label.String = "Vertically integrated density (birds km^{-2})";
set(cb, FontName="Helvetica", FontSize=13, Color="k")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "column_density_" + night_tag + ".png"), Resolution=300)
close(fig)

%% Animation of all times -------------------------------------------------
video = VideoWriter(fullfile(figure_dir, "column_density_" + night_tag + ".mp4"), "MPEG-4");
video.FrameRate = 2;
open(video)
fig = bmm4d.plot.create_export_figure(2400, 1800, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    ax = axes(fig);
    bmm4d.plot.plot_night_map(ax, g, background, column_density(:, :, i_t), c_lim, states, radar_selected, map_bounds);
    title(ax, frame_labels(i_t), FontSize=14, FontWeight="bold")
    cb = colorbar(ax, Color="k");
    cb.Label.String = "Vertically integrated density (birds km^{-2})";
    writeVideo(video, getframe(fig))
end
close(video)
close(fig)
fprintf("[UN/D] exported figures to %s\n", figure_dir)
