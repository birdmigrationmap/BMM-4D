%% Artefacts F: video of the us_east night with the fixes, and the five spots (dev)
% The field of AR_E_us_east.m (prod/us_east_fix/inversion_fix.mat) drawn as
% the round-0 video of studies/sandbox/SB_us_east_video.m (column density,
% night mask, 5 min advection morph between the 30 min frames) with the
% same colour scale, so the two videos can be compared frame by frame.
% Then the five places of the screenshots, round 0 against the fixes:
% High Plains 05:30, Appalachians 05:30, Missouri 05:30, central Texas
% 02:30, KDMX 06:30.
% Outputs: figures/us_east_night_<output_name>.mp4 (+ middle frame .png),
% figures/AR_F_spots_<output_name>.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
c = case_settings("us_east");
% The reconstruction to draw (AR_E_us_east.m): "us_east_fix" (per-volume
% aspect fit, section 3) or "us_east_shared" (shared aspect function, section 4)
fix_folder = "us_east_fix";
output_name = "fix";
c_fix = ar_test_case(fix_folder, [-104 -66 24 50], c.times, "read_options", {"shape", false});
load(fullfile(c_fix.data_dir, "case.mat"), "g", "radar_scans", "radar_selected")
fixed = load(fullfile(c_fix.data_dir, "inversion_fix.mat"), "active_field", "active_id");
round0 = load(fullfile(c.data_dir, "cut_levels.mat"), "first_field", "last_field");
inv = c.inversion;

%% Settings ---------------------------------------------------------------
n_sub = 6;                                          % 30 min -> 5 min
frame_rate = 12;
video_width_px = 1920;

%% Active cells: the same as round 0 ----------------------------------------------------
column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
active_mask = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask(:));
assert(isequal(active_id, fixed.active_id) && numel(active_id) == size(round0.first_field, 1), ...
    "The active cells of the fixed and round-0 fields differ.")
n_time = size(fixed.active_field, 2);
cell_air = zeros(g.sz, "single");
cell_air(active_id) = bmm4d.grid.air_fraction(g, active_id);
active_column = any(active_mask, 3);
to_column = @(frame) sum(frame .* cell_air, 3) * g.dz / 1000;     % birds km^-2
dt_seconds = seconds(diff(c.times));

%% Colour scale of the round-0 video (both of its versions) -----------------------------
column_values = [];
for version = ["first_field" "last_field"]
    for i_t = 1:n_time
        frame = zeros(g.sz, "single");
        frame(active_id) = max(round0.(version)(:, i_t), 0);
        column = to_column(frame);
        column_values = [column_values; column(active_column)]; %#ok<AGROW>
    end
end
c_lim = [0 quantile(column_values, 0.995)];
clear column_values

%% Layer winds, background, map ---------------------------------------------------------
n_radar = size(radar_scans, 1);
interpolation_weights = bmm4d.prior.build_horizontal_idw_weights(g, radar_scans, ...
    n_neighbors=inv.velocity_idw_neighbors, power=inv.velocity_idw_power);
uv_profiles = nan(g.sz(3), 2, n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        uv_profiles(:, :, i_r, i_t) = bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g);
    end
end
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
map_bounds = c.figure_bounds;
aspect = diff(map_bounds(3:4)) / (diff(map_bounds(1:2)) * cosd(mean(map_bounds(3:4))));
map_width_px = video_width_px - 230;
video_height_px = 2 * round((map_width_px * aspect + 100) / 2);

%% Video --------------------------------------------------------------------------------
% Full-grid frames; inactive cells take the nearest active value of their
% column, so the advection near terrain steps does not pull in zeros
fields = zeros([g.sz n_time], "single");
for i_t = 1:n_time
    frame = nan(g.sz, "single");
    frame(active_id) = max(fixed.active_field(:, i_t), 0);
    frame = fillmissing(frame, "nearest", 3);
    frame(isnan(frame)) = 0;
    fields(:, :, :, i_t) = frame;
end
video = VideoWriter(fullfile(figure_dir, "us_east_night_" + output_name + ".mp4"), "MPEG-4");
video.FrameRate = frame_rate;
video.Quality = 95;
open(video)
fig = figure(Visible="off", Units="pixels", Position=[0 0 video_width_px video_height_px], Color="w");
middle = round((n_time - 1) * n_sub / 2);
k_frame = 0;
for i_t = 1:n_time - 1
    profiles = 0.5 * (uv_profiles(:, :, :, i_t) + uv_profiles(:, :, :, i_t + 1));
    uv = zeros([g.sz 2], "single");
    for i_c = 1:2
        profile_matrix = reshape(permute(profiles(:, i_c, :), [3 1 2]), n_radar, g.sz(3));
        uv(:, :, :, i_c) = reshape(interpolation_weights * fillmissing(profile_matrix, "constant", 0), g.sz);
    end
    for i_s = 0:n_sub - 1 + (i_t == n_time - 1)
        tau = i_s / n_sub;
        if tau == 0
            column = to_column(fields(:, :, :, i_t));
        elseif tau == 1
            column = to_column(fields(:, :, :, i_t + 1));
        else
            morph = (1 - tau) * advect_layers(fields(:, :, :, i_t), uv, tau * dt_seconds(i_t), g) + ...
                tau * advect_layers(fields(:, :, :, i_t + 1), uv, -(1 - tau) * dt_seconds(i_t), g);
            column = to_column(morph);
        end
        column(~active_column) = nan;
        moment = c.times(i_t) + tau * (c.times(i_t + 1) - c.times(i_t));
        elevation = bmm4d.util.sun_elevation(lat_grid, lon_grid, moment);
        night_alpha = min(1, max(0, (c.sun_limit_deg + 1 - elevation) / 2));
        clf(fig)
        ax = axes(fig, Units="pixels", Position=[60 40 map_width_px video_height_px - 80]);
        bmm4d.plot.plot_night_map(ax, g, background, column, c_lim, states, radar_selected, map_bounds, night_alpha);
        title(ax, "us\_east, " + strrep(fix_folder, "_", "\_") + ", " + string(moment, "dd MMM yyyy HH:mm") + " UTC", ...
            FontName="Helvetica", FontSize=20, FontWeight="bold")
        cb = colorbar(ax, FontName="Helvetica", FontSize=14);
        cb.Label.String = "Birds in the column (km^{-2})";
        image_rgb = print(fig, "-RGBImage", "-r0");
        image_rgb = image_rgb(1:2 * floor(end / 2), 1:2 * floor(size(image_rgb, 2) / 2), :);
        writeVideo(video, image_rgb)
        k_frame = k_frame + 1;
        if k_frame == middle
            imwrite(image_rgb, fullfile(figure_dir, "us_east_night_" + output_name + ".png"))
        end
    end
    fprintf("fix: interval %d/%d done\n", i_t, n_time - 1)
end
close(video)
close(fig)
clear fields morph uv

%% The five spots: round 0 against the fixes ---------------------------------------------
spots = struct("name", {"High Plains", "Appalachians", "Missouri", "Central Texas", "KDMX"}, ...
    "bounds", {[-104.5 -96 38 46], [-88 -79 34 38.5], [-97 -88 35.5 41], [-101.5 -95 27.5 32.8], [-96.5 -91 39.8 43.8]}, ...
    "time", {datetime("2026-10-22 05:30"), datetime("2026-10-22 05:30"), datetime("2026-10-22 05:30"), ...
    datetime("2026-10-22 02:30"), datetime("2026-10-22 06:30")});
versions = {round0.first_field, "round 0"; fixed.active_field, "with the fixes"};
fig = bmm4d.plot.create_export_figure(2 * 900 + 200, numel(spots) * 760, DefaultAxesFontName="Helvetica", ...
    DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, numel(spots), 2, TileSpacing="tight", Padding="tight");
for i_p = 1:numel(spots)
    i_t = find(c.times == spots(i_p).time);
    maps = cell(1, 2);
    for i_v = 1:2
        frame = zeros(g.sz, "single");
        frame(active_id) = max(versions{i_v, 1}(:, i_t), 0);
        maps{i_v} = to_column(frame);
        maps{i_v}(~active_column) = nan;
    end
    inside = lon_grid >= spots(i_p).bounds(1) & lon_grid <= spots(i_p).bounds(2) & ...
        lat_grid >= spots(i_p).bounds(3) & lat_grid <= spots(i_p).bounds(4);
    limit = quantile([maps{1}(inside); maps{2}(inside)], 0.99);
    for i_v = 1:2
        ax = nexttile(layout);
        imagesc(ax, g.lon, g.lat, maps{i_v}', AlphaData=~isnan(maps{i_v}'));
        set(ax, YDir="normal", Color=[0.85 0.92 0.97], FontSize=13)
        axis(ax, spots(i_p).bounds)
        colormap(ax, parula)
        clim(ax, [0 limit])
        hold(ax, "on")
        for i_s = 1:numel(states)
            plot(ax, states(i_s).Lon, states(i_s).Lat, "-", Color=[0.3 0.3 0.3], LineWidth=0.5)
        end
        plot(ax, radar_selected.lon, radar_selected.lat, "k.", MarkerSize=12)
        title(ax, spots(i_p).name + ", " + versions{i_v, 2} + ", " + string(c.times(i_t), "HH:mm") + " UTC", ...
            FontSize=14, FontWeight="bold")
        cb = colorbar(ax);
        cb.FontSize = 12;
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_F_spots_" + output_name + ".png"), Resolution=300);
close(fig)

%% Local functions ----------------------------------------------------------
function moved = advect_layers(frame, uv, shift_seconds, g)
% Move every layer of FRAME along its wind UV by SHIFT_SECONDS: the value at
% x is taken from the departure point x - s u, found in two stages (wind at
% the half-way point). Bilinear interpolation; outside the grid 0.
[x, y] = ndgrid(g.x, g.y);
moved = zeros(size(frame), "like", frame);
for i_z = 1:size(frame, 3)
    if ~any(frame(:, :, i_z), "all")
        continue
    end
    u = double(uv(:, :, i_z, 1));
    v = double(uv(:, :, i_z, 2));
    x_half = x - 0.5 * shift_seconds * u;
    y_half = y - 0.5 * shift_seconds * v;
    u_mid = interpn(g.x, g.y, u, x_half, y_half, "linear", 0);
    v_mid = interpn(g.x, g.y, v, x_half, y_half, "linear", 0);
    moved(:, :, i_z) = interpn(g.x, g.y, double(frame(:, :, i_z)), x - shift_seconds * u_mid, y - shift_seconds * v_mid, "linear", 0);
end
end
