%% Sandbox: videos of the us_east night, round 0 and with the cut levels
% Dev study. The column density of the us_east night (101 radars, 26 frames)
% from SB_us_east_cut_levels (prod/us_east/cut_levels.mat): round 0
% (first_field, no levels) and after the cut-level rounds (last_field), drawn
% as bmm4d.export_maps does (air part of the lowest cell, light land and
% water, cells hidden while the sun is up, fading over +-1 deg) on the full
% grid. Between the 30 min frames, N_SUB frames come from the advection
% morph of SB_morph_animation: both frames moved along the interval's layer
% wind toward the moment and blended, so birds travel with the wind instead
% of cross-fading. Both videos share one colour scale to be compared.
% Output: figures/us_east_night_round0.mp4 and _cut_levels.mp4 (+ middle PNG).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("us_east");
load(fullfile(c.data_dir, "case.mat"), "g", "radar_scans", "radar_selected")
saved = load(fullfile(c.data_dir, "cut_levels.mat"), "first_field", "last_field");
inv = c.inversion;

%% Settings ---------------------------------------------------------------
n_sub = 6;                                          % 30 min -> 5 min
frame_rate = 12;
video_width_px = 1920;
versions = ["first_field" "last_field"];
version_names = ["round0" "cut_levels"];
version_titles = ["round 0 (no cut levels)" "cut levels"];

%% Active cells (as sandbox_load) -----------------------------------------------------
column_mask = [];
if ~isempty(c.domain)
    column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
end
active_mask = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask(:));
assert(numel(active_id) == size(saved.last_field, 1), "The active cells do not match the saved field.")
n_time = size(saved.last_field, 2);
cell_air = zeros(g.sz, "single");
cell_air(active_id) = bmm4d.grid.air_fraction(g, active_id);
active_column = any(active_mask, 3);
to_column = @(frame) sum(frame .* cell_air, 3) * g.dz / 1000;     % birds km^-2
dt_seconds = seconds(diff(c.times));

%% Layer wind profiles at the radars, as the inversion's advection term ------------------
n_radar = size(radar_scans, 1);
interpolation_weights = bmm4d.prior.build_horizontal_idw_weights(g, radar_scans, ...
    n_neighbors=inv.velocity_idw_neighbors, power=inv.velocity_idw_power);
uv_profiles = nan(g.sz(3), 2, n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        uv_profiles(:, :, i_r, i_t) = bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g);
    end
end

%% Shared colour scale and map background ----------------------------------------------
column_values = [];
for i_v = 1:numel(versions)
    field = saved.(versions(i_v));
    for i_t = 1:n_time
        frame = zeros(g.sz, "single");
        frame(active_id) = max(field(:, i_t), 0);
        column = to_column(frame);
        column_values = [column_values; column(active_column)]; %#ok<AGROW>
    end
end
c_lim = [0 quantile(column_values, 0.995)];
clear column_values
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
land_color = [0.96 0.96 0.94];
water_color = [0.85 0.92 0.97];
background = permute(reshape(land_color, 1, 1, 3) .* ~g.mask_water + reshape(water_color, 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
map_bounds = c.figure_bounds;
aspect = diff(map_bounds(3:4)) / (diff(map_bounds(1:2)) * cosd(mean(map_bounds(3:4))));
map_width_px = video_width_px - 230;                % room for tick labels and colorbar
video_height_px = 2 * round((map_width_px * aspect + 100) / 2);

%% Videos ------------------------------------------------------------------------------
middle = round((n_time - 1) * n_sub / 2);
for i_v = 1:numel(versions)
    % Full-grid frames; inactive cells take the nearest active value of
    % their column, so the advection near terrain steps does not pull in zeros
    fields = zeros([g.sz n_time], "single");
    for i_t = 1:n_time
        frame = nan(g.sz, "single");
        frame(active_id) = max(saved.(versions(i_v))(:, i_t), 0);
        frame = fillmissing(frame, "nearest", 3);
        frame(isnan(frame)) = 0;
        fields(:, :, :, i_t) = frame;
    end

    video = VideoWriter(fullfile(figure_dir, "us_east_night_" + version_names(i_v) + ".mp4"), "MPEG-4");
    video.FrameRate = frame_rate;
    video.Quality = 95;
    open(video)
    fig = figure(Visible="off", Units="pixels", Position=[0 0 video_width_px video_height_px], Color="w");
    k_frame = 0;
    for i_t = 1:n_time - 1
        % Interval wind: mean of the two frames' profiles, spread to the grid
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
            title(ax, "us\_east, " + version_titles(i_v) + ", " + string(moment, "dd MMM yyyy HH:mm") + " UTC", ...
                FontName="Helvetica", FontSize=20, FontWeight="bold")
            cb = colorbar(ax, FontName="Helvetica", FontSize=14);
            cb.Label.String = "Birds in the column (km^{-2})";
            image_rgb = print(fig, "-RGBImage", "-r0");
            image_rgb = image_rgb(1:2 * floor(end / 2), 1:2 * floor(size(image_rgb, 2) / 2), :);
            writeVideo(video, image_rgb)
            k_frame = k_frame + 1;
            if k_frame == middle
                imwrite(image_rgb, fullfile(figure_dir, "us_east_night_" + version_names(i_v) + ".png"))
            end
        end
        fprintf("%s: interval %d/%d done\n", version_names(i_v), i_t, n_time - 1)
    end
    close(video)
    close(fig)
end

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
