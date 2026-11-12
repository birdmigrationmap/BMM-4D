%% Sandbox: smoother animations by advecting between frames
% Dev study. Between two reconstructed frames t and t+1 (dt apart), the frame
% at fraction tau of the interval is
%   rho(tau) = (1 - tau) * A(rho_t, +tau dt) + tau * A(rho_t+1, -(1 - tau) dt),
% where A(rho, s) moves each layer of rho along the layer wind by s seconds
% (semi-Lagrangian: the value at x is taken from x - s u, the departure point
% found in two stages for curving winds; bilinear interpolation). Both
% frames are moved toward tau, so the result equals the frames at tau = 0 and
% 1, structures travel with the wind instead of fading out in one place and in
% at another (cross-fade), and change that the wind does not explain (takeoff,
% landing) fades in. The wind is the interval's mean layer wind of the radar
% profiles, interpolated as in the inversion's advection term
% (bmm4d.advect.build_radar_advection_steps). Cost: two interpolations per
% layer and intermediate frame.
% Output: figures/morph_<case>.mp4 and .png, cross-fade (top) against
% advection morph (bottom), column map and a west-east section.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("case_name", "var")
    case_name = "night_20261021";
end
c = case_settings(case_name);
load(fullfile(c.data_dir, "case.mat"), "g", "radar_scans", "radar_selected")
load(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id")
inv = c.inversion;

%% Settings ---------------------------------------------------------------
n_sub = 6;                                          % frames per interval: 30 min -> 5 min
section_radar = "KLSX";
dt_seconds = seconds(diff(c.times));

%% Layer wind per interval, as in the inversion ------------------------------------------
n_radar = size(radar_scans, 1);
n_time = size(radar_scans, 2);
interpolation_weights = bmm4d.prior.build_horizontal_idw_weights(g, radar_scans, ...
    n_neighbors=inv.velocity_idw_neighbors, power=inv.velocity_idw_power);
uv_profiles = nan(g.sz(3), 2, n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        uv_profiles(:, :, i_r, i_t) = bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g);
    end
end
uv = cell(n_time - 1, 1);                           % [nx ny nz 2], m/s, east and north
for i_t = 1:n_time - 1
    profiles = 0.5 * (uv_profiles(:, :, :, i_t) + uv_profiles(:, :, :, i_t + 1));
    uv{i_t} = zeros([g.sz 2]);
    for i_c = 1:2
        profile_matrix = reshape(permute(profiles(:, i_c, :), [3 1 2]), n_radar, g.sz(3));
        uv{i_t}(:, :, :, i_c) = reshape(interpolation_weights * fillmissing(profile_matrix, "constant", 0), g.sz);
    end
end

%% Full-grid fields ---------------------------------------------------------
fields = zeros([g.sz n_time], "single");
for i_t = 1:n_time
    frame = zeros(g.sz, "single");
    frame(active_id) = max(active_field(:, i_t), 0);
    fields(:, :, :, i_t) = frame;
end
active = false(g.sz);
active(active_id) = true;
% Cells below ground (and other inactive cells) take the nearest active value of
% their column, so that advection near terrain steps does not pull in zeros
filled = fields;
for i_t = 1:n_time
    frame = fields(:, :, :, i_t);
    frame(~active) = nan;
    frame = fillmissing(frame, "nearest", 3);
    frame(isnan(frame)) = 0;
    filled(:, :, :, i_t) = frame;
end

%% Animation: cross-fade against advection morph -------------------------------------------
air = bmm4d.grid.air_fraction(g, active_id);
air_grid = zeros(g.sz);
air_grid(active_id) = air;
to_column = @(frame) sum(double(frame) .* air_grid, 3) * g.dz / 1000;      % birds km^-2
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
active_column = any(active, 3);
i_radar = find(radar_selected.name == section_radar, 1);
[~, y_id] = min(abs(g.lat - radar_selected.lat(i_radar)));
distance_km = (g.lon - radar_selected.lon(i_radar)) * 111.32 * cosd(radar_selected.lat(i_radar));
section_limit = prctile(fields(:, y_id, :, :), 99, "all");
video = VideoWriter(fullfile(figure_dir, "morph_" + case_name + ".mp4"), "MPEG-4");
video.FrameRate = 12;
open(video)
fig = bmm4d.plot.create_export_figure(2400, 1600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
morph_seconds = 0;
for i_t = 1:n_time - 1
    for i_s = 0:n_sub - 1
        tau = i_s / n_sub;
        crossfade = (1 - tau) * fields(:, :, :, i_t) + tau * fields(:, :, :, i_t + 1);
        tic
        morph = (1 - tau) * advect_layers(filled(:, :, :, i_t), uv{i_t}, tau * dt_seconds(i_t), g) + ...
            tau * advect_layers(filled(:, :, :, i_t + 1), uv{i_t}, -(1 - tau) * dt_seconds(i_t), g);
        morph_seconds = morph_seconds + toc;
        morph(~active) = 0;
        clf(fig)
        t = tiledlayout(fig, 2, 2, TileSpacing="tight", Padding="tight");
        title(t, string(c.times(i_t) + tau * (c.times(i_t + 1) - c.times(i_t)), "dd MMM yyyy HH:mm") + " UTC", FontSize=16, FontWeight="bold")
        frames = {crossfade, morph};
        names = ["cross-fade" "advection morph"];
        for k = 1:2
            ax = nexttile(t);
            column = nan(g.sz(1:2));
            values = to_column(frames{k});
            column(active_column) = log10(max(values(active_column), 1));
            bmm4d.plot.plot_night_map(ax, g, background, column, log10([5 300]), states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
            colormap(ax, parula)
            hold(ax, "on")
            plot(ax, radar_selected.lon(i_radar) + [-3 3], g.lat(y_id) * [1 1], "k-", LineWidth=1.2)
            ylabel(ax, names(k), FontWeight="bold", FontSize=14)
            ax = nexttile(t);
            section = squeeze(frames{k}(:, y_id, :))';
            section(~squeeze(active(:, y_id, :))') = nan;
            imagesc(ax, distance_km, squeeze(g.z), section, AlphaData=~isnan(section))
            set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12, XLim=[-200 200], YLim=[0 3000])
            clim(ax, [0 section_limit])
            colormap(ax, parula)
            title(ax, section_radar + " west-east", FontSize=13)
            ylabel(ax, "Height (m a.s.l.)")
        end
        drawnow
        writeVideo(video, getframe(fig))
        if i_t == round(n_time / 2) && i_s == n_sub / 2
            bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "morph_" + case_name + ".png"), Resolution=150)
        end
    end
end
close(video)
close(fig)
fprintf("Morph: %.3f s per intermediate frame (%d frames)\n", morph_seconds / ((n_time - 1) * n_sub), (n_time - 1) * n_sub)

%% Local functions ----------------------------------------------------------
function moved = advect_layers(frame, uv, shift_seconds, g)
% Move every layer of FRAME along its wind UV by SHIFT_SECONDS: the value at
% x is taken from the departure point x - s u, found in two stages (wind at
% the half-way point) so that curving winds are followed. Bilinear
% interpolation; outside the grid the density is 0.
[x, y] = ndgrid(g.x, g.y);
moved = zeros(size(frame), "like", frame);
for i_z = 1:size(frame, 3)
    u = uv(:, :, i_z, 1);
    v = uv(:, :, i_z, 2);
    u_mid = interpn(g.x, g.y, u, x - 0.5 * shift_seconds * u, y - 0.5 * shift_seconds * v, "linear", 0);
    v_mid = interpn(g.x, g.y, v, x - 0.5 * shift_seconds * u, y - 0.5 * shift_seconds * v, "linear", 0);
    moved(:, :, i_z) = interpn(g.x, g.y, double(frame(:, :, i_z)), x - shift_seconds * u_mid, y - shift_seconds * v_mid, "linear", 0);
end
end
