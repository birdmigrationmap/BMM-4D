%% Sandbox radars on 2 nights: production run, video and sections
% The sandbox night (21 October 2026, autumn) and one dry spring night of
% SB_night_screen (12 May 2021), each as case night_YYYYMMDD (examples/case_settings: the sandbox
% radars and grid, every 30 min from 1 h after sunset to 30 min before
% sunrise), run through the production pipeline with its defaults (gate
% estimator, no cuts, lambda_profile = 1). A night is skipped when its
% inversion exists for the same frames.
% For each night, a video (one frame per time) and its middle frame as PNG:
% - column density, and column / 150 km regional mean (radar imprint);
% - radar-distance curve of the frame (1 km bins) against all frames;
% - west-east sections through KDMX (41.7 N) and KLSX (38.7 N), south-north
%   sections through KDMX (93.7 W) and KILX (89.3 W);
% - mean profile near the radars (< 25 km) and between them (100-150 km).
% Outputs: <figure_dir>/night_YYYYMMDD/night_review.mp4, .png, every frame
% in night_review_frames/, and the
% production maps (export_maps); scores in studies/sandbox/figures/nights.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
night_dates = [datetime(2026, 10, 21); datetime(2021, 5, 12)];
scores = table();

for i_n = 1:numel(night_dates)
    c = case_settings("night_" + string(night_dates(i_n), "yyyyMMdd"));
    fprintf("[nights] %s: %s to %s\n", c.name, string(c.times(1)), string(c.times(end)))

    %% Production pipeline (stages already done are skipped)
    manifest_file = fullfile(c.data_dir, "manifest.mat");
    done = isfile(fullfile(c.data_dir, "inversion.mat")) && isfile(manifest_file);
    if done
        manifest = load(manifest_file, "times");
        done = isequal(manifest.times, c.times);
    end
    if ~done
        bmm4d.get_volumes(c);
        bmm4d.clean_volumes(c);
        bmm4d.static_clutter(c);
        bmm4d.build_observations(c);
        bmm4d.invert(c);
        bmm4d.export_maps(c);
    end

    %% Review video
    scores = [scores; night_review(c)]; %#ok<AGROW>
    writetable(scores, fullfile(figure_dir, "nights.csv"))
end
disp(scores)

%% Local functions ----------------------------------------------------------
function score = night_review(c)
% Video of one night (see header) and its scores: radar-distance bias
% within 200 km (bmm4d.plot.radar_distance_bias) and the column below 1 km
% near and between radars.
load(fullfile(c.data_dir, "case.mat"), "g", "radar_selected")
load(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id")
n_time = size(active_field, 2);
[ix, iy, iz] = ind2sub(g.sz, active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
cell_air = bmm4d.grid.air_fraction(g, active_id);   % columns count the air part of the lowest cell
is_active = false(g.sz(1:2));
is_active(column_id) = true;
z_layer = squeeze(g.z);

% Distance to the nearest radar, regional mean weights
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for k = 1:height(radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(radar_selected.lat(k)) + ...
        cosd(lat_grid) .* cosd(radar_selected.lat(k)) .* cosd(lon_grid - radar_selected.lon(k)))));
end
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
fine_edges = 0:1:250;
fine_center = fine_edges(1:end - 1)' + 0.5;
fine_bin = discretize(distance_km, fine_edges);
cell_distance = distance_km(column_id);
near = cell_distance < 25;
between = cell_distance >= 100 & cell_distance < 150;

% Per frame: column, relative map, curve, profiles
column = nan([g.sz(1:2) n_time]);
relative = nan([g.sz(1:2) n_time]);
curve = nan(numel(fine_center), n_time);
profile_near = zeros(g.sz(3), n_time);
profile_between = zeros(g.sz(3), n_time);
for i_t = 1:n_time
    field = double(active_field(:, i_t));
    map = reshape(accumarray(column_id, field .* cell_air, [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    map = max(map, 0) .* is_active;
    regional = imgaussfilt(map, sigma_cells) ./ active_weight;
    ok = is_active & regional > 5;
    ratio = nan(g.sz(1:2));
    ratio(ok) = map(ok) ./ regional(ok);
    use = ok & ~isnan(fine_bin);
    curve(:, i_t) = accumarray(fine_bin(use), ratio(use), [numel(fine_center) 1], @mean, nan);
    map(~is_active) = nan;
    column(:, :, i_t) = map;
    relative(:, :, i_t) = ratio;
    profile_near(:, i_t) = accumarray(iz(near), field(near), [g.sz(3) 1], @mean, nan);
    profile_between(:, i_t) = accumarray(iz(between), field(between), [g.sz(3) 1], @mean, nan);
end
[~, ~, bias_frame] = bmm4d.plot.radar_distance_bias(active_field, active_id, g, radar_selected);
low = z_layer < 1000;
score = table(string(c.name), numel(c.times), mean(bias_frame, "omitnan"), ...
    mean(sum(profile_near(low, :), 1) * g.dz / 1000), mean(sum(profile_between(low, :), 1) * g.dz / 1000), ...
    max(column, [], "all"), VariableNames=["night" "frames" "bias" "below_1km_near" "below_1km_between" "max_column"]);

% Video
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
column_limit = max(50, round(prctile(column(:), 99), -1));
% Sections: west-east through KDMX and KLSX, north-south through KDMX and KILX
sections = struct("radar", {"KDMX", "KDMX", "KLSX", "KILX"}, "direction", {"we", "ns", "we", "ns"});
out_dir = c.figure_dir;
frame_dir = fullfile(out_dir, "night_review_frames");
if ~isfolder(frame_dir)
    mkdir(frame_dir);
end
video = VideoWriter(fullfile(out_dir, "night_review.mp4"), "MPEG-4");
video.FrameRate = 1;
open(video)
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, 4, 3, TileSpacing="tight", Padding="tight");
    title(t, sprintf("%s | %s UTC | frame bias %.3f (night %.3f)", c.name, string(c.times(i_t), "yyyy-MM-dd HH:mm"), ...
        bias_frame(i_t), score.bias), FontSize=15, Interpreter="none")
    ax = nexttile(t, 1, [2 1]);
    bmm4d.plot.plot_night_map(ax, g, background, column(:, :, i_t), [0 column_limit], states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
    for s_k = sections
        k = radar_selected.name == s_k.radar;
        if s_k.direction == "we"
            yline(ax, radar_selected.lat(k), "w--")
        else
            xline(ax, radar_selected.lon(k), "w:", LineWidth=1.5)
        end
    end
    cb = colorbar(ax); cb.Label.String = "Column density (birds km^{-2})";
    title(ax, "Reconstruction (lines: sections)", FontWeight="normal")
    ax = nexttile(t, 2, [2 1]);
    bmm4d.plot.plot_night_map(ax, g, background, log2(relative(:, :, i_t)), [-1 1], states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, signed_map); cb = colorbar(ax); cb.Label.String = "log_2(column / 150 km regional mean)";
    title(ax, "Radar imprint", FontWeight="normal")
    ax = nexttile(t, 3, [2 1]); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, fine_center, curve, Color=[0.8 0.8 0.8])
    plot(ax, fine_center, curve(:, i_t), "k-", LineWidth=1.5)
    yline(ax, 1, "k--"); xlim(ax, [0 250]); ylim(ax, [0.4 1.6]); set(ax, FontSize=12)
    xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column / regional mean")
    title(ax, "Radar-distance curve (grey: other frames)", FontWeight="normal")
    density_limit = max(20, round(prctile(active_field(:, i_t), 99.5), -1));
    field = nan(g.sz);
    field(active_id) = active_field(:, i_t);
    section_tiles = [7 8 10 11];
    for i_s = 1:numel(sections)
        k = find(radar_selected.name == sections(i_s).radar, 1);
        ax = nexttile(t, section_tiles(i_s));
        hold(ax, "on")
        if sections(i_s).direction == "we"
            [~, line_id] = min(abs(g.lat - radar_selected.lat(k)));
            section = squeeze(field(:, line_id, :))';
            imagesc(ax, g.lon, z_layer, section, AlphaData=~isnan(section))
            on_line = abs(radar_selected.lat - g.lat(line_id)) < 0.6;
            plot(ax, radar_selected.lon(on_line), 250 + 0 * radar_selected.lon(on_line), "kv", MarkerFaceColor="k")
            xlim(ax, c.figure_bounds(1:2))
            section_title = sprintf("W-E through %s (%.1f°N)", sections(i_s).radar, g.lat(line_id));
        else
            [~, line_id] = min(abs(g.lon - radar_selected.lon(k)));
            section = squeeze(field(line_id, :, :))';
            imagesc(ax, g.lat, z_layer, section, AlphaData=~isnan(section))
            on_line = abs(radar_selected.lon - g.lon(line_id)) < 0.8;
            plot(ax, radar_selected.lat(on_line), 250 + 0 * radar_selected.lat(on_line), "kv", MarkerFaceColor="k")
            xlim(ax, c.figure_bounds(3:4))
            section_title = sprintf("S-N through %s (%.1f°W)", sections(i_s).radar, -g.lon(line_id));
        end
        set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=11); clim(ax, [0 density_limit]); colormap(ax, parula)
        ylim(ax, [0 3000])
        if i_s == 2 || i_s == 4
            cb = colorbar(ax); cb.Label.String = "Density (birds km^{-3})";
        end
        if i_s == 1 || i_s == 3
            ylabel(ax, "Height (m a.s.l.)")
        end
        title(ax, section_title + ", ▼ radars", FontWeight="normal")
    end
    ax = nexttile(t, 9, [2 1]); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, profile_near(:, i_t), z_layer, "-", Color=[0.1 0.3 0.6], LineWidth=2, DisplayName="< 25 km from a radar")
    plot(ax, profile_between(:, i_t), z_layer, "-", Color=[0.85 0.33 0.1], LineWidth=2, DisplayName="100-150 km (between radars)")
    ylim(ax, [0 3000]); set(ax, FontSize=12)
    xlabel(ax, "Mean density (birds km^{-3})"); ylabel(ax, "Height (m a.s.l.)")
    legend(ax, Location="northeast"); title(ax, "Mean profile", FontWeight="normal")
    drawnow
    writeVideo(video, getframe(fig))
    bmm4d.plot.export_fixed_figure(fig, fullfile(frame_dir, sprintf("frame_%02d_%s.png", i_t, string(c.times(i_t), "HHmm"))), Resolution=100)
    if i_t == ceil(n_time / 2)
        bmm4d.plot.export_fixed_figure(fig, fullfile(out_dir, "night_review.png"), Resolution=150)
    end
end
close(video)
close(fig)
end
