function skin = ar_check_figures(S, results, section_lats, section_lon, section_time, file_prefix)
%AR_CHECK_FIGURES Maps and sections to look for artefacts of a vertical rule (dev).
%   SKIN = AR_CHECK_FIGURES(S, RESULTS, SECTION_LATS, SECTION_LON,
%   SECTION_TIME, FILE_PREFIX) for the solves RESULTS (ar_solve) of case S:
%   - <prefix>_maps.png: column maps at the first, middle and last frame
%     (ar_plot_maps), the last row the ratio of the last result to the
%     first (log2 scale);
%   - <prefix>_sections.png: west-east sections at SECTION_LATS, between
%     SECTION_LON, at SECTION_TIME, one column per latitude and one row per
%     result, density by height above sea level, the 30" ground in gray
%     and the grid ground in red;
%   - SKIN: per result, the median ratio of the lowest active cell to the
%     cell 300 m above it, over ground above 1000 m and below 450 m (a
%     layer piled into the lowest cell over the ridges shows as a ratio
%     well above the valleys').
% imagesc wants [y x] / [z x], so map and section slices are transposed in
% the plotting calls.

g = S.g;
z = squeeze(g.z);
ground = max(g.dem, 0);
font = 13;
n_time = numel(S.c.times);
frames = [1 ceil(n_time / 2) n_time];

% Column maps, then the ratio of the last result to the first
ratio_result = struct("name", "Column ratio to " + results{1}.name, ...
    "column", log2(max(results{end}.column, 1) ./ max(results{1}.column, 1)));
ar_plot_maps(S, results, frames, file_prefix + "_maps.png", ...
    bounds=[S.g.lon(1) S.g.lon(end) S.g.lat(1) S.g.lat(end)], extra=ground, extra_title="Ground (m)", extra_limits=[0 1600]);
fig = bmm4d.plot.create_export_figure(3600, 900, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 1, numel(frames), TileSpacing="compact", Padding="compact");
for i_f = 1:numel(frames)
    ax = nexttile(layout);
    map = ratio_result.column(:, :, frames(i_f));
    imagesc(ax, g.lon, g.lat, map', AlphaData=~isnan(map'));
    hold(ax, "on")
    contour(ax, g.lon, g.lat, ground', [1000 1000], "k", LineWidth=0.8)
    set(ax, YDir="normal", Color="w", FontSize=font, DataAspectRatio=[1 cosd(mean(g.lat)) 1])
    colormap(ax, interp1([-1 0 1], [0.2 0.4 0.8; 1 1 1; 0.8 0.3 0.2], linspace(-1, 1, 256)))
    clim(ax, [-1 1])
    title(ax, ratio_result.name + ", " + string(S.c.times(frames(i_f)), "HH:mm"), FontSize=font, FontWeight="bold")
end
cb = colorbar(ax, Color="k");
cb.Label.String = "log_2 column ratio (black: ground 1000 m)";
bmm4d.plot.export_fixed_figure(fig, file_prefix + "_ratio.png", Resolution=300);
close(fig)

% Sections
dem_data = load(bmm4d.util.data_path("reference", "DEM", "dem.mat"), "dem", "dem_lat", "dem_lon");
lon_id = find(dem_data.dem_lon >= section_lon(1) & dem_data.dem_lon <= section_lon(2));
x_id = find(g.lon >= section_lon(1) & g.lon <= section_lon(2));
i_t = find(S.c.times == section_time);
fields = cell(numel(results), numel(section_lats));
for i_r = 1:numel(results)
    frame = nan(g.sz);
    frame(S.active_id) = results{i_r}.field(:, i_t);
    for i_l = 1:numel(section_lats)
        [~, y_id] = min(abs(g.lat - section_lats(i_l)));
        fields{i_r, i_l} = squeeze(frame(x_id, y_id, :));
    end
end
limit = quantile(cell2mat(cellfun(@(f) f(:), fields(:), UniformOutput=false)), 0.99);
fig = bmm4d.plot.create_export_figure(3600, 450 * numel(results) + 200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, numel(results), numel(section_lats), TileSpacing="compact", Padding="compact");
for i_r = 1:numel(results)
    for i_l = 1:numel(section_lats)
        ax = nexttile(layout);
        [~, y_id] = min(abs(g.lat - section_lats(i_l)));
        [~, dem_row] = min(abs(dem_data.dem_lat - section_lats(i_l)));
        imagesc(ax, g.lon(x_id), z / 1000, fields{i_r, i_l}', AlphaData=~isnan(fields{i_r, i_l}'));
        hold(ax, "on")
        dem_profile = max(double(dem_data.dem(dem_row, lon_id)), 0);
        area(ax, dem_data.dem_lon(lon_id), dem_profile / 1000, FaceColor=[0.55 0.55 0.55], EdgeColor="none")
        plot(ax, g.lon(x_id), ground(x_id, y_id) / 1000, Color=[0.75 0.1 0.1], LineWidth=1)
        set(ax, YDir="normal", Color="w", FontSize=font - 2, Layer="top")
        colormap(ax, parula)
        clim(ax, [0 limit])
        ylim(ax, [0 3])
        xlim(ax, section_lon)
        if i_r == 1
            title(ax, sprintf("%.1f° N", section_lats(i_l)), FontSize=font, FontWeight="bold")
        end
        if i_l == 1
            ylabel(ax, [results{i_r}.name; "Height (km asl)"], FontSize=font - 1)
        end
    end
end
cb = colorbar(ax, Color="k");
cb.Layout.Tile = "east";
cb.Label.String = "Birds km^{-3}, " + string(section_time, "HH:mm") + " UTC";
bmm4d.plot.export_fixed_figure(fig, file_prefix + "_sections.png", Resolution=300);
close(fig)

% Lowest cell against the cell 300 m above it
[has_active, lowest] = max(S.active_mask_3d, [], 3);
nxy = prod(g.sz(1:2));
columns = find(has_active & lowest + 3 <= g.sz(3));
low_id = columns + (lowest(columns) - 1) * nxy;
up_id = columns + (lowest(columns) + 2) * nxy;
position = zeros(prod(g.sz), 1);
position(S.active_id) = 1:numel(S.active_id);
keep = position(low_id) > 0 & position(up_id) > 0;
columns = columns(keep);
low_id = low_id(keep);
up_id = up_id(keep);
ridge = ground(columns) > 1000;
valley = ground(columns) < 450;
skin = table();
for i_r = 1:numel(results)
    ratio = mean(max(results{i_r}.field(position(low_id), :), 0.1), 2) ./ mean(max(results{i_r}.field(position(up_id), :), 0.1), 2);
    skin = [skin; table(results{i_r}.name, median(ratio(ridge)), median(ratio(valley)), ...
        VariableNames=["solve" "lowest_over_300m_ridge" "lowest_over_300m_valley"])]; %#ok<AGROW>
end
end
