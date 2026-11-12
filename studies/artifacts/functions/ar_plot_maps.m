function ar_plot_maps(S, results, frames, file, varargin)
%AR_PLOT_MAPS Column maps of several solves side by side (dev).
%   AR_PLOT_MAPS(S, RESULTS, FRAMES, FILE) draws one row per result
%   (ar_solve) and one column per frame index in FRAMES, with a shared
%   colour scale per frame (99.5th percentile over the results), state
%   lines and radars, and exports FILE (PNG, 300 dpi). Options:
%   - "bounds": [lon_min lon_max lat_min lat_max] of the maps (default
%     the grid);
%   - "extra": [nx ny] map drawn as an extra last column (e.g. the DEM);
%   - "extra_title", "extra_limits".

p = inputParser;
addParameter(p, "bounds", [S.g.lon(1) S.g.lon(end) S.g.lat(1) S.g.lat(end)]);
addParameter(p, "extra", []);
addParameter(p, "extra_title", "");
addParameter(p, "extra_limits", []);
parse(p, varargin{:});
opt = p.Results;

n_row = numel(results);
n_col = numel(frames) + ~isempty(opt.extra);
panel_px = 700;
aspect = diff(opt.bounds(3:4)) / (diff(opt.bounds(1:2)) * cosd(mean(opt.bounds(3:4))));
fig = bmm4d.plot.create_export_figure(panel_px * n_col + 200, round(panel_px * aspect * n_row) + 150, ...
    DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, n_row, n_col, TileSpacing="tight", Padding="tight");
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
for i_row = 1:n_row
    for i_col = 1:n_col
        ax = nexttile(layout);
        if i_col <= numel(frames)
            i_t = frames(i_col);
            values = cellfun(@(r) r.column(:, :, i_t), results, UniformOutput=false);
            limit = quantile(vertcat(values{:}), 0.995, "all");
            map = results{i_row}.column(:, :, i_t);
            imagesc(ax, S.g.lon, S.g.lat, map', AlphaData=~isnan(map'));
            colormap(ax, parula)
            clim(ax, [0 limit])
            title_text = [results{i_row}.name; string(S.c.times(i_t), "HH:mm") + " UTC"];
        else
            imagesc(ax, S.g.lon, S.g.lat, opt.extra');
            colormap(ax, gray)
            if ~isempty(opt.extra_limits)
                clim(ax, opt.extra_limits)
            end
            title_text = opt.extra_title;
        end
        set(ax, YDir="normal", Color=[0.85 0.92 0.97], FontSize=13)
        axis(ax, opt.bounds)
        hold(ax, "on")
        for i_s = 1:numel(states)
            plot(ax, states(i_s).Lon, states(i_s).Lat, "-", Color=[0.3 0.3 0.3], LineWidth=0.5)
        end
        plot(ax, S.radar_selected.lon, S.radar_selected.lat, "k.", MarkerSize=12)
        title(ax, title_text, FontSize=14, FontWeight="bold")
        cb = colorbar(ax);
        cb.FontSize = 12;
    end
end
bmm4d.plot.export_fixed_figure(fig, file, Resolution=300);
close(fig)
end
