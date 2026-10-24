function plot_night_map(ax, g, background, column_density, c_lim, states, radar_selected, map_bounds, alpha)
%PLOT_NIGHT_MAP Draw one column-density map in lon/lat with state outlines.
%   COLUMN_DENSITY is [x y] on g; it is transposed here because image()
%   expects rows along y. NaN cells (outside the active domain) show the
%   land/water background. ALPHA ([x y], 0-1) fades the density, e.g.
%   beyond the range of the radar data.

hold(ax, "on")
image(ax, g.lon, g.lat, background)
imagesc(ax, g.lon, g.lat, column_density', AlphaData=alpha' .* ~isnan(column_density'))
for i_state = 1:numel(states)
    plot(ax, states(i_state).Lon, states(i_state).Lat, Color=[0.35 0.35 0.35], LineWidth=0.4)
end
plot(ax, radar_selected.lon, radar_selected.lat, "k.", MarkerSize=6)
colormap(ax, parula(256))
clim(ax, c_lim)
axis(ax, map_bounds)
set(ax, YDir="normal", DataAspectRatio=[1 cosd(mean(map_bounds(3:4))) 1], ...
    FontName="Helvetica", FontSize=13, Layer="top", Box="on")
end
