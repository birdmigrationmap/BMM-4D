function [density_image, time_text] = plot_satellite_density(ax, field, g, radar, background, c_lim, fade_density, time_label)
% Plot a density overlay using one common satellite view and opacity ramp.

hold(ax, "on")
image(ax, background.lon, background.lat, background.image)
alpha_data = min(max((field - fade_density) ./ (5000 .* g.dz ./ 1000), 0), 1);
% Transpose [x y] fields only for latitude-by-longitude image display.
density_image = imagesc(ax, g.lon, g.lat, field', AlphaData=alpha_data');
scatter(ax, radar.lon, radar.lat, 45, "w", "filled", MarkerEdgeColor="k", LineWidth=0.6)
colormap(ax, parula(256))
clim(ax, c_lim)
set(ax, FontName="Helvetica", YDir="normal", Color=[0.88 0.94 0.97], XColor="none", YColor="none")
axis(ax, "equal")
xlim(ax, background.lonlim)
ylim(ax, background.latlim)
time_text = text(ax, background.lonlim(1) + 0.15, background.latlim(2) - 0.15, time_label, ...
    FontName="Helvetica", FontSize=18, FontWeight="bold", Color="w", ...
    BackgroundColor=[0.15 0.15 0.15], Margin=3, VerticalAlignment="top", Interpreter="none");
end
