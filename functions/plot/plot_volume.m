function plot_volume(field, varargin)
p = inputParser;
addRequired(p, "field");
addParameter(p, "ax", []);
addParameter(p, "ttl", "");
addParameter(p, "g", []);
addParameter(p, "radar", []);
addParameter(p, "c_lim", []);
addParameter(p, "r_lim", []);
addParameter(p, "z_slice", []);
addParameter(p, "show_colorbar", []);
addParameter(p, "x_label", []);
addParameter(p, "y_label", []);
parse(p, field, varargin{:});

ttl = p.Results.ttl;
g = p.Results.g;
c_lim = p.Results.c_lim;
r_lim = p.Results.r_lim;
show_colorbar = p.Results.show_colorbar;
if isempty(show_colorbar)
    show_colorbar = isempty(c_lim);
end
radar_xy_km = [];

hold on; box on; grid(gca, "on")

if isempty(g)
    sz = size(field);
    x_plot = (1:sz(1))' - sz(1) / 2;
    y_plot = (1:sz(2))' - sz(2) / 2;
    z_plot = (1:sz(3))';
    x_label = "x (km)";
    y_label = "y (km)";
    zlabel("Altitude (m)")
else
    x_plot = g.x(:) - g.x(end) / 2;
    y_plot = g.y(:) - g.y(end) / 2;
    z_plot = g.z(:);
    radar_xy_km = get_plot_radar_xy(g, p.Results.radar);
    x_label = "x";
    y_label = "y";
    zlabel("z")
end

if ~isempty(p.Results.x_label)
    x_label = p.Results.x_label;
end
if ~isempty(p.Results.y_label)
    y_label = p.Results.y_label;
end

field_plot = permute(field, [2 1 3]);
if isempty(p.Results.z_slice)
    z_slice = z_plot(min(4, numel(z_plot)));
else
    z_slice = p.Results.z_slice;
end

[z_plot_slice, field_plot_slice] = extend_top_slice_for_plot(z_plot, field_plot, z_slice);

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end

h = slice(x_plot / 1000, y_plot / 1000, z_plot_slice, field_plot_slice, 0, 0, z_slice);
set(h, "EdgeColor", "none", "FaceAlpha", 0.95)

if ~isempty(radar_xy_km)
    scatter3(radar_xy_km(:, 1), radar_xy_km(:, 2), zeros(size(radar_xy_km, 1), 1), ...
        70, "w", "filled", "MarkerEdgeColor", "k", "LineWidth", 0.6)
end

axis tight
% axis vis3d
view(3)
if ~isempty(c_lim)
    clim(c_lim)
end
if show_colorbar
    cb = colorbar;
    cb.FontSize = 6;
end

if ~isempty(r_lim)
    xlim([-r_lim r_lim])
    ylim([-r_lim r_lim])
    zlim([0 max(z_plot_slice)])
elseif ~isempty(g)
    zlim([0 max(z_plot_slice)])
end

set_axis_ticks_minmax_zero(gca, "x", xlim(gca))
set_axis_ticks_minmax_zero(gca, "y", ylim(gca))
set_axis_ticks_minmax_zero(gca, "z", zlim(gca))

xlabel(x_label)
ylabel(y_label)
set(gca, "FontSize", 6, "XColor", "k", "YColor", "k", "ZColor", "k", "GridColor", "k")
title(ttl, "FontSize", 6, "FontWeight", "normal", "Color", "k")

end

function [z_plot_out, field_plot_out] = extend_top_slice_for_plot(z_plot, field_plot, z_slice)

z_plot_out = z_plot;
field_plot_out = field_plot;

if z_slice > max(z_plot)
    z_plot_out(end + 1) = z_slice;
    field_plot_out(:, :, end + 1) = field_plot(:, :, end);
end
end
