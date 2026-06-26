function plot_xy(field, varargin)

p = inputParser;
addRequired(p, "field");
addParameter(p, "ax", []);
addParameter(p, "ttl", "");
addParameter(p, "grid", []);
addParameter(p, "radar", []);
addParameter(p, "c_lim", []);
addParameter(p, "r_lim", []);
addParameter(p, "show_colorbar", []);
addParameter(p, "radar_marker_size", 70);
addParameter(p, "x_label", "x (km)");
addParameter(p, "y_label", "y (km)");
addParameter(p, "use_geographic_axis", false);
parse(p, field, varargin{:});

ttl = p.Results.ttl;
g = p.Results.grid;
c_lim = p.Results.c_lim;
r_lim = p.Results.r_lim;
show_colorbar = p.Results.show_colorbar;
use_geographic_axis = p.Results.use_geographic_axis;
if isempty(show_colorbar)
    show_colorbar = isempty(c_lim);
end
radar_xy = [];

if isempty(g)
    x_plot = (1:size(field, 1))' - size(field, 1) / 2;
    y_plot = (1:size(field, 2))' - size(field, 2) / 2;
else
    if use_geographic_axis && isfield(g, "lon") && isfield(g, "lat")
        x_plot = g.lon(:);
        y_plot = g.lat(:);
        if ~isempty(p.Results.radar)
            radar_xy = [p.Results.radar.lon p.Results.radar.lat];
        end
    else
        x_plot = g.x(:) - g.x(end) / 2;
        y_plot = g.y(:) - g.y(end) / 2;
        radar_xy = get_plot_radar_xy(g, p.Results.radar);
    end
end

field_xy = field';

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end
hold on; box on; grid on

if use_geographic_axis && ~isempty(g)
    imagesc(x_plot, y_plot, field_xy)
else
    imagesc(x_plot / 1000, y_plot / 1000, field_xy)
end
if ~isempty(radar_xy)
    scatter(radar_xy(:, 1), radar_xy(:, 2), p.Results.radar_marker_size, ...
        "w", "filled", "MarkerEdgeColor", "k", "LineWidth", 0.6)
end

set(gca, "YDir", "normal")
axis image
if ~isempty(c_lim)
    clim(c_lim)
end
if show_colorbar
    cb = colorbar;
    cb.FontSize = 6;
end

if use_geographic_axis && ~isempty(g)
    xlim([min(x_plot) max(x_plot)])
    ylim([min(y_plot) max(y_plot)])
elseif ~isempty(r_lim)
    xlim([-r_lim r_lim])
    ylim([-r_lim r_lim])
elseif ~isempty(g)
    r_lim = max(abs([x_plot(:); y_plot(:)])) / 1000;
    xlim([-r_lim r_lim])
    ylim([-r_lim r_lim])
end

if ~use_geographic_axis || isempty(g)
    set_axis_ticks_minmax_zero(gca, "x", xlim(gca))
    set_axis_ticks_minmax_zero(gca, "y", ylim(gca))
end

xlabel(p.Results.x_label)
ylabel(p.Results.y_label)
set(gca, "FontSize", 6, "XColor", "k", "YColor", "k", "GridColor", "k")
xlabel(p.Results.x_label, "Color", "k")
ylabel(p.Results.y_label, "Color", "k")
title(ttl, "FontSize", 6, "FontWeight", "normal", "Color", "k")
end
