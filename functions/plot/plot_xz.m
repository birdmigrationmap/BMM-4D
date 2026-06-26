function plot_xz(field, varargin)

field_size_raw = size(field);
field = squeeze(field);

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
addParameter(p, "x_label", []);
addParameter(p, "y_label", []);
addParameter(p, "use_geographic_axis", false);
parse(p, field, varargin{:});

ttl = p.Results.ttl;
g = p.Results.grid;
c_lim = p.Results.c_lim;
r_lim = p.Results.r_lim;
horizontal_label = "x (km)";
y_label = "z (m)";
show_colorbar = p.Results.show_colorbar;
use_geographic_axis = p.Results.use_geographic_axis;
if isempty(show_colorbar)
    show_colorbar = isempty(c_lim);
end
radar_x = [];

if isempty(g)
    x_plot = (1:size(field, 1))' - size(field, 1) / 2;
    z_plot = (1:size(field, 2))';
    z_max = max(z_plot);
else
    [slice_axis, x_plot] = infer_vertical_slice_axis(field_size_raw, size(field), g);
    z_plot = g.z(:);
    if numel(z_plot) > 1
        z_max = z_plot(end) + median(diff(z_plot)) / 2;
    else
        z_max = z_plot(end);
    end
    if use_geographic_axis && isfield(g, "lon") && isfield(g, "lat")
        if slice_axis == "y"
            x_plot = g.lat(:);
            horizontal_label = "latitude (deg)";
            if ~isempty(p.Results.radar)
                radar_x = p.Results.radar.lat;
            end
        else
            x_plot = g.lon(:);
            horizontal_label = "longitude (deg)";
            if ~isempty(p.Results.radar)
                radar_x = p.Results.radar.lon;
            end
        end
    else
        radar_xy_km = get_plot_radar_xy(g, p.Results.radar);
        if slice_axis == "y"
            radar_x = radar_xy_km(:, 2);
            horizontal_label = "y (km)";
        else
            radar_x = radar_xy_km(:, 1);
        end
    end
end

field_xz = field';

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end
hold on; box on; grid on

if use_geographic_axis && ~isempty(g)
    imagesc(x_plot, z_plot, field_xz)
else
    imagesc(x_plot / 1000, z_plot, field_xz)
end
if ~isempty(radar_x)
    scatter(radar_x, zeros(size(radar_x)), p.Results.radar_marker_size, ...
        "w", "filled", "MarkerEdgeColor", "k", "LineWidth", 0.6)
end

set(gca, "YDir", "normal")
axis tight
if ~isempty(c_lim)
    clim(c_lim)
end
if show_colorbar
    cb = colorbar;
    cb.FontSize = 6;
end

if use_geographic_axis && ~isempty(g)
    xlim([min(x_plot) max(x_plot)])
elseif ~isempty(r_lim)
    xlim([-r_lim r_lim])
elseif ~isempty(g)
    r_lim = max(abs(x_plot)) / 1000;
    xlim([-r_lim r_lim])
end
ylim([0 z_max])

if ~use_geographic_axis || isempty(g)
    set_axis_ticks_minmax_zero(gca, "x", xlim(gca))
end
set_axis_ticks_minmax_zero(gca, "y", ylim(gca))

if ~isempty(p.Results.x_label)
    horizontal_label = p.Results.x_label;
end
if ~isempty(p.Results.y_label)
    y_label = p.Results.y_label;
end

xlabel(horizontal_label)
ylabel(y_label)
set(gca, "FontSize", 6, "XColor", "k", "YColor", "k", "GridColor", "k")
xlabel(horizontal_label, "Color", "k")
ylabel(y_label, "Color", "k")
title(ttl, "FontSize", 6, "FontWeight", "normal", "Color", "k")
end

function [slice_axis, horizontal_values] = infer_vertical_slice_axis(field_size_raw, field_size_squeezed, g)

slice_axis = "x";
horizontal_values = g.x(:) - g.x(end) / 2;

if numel(field_size_raw) >= 3
    if field_size_raw(1) == 1 && field_size_raw(2) == numel(g.y) && field_size_raw(3) == numel(g.z)
        slice_axis = "y";
        horizontal_values = g.y(:) - g.y(end) / 2;
        return
    end
    if field_size_raw(2) == 1 && field_size_raw(1) == numel(g.x) && field_size_raw(3) == numel(g.z)
        return
    end
end

if isequal(field_size_squeezed, [numel(g.y) numel(g.z)]) && ~isequal(field_size_squeezed, [numel(g.x) numel(g.z)])
    slice_axis = "y";
    horizontal_values = g.y(:) - g.y(end) / 2;
end
end
