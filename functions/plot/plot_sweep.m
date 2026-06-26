function [Z, x_plot, y_plot] = plot_sweep(sweep_data, sweep_az, sweep_range, varargin)
%PLOT_SWEEP Plot a raw radar sweep in polar or local Cartesian coordinates.
% Use coordinate_system="polar" for azimuth-range diagnostics and
% coordinate_system="xy" for a local Cartesian map around the radar.

p = inputParser;
addRequired(p, "sweep_data");
addRequired(p, "sweep_az");
addRequired(p, "sweep_range");
addParameter(p, "ax", []);
addParameter(p, "coordinate_system", "xy");
addParameter(p, "max_range", []);
addParameter(p, "cart_resolution", 1000);
addParameter(p, "interpolation", "nearest");
addParameter(p, "c_lim", []);
addParameter(p, "ttl", "");
addParameter(p, "x_label", "");
addParameter(p, "y_label", "");
addParameter(p, "show_colorbar", []);
addParameter(p, "range_unit", "m");
addParameter(p, "alpha_data", []);
parse(p, sweep_data, sweep_az, sweep_range, varargin{:});

coordinate_system = string(p.Results.coordinate_system);
range_unit = string(p.Results.range_unit);
interpolation = string(p.Results.interpolation);
x_label = string(p.Results.x_label);
y_label = string(p.Results.y_label);
ttl = string(p.Results.ttl);
max_range = p.Results.max_range;
c_lim = p.Results.c_lim;
show_colorbar = p.Results.show_colorbar;
if isempty(show_colorbar)
    show_colorbar = isempty(c_lim);
end

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end
hold on; box on

if coordinate_system == "polar"
    id_range = true(size(sweep_range));
    if ~isempty(max_range)
        id_range = sweep_range <= max_range;
    end
    Z = sweep_data(id_range, :);
    x_plot = sweep_az(:)';
    y_plot = sweep_range(id_range);
    if range_unit == "km"
        y_plot = y_plot ./ 1000;
    end
    alpha_data = p.Results.alpha_data;
    if ~isempty(alpha_data)
        alpha_data = alpha_data(id_range, :);
        imagesc(x_plot, y_plot, Z, AlphaData=alpha_data);
    else
        imagescnan(x_plot, y_plot, Z);
    end
    set(gca, "YDir", "normal")
    xlim([0 365])
    ylim([min(y_plot) max(y_plot)])
    if x_label == ""
        xlabel("Azimuth (deg)")
    else
        xlabel(x_label)
    end
    if y_label == ""
        if range_unit == "km"
            ylabel("Range (km)")
        else
            ylabel("Range (m)")
        end
    else
        ylabel(y_label)
    end
elseif coordinate_system == "xy"
    if isempty(max_range) || ~isfinite(max_range)
        max_range = max(sweep_range);
    end
    [Z, x_plot, y_plot] = mat2cart( ...
        sweep_data, sweep_az, sweep_range, ...
        p.Results.cart_resolution, max_range, interpolation);
    imagescnan(x_plot, y_plot, Z);
    axis xy equal tight
    if x_label ~= ""
        xlabel(x_label)
    end
    if y_label ~= ""
        ylabel(y_label)
    end
else
    error("plot_sweep:UnknownCoordinateSystem", ...
        "Unknown coordinate_system '%s'. Use 'polar' or 'xy'.", coordinate_system);
end

if ~isempty(c_lim)
    clim(c_lim)
end
if show_colorbar
    cb = colorbar;
    cb.FontSize = 6;
end
title(ttl, "FontSize", 6, "FontWeight", "normal", "Color", "k")
set(gca, "FontSize", 6, "XColor", "k", "YColor", "k")
end
