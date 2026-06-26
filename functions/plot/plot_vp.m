function plot_vp(vp, varargin)
p = inputParser;
addRequired(p, "vp");
addParameter(p, "ax", []);
addParameter(p, "grid", []);
addParameter(p, "z", []);
addParameter(p, "std_vp", []);
addParameter(p, "ttl", "");
addParameter(p, "x_lim", []);
addParameter(p, "x_label", "Density");
addParameter(p, "y_label", "Altitude (m)");
parse(p, vp, varargin{:});

if ~isempty(p.Results.grid)
    z_plot = p.Results.grid.z(:);
elseif ~isempty(p.Results.z)
    z_plot = p.Results.z(:);
else
    z_plot = (1:numel(vp))';
end

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end

hold on; box on; grid on

if ~isempty(p.Results.std_vp)
    std_vp = p.Results.std_vp(:);
    x_patch = [vp(:) - std_vp; flipud(vp(:) + std_vp)];
    y_patch = [z_plot; flipud(z_plot)];
    patch(x_patch, y_patch, [0 0.45 0.74], ...
        "FaceAlpha", 0.15, "EdgeColor", "none");
end

plot(vp, z_plot, "k-", "LineWidth", 1.8)

if ~isempty(p.Results.x_lim)
    xlim(p.Results.x_lim)
end
ylim([min(z_plot) max(z_plot)])

xlabel(p.Results.x_label)
ylabel(p.Results.y_label)
title(p.Results.ttl)
end
