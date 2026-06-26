function plot_volume_adv(fields, varargin)
p = inputParser;
addRequired(p, "fields");
addParameter(p, "ax", []);
addParameter(p, "ttl", "");
addParameter(p, "grid", []);
addParameter(p, "radar", []);
addParameter(p, "c_lim", []);
addParameter(p, "r_lim", []);
addParameter(p, "z_slice", []);
addParameter(p, "frame_pause", 0.2);
addParameter(p, "x_label", "x (km)");
addParameter(p, "y_label", "y (km)");
parse(p, fields, varargin{:});

ttl = p.Results.ttl;
g = p.Results.grid;
c_lim = p.Results.c_lim;
r_lim = p.Results.r_lim;
n_time = size(fields, 4);

if isempty(g)
    sz = size(fields);
    x_plot = (1:sz(1))' - sz(1) / 2;
    y_plot = (1:sz(2))' - sz(2) / 2;
    z_plot = (1:sz(3))';
    time_values = [];
else
    x_plot = g.x(:) - g.x(end) / 2;
    y_plot = g.y(:) - g.y(end) / 2;
    z_plot = g.z(:);
    radar_xy_km = get_plot_radar_xy(g, p.Results.radar);
    if isfield(g, "time") && numel(g.time) == n_time
        time_values = g.time(:);
    else
        time_values = [];
    end
end

if isempty(g)
    radar_xy_km = [];
end

if isempty(p.Results.z_slice)
    z_slice = z_plot(min(4, numel(z_plot)));
else
    z_slice = p.Results.z_slice;
end

[z_plot_slice, ~] = extend_top_slice_for_plot(z_plot, permute(fields(:, :, :, 1), [2 1 3]), z_slice);

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end

for i_t = 1:n_time
    cla
    hold on; box on; grid on

    field_plot = permute(fields(:, :, :, i_t), [2 1 3]);
    [~, field_plot_slice] = extend_top_slice_for_plot(z_plot, field_plot, z_slice);
    h = slice(x_plot / 1000, y_plot / 1000, z_plot_slice, field_plot_slice, 0, 0, z_slice);
    set(h, "EdgeColor", "none", "FaceAlpha", 0.95)

    if ~isempty(radar_xy_km)
        scatter3(radar_xy_km(:, 1), radar_xy_km(:, 2), zeros(size(radar_xy_km, 1), 1), ...
            250, "w", "filled", "MarkerEdgeColor", "k")
    end

    axis tight
    view(3)
    if ~isempty(c_lim)
        clim(c_lim)
    end
    colorbar

    if ~isempty(r_lim)
        xlim([-r_lim r_lim])
        ylim([-r_lim r_lim])
        zlim([0 max(z_plot_slice)])
    else
        zlim([0 max(z_plot_slice)])
    end

    xlabel(p.Results.x_label)
    ylabel(p.Results.y_label)
    zlabel("Altitude (m)")

    if isempty(time_values)
        title(ttl + " | t = " + i_t)
    else
        title(ttl + " | " + string(time_values(i_t)))
    end

    drawnow
    pause(p.Results.frame_pause)
end

function [z_plot_out, field_plot_out] = extend_top_slice_for_plot(z_plot, field_plot_in, z_slice)

field_plot_out = field_plot_in;
z_plot_out = z_plot;

if z_slice > max(z_plot)
    z_plot_out(end + 1) = z_slice;
    field_plot_out(:, :, end + 1) = field_plot_out(:, :, end);
end
end
end
