function plot_obs(obs_vector, elevation_kernels, info, varargin)
%PLOT_OBS Plot processed radar observations as 3D samples or one PPI panel.

p = inputParser;
addRequired(p, "obs_vector");
addRequired(p, "elevation_kernels");
addRequired(p, "info");
addParameter(p, "ax", []);
addParameter(p, "plot_type", "scatter3");
addParameter(p, "i_elev", []);
addParameter(p, "c_lim", []);
addParameter(p, "r_lim", []);
addParameter(p, "ttl", "");
addParameter(p, "x_label", "x (km)");
addParameter(p, "y_label", "y (km)");
addParameter(p, "show_colorbar", true);
addParameter(p, "min_value", []);
parse(p, obs_vector, elevation_kernels, info, varargin{:});

c_lim = p.Results.c_lim;
r_lim = p.Results.r_lim;
ttl = string(p.Results.ttl);
plot_type = string(p.Results.plot_type);

if isgraphics(p.Results.ax, "axes")
    axes(p.Results.ax)
end

if plot_type == "scatter3"
    hold on; box on; grid(gca, "on")
    [obs, ~, azimuth_deg, elevation_deg, range_m, beam_height_m] = ...
        bmm4d.forward.format_data(obs_vector, elevation_kernels, info);

    w = range_m .* obs;
    w(~isfinite(obs)) = 0;
    w(obs < 1) = 0;

    if any(w > 0)
        id = datasample((1:numel(obs))', min(nnz(w > 0), 100000), ...
            "Replace", false, "Weights", w);
    else
        id = find(isfinite(obs) & obs > 0);
        if isempty(id)
            id = (1:numel(obs))';
        end
    end

    [x, y] = sph2cart(cmp2pol(azimuth_deg(id)), deg2rad(elevation_deg(id)), range_m(id));

    scatter3(x / 1000, y / 1000, beam_height_m(id), 3, obs(id), "filled")
    scatter3(0, 0, 0, 36, "w", "filled", "MarkerEdgeColor", "k", "LineWidth", 0.5)

    view(3)
    if ~isempty(c_lim)
        clim(c_lim)
    end
    if p.Results.show_colorbar
        cb = colorbar;
        cb.FontSize = 6;
    end

    if ~isempty(r_lim)
        xlim([-r_lim r_lim])
        ylim([-r_lim r_lim])
    end
    zlim([0 2000])

    xlabel(p.Results.x_label)
    ylabel(p.Results.y_label)
    zlabel("Altitude (m)")
    set(gca, "FontSize", 6, "XColor", "k", "YColor", "k", "ZColor", "k", "GridColor", "k")
elseif plot_type == "ppi"
    i_elev = p.Results.i_elev;
    obs_matrix = reshape(full(obs_vector), info.obs_grid_size);
    ppi = obs_matrix(:, :, 1, i_elev);
    elevation_kernel = elevation_kernels{i_elev};
    min_value = p.Results.min_value;
    if isempty(min_value)
        if ~isempty(c_lim)
            min_value = max(eps, 0.02 .* c_lim(2));
        else
            positive_values = ppi(isfinite(ppi) & ppi > 0);
            if isempty(positive_values)
                min_value = eps;
            else
                min_value = max(eps, quantile(positive_values, 0.05));
            end
        end
    end

    ppi(~isfinite(ppi) | ppi <= min_value) = nan;

    ppi = fillmissing2([ppi ppi(:, 1)], "linear");
    ppi = ppi(:, 1:end-1);

    [Z, x, y] = mat2cart( ...
        ppi, ...
        elevation_kernel.az, ...
        elevation_kernel.range, ...
        1000, ...
        max(elevation_kernel.x), ...
        "linear");

    imagescnan(x / 1000, y / 1000, Z);
    hold on
    scatter(0, 0, 36, "w", "filled", "MarkerEdgeColor", "k", "LineWidth", 0.5)

    axis xy equal tight
    if ~isempty(c_lim)
        clim(c_lim)
    end
    if p.Results.show_colorbar
        cb = colorbar;
        cb.FontSize = 6;
    end

    r_lim = max(abs([x(:); y(:)])) / 1000;
    xlim([-r_lim r_lim])
    ylim([-r_lim r_lim])
    bmm4d.plot.set_axis_ticks_minmax_zero(gca, "x", xlim(gca))
    bmm4d.plot.set_axis_ticks_minmax_zero(gca, "y", ylim(gca))

    xlabel(p.Results.x_label)
    ylabel(p.Results.y_label)
    set(gca, "FontSize", 6, "XColor", "k", "YColor", "k")
else
    error("plot_obs:UnknownPlotType", ...
        "Unknown plot_type '%s'. Use 'scatter3' or 'ppi'.", plot_type);
end

title(ttl, "FontSize", 6, "FontWeight", "normal", "Color", "k")
end
