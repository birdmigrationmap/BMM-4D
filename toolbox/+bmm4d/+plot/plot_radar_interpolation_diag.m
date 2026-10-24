function plot_radar_interpolation_diag(field_xyz, radar_profiles, radar, varargin)

p = inputParser;
addRequired(p, "field_xyz");
addRequired(p, "radar_profiles");
addRequired(p, "radar");
addParameter(p, "g", []);
addParameter(p, "z_ids", []);
addParameter(p, "field_time", []);
addParameter(p, "time_values", []);
addParameter(p, "time_idx", []);
addParameter(p, "figure_name", "Radar interpolation diagnostics");
addParameter(p, "figure_title", "");
addParameter(p, "xy_title_prefix", "XY slice");
addParameter(p, "xz_title_prefix", "XZ slice");
addParameter(p, "tz_title", "Time-altitude slice");
addParameter(p, "tx_title", "Time-x slice");
addParameter(p, "selected_radar_ids", []);
parse(p, field_xyz, radar_profiles, radar, varargin{:});

field_xyz = squeeze(field_xyz);

g = p.Results.g;
if isempty(g)
    error("plot_radar_interpolation_diag:MissingGrid", "Pass g to plot_radar_interpolation_diag.")
end

if isempty(radar)
    radar_xy = [mean(g.x) mean(g.y)];
else
    radar_xy = [g.f_X(radar.lon) g.f_Y(radar.lat)];
end

if isempty(p.Results.z_ids)
    z_ids = unique(max(1, min(g.sz(3), round(1 + (g.sz(3) - 1) .* [0.25 0.75]))));
else
    z_ids = p.Results.z_ids(:)';
end
z_ids = z_ids(1:min(2, numel(z_ids)));
if isscalar(z_ids)
    z_ids(2) = z_ids(1);
end

n_radar = size(radar_xy, 1);
if isempty(p.Results.selected_radar_ids)
    [~, radar_min_y] = min(radar_xy(:, 2));
    [~, radar_max_y] = max(radar_xy(:, 2));
    selected_radar_ids = unique([radar_min_y radar_max_y], "stable");
    if isscalar(selected_radar_ids) && n_radar > 1
        candidate_ids = setdiff(1:n_radar, selected_radar_ids, "stable");
        selected_radar_ids(2) = candidate_ids(1);
    end
else
    selected_radar_ids = p.Results.selected_radar_ids(:)';
end
selected_radar_ids = selected_radar_ids(1:min(2, numel(selected_radar_ids)));
if isscalar(selected_radar_ids)
    selected_radar_ids(2) = selected_radar_ids(1);
end

radar_x_km = (radar_xy(:, 1) - g.x(end) / 2) ./ 1000;
radar_y_km = (radar_xy(:, 2) - g.y(end) / 2) ./ 1000;
y_cut_km = radar_y_km(selected_radar_ids);
x_cut_km = radar_x_km(selected_radar_ids);
z_cut_m = g.z(z_ids);
has_time = ~isempty(p.Results.field_time);
if has_time
    field_time = squeeze(p.Results.field_time);
    if isempty(p.Results.time_values)
        time_values = 1:size(field_time, 4);
    else
        time_values = p.Results.time_values(:)';
    end
    if isempty(p.Results.time_idx)
        time_idx = max(1, ceil(numel(time_values) / 2));
    else
        time_idx = p.Results.time_idx;
    end
end

figure("Name", p.Results.figure_name, "Color", "w");
t = tiledlayout(3, 2, "TileSpacing", "tight", "Padding", "tight");
title(t, p.Results.figure_title)

for i_panel = 1:2
    z_idx = z_ids(i_panel);
    field_xy = field_xyz(:, :, z_idx);

    nexttile
    bmm4d.plot.plot_xy(field_xy, g=g, ttl=sprintf("XY | z = %.0f m", g.z(z_idx)), show_colorbar=true);
    x_lim_xy = xlim(gca);
    y_lim_xy = ylim(gca);
    for i_cut = 1:numel(y_cut_km)
        plot(xlim(gca), y_cut_km(i_cut) .* [1 1], "k--", "LineWidth", 1.2)
    end
    scatter(radar_x_km, radar_y_km, 90, radar_profiles(z_idx, :), "filled", "MarkerEdgeColor", "k")
    xlim(gca, x_lim_xy)
    ylim(gca, y_lim_xy)
end

for i_panel = 1:2
    radar_id = selected_radar_ids(i_panel);
    [~, y_idx] = min(abs(g.y(:) - radar_xy(radar_id, 2)));
    xz_slice = squeeze(field_xyz(:, y_idx, :));

    nexttile
    bmm4d.plot.plot_xz(xz_slice, g=g, ttl=sprintf("XZ | y = %.1f km", radar_y_km(radar_id)), show_colorbar=true);
    x_lim_xz = xlim(gca);
    y_lim_xz = ylim(gca);
    for i_cut = 1:numel(z_cut_m)
        plot(xlim(gca), z_cut_m(i_cut) .* [1 1], "k--", "LineWidth", 1.2)
    end
    scatter(radar_x_km(radar_id) .* ones(g.sz(3), 1), g.z(:), 70, radar_profiles(:, radar_id), ...
        "filled", "MarkerEdgeColor", "k")
    xlim(gca, x_lim_xz)
    ylim(gca, y_lim_xz)
end

if has_time
    radar_id = selected_radar_ids(1);
    z_idx = z_ids(1);
    [~, x_idx] = min(abs(g.x(:) - radar_xy(radar_id, 1)));
    [~, y_idx] = min(abs(g.y(:) - radar_xy(radar_id, 2)));
    time_lim = [min(time_values) max(time_values)];

    tz_slice = squeeze(field_time(x_idx, y_idx, :, :));
    tx_slice = squeeze(field_time(:, y_idx, z_idx, :));

    nexttile
    hold on; box on; grid(gca, "on")
    imagesc(time_values, g.z(:), tz_slice)
    set(gca, "YDir", "normal")
    colorbar
    x_lim_tz = time_lim;
    y_lim_tz = ylim(gca);
    plot(time_values(time_idx) .* [1 1], ylim(gca), "k--", "LineWidth", 1.2)
    for i_cut = 1:numel(z_cut_m)
        plot(xlim(gca), z_cut_m(i_cut) .* [1 1], "k--", "LineWidth", 1.2)
    end
    xlim(gca, x_lim_tz)
    ylim(gca, y_lim_tz)
    xlabel("Time")
    ylabel("Altitude (m)")
    title("Time-z")

    nexttile
    hold on; box on; grid(gca, "on")
    imagesc(time_values, g.x(:) ./ 1000 - g.x(end) ./ 2000, tx_slice)
    set(gca, "YDir", "normal")
    colorbar
    x_lim_tx = time_lim;
    y_lim_tx = ylim(gca);
    plot(time_values(time_idx) .* [1 1], ylim(gca), "k--", "LineWidth", 1.2)
    for i_cut = 1:numel(x_cut_km)
        plot(xlim(gca), x_cut_km(i_cut) .* [1 1], "k--", "LineWidth", 1.2)
    end
    xlim(gca, x_lim_tx)
    ylim(gca, y_lim_tx)
    xlabel("Time")
    ylabel("x (km)")
    title("Time-x")
end
end
