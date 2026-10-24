function fig = plot_xyt(field_time, varargin)

p = inputParser;
addRequired(p, "field_time");
addParameter(p, "g", []);
addParameter(p, "c_lim", []);
addParameter(p, "frame_ids", []);
addParameter(p, "time_labels", []);
addParameter(p, "figure_name", "XYT");
addParameter(p, "radar", []);
addParameter(p, "radar_marker_size", 80);
parse(p, field_time, varargin{:});

field_time = reshape(field_time, size(field_time, 1), size(field_time, 2), []);
n_frame = size(field_time, 3);

frame_ids = p.Results.frame_ids;
if isempty(frame_ids)
    frame_ids = 1:n_frame;
end

time_labels = p.Results.time_labels;
if isempty(time_labels)
    time_labels = string(frame_ids);
else
    time_labels = string(time_labels);
end

c_lim = p.Results.c_lim;
if isempty(c_lim)
    c_lim = [0 quantile(reshape(field_time(:, :, frame_ids), [], 1), 0.99)];
    c_lim(2) = max(c_lim(2), eps);
end

radar_xy_km = [];
if ~isempty(p.Results.g)
    radar_xy_km = bmm4d.plot.get_plot_radar_xy(p.Results.g, p.Results.radar);
end

fig = bmm4d.plot.create_export_figure(2800, 900, "Name", p.Results.figure_name);
tiledlayout("flow", TileSpacing="tight", Padding="tight");
for i_frame = 1:numel(frame_ids)
    i_t = frame_ids(i_frame);
    nexttile;
    bmm4d.plot.plot_xy(field_time(:, :, i_t), ...
        g=p.Results.g, ...
        radar=p.Results.radar, ...
        c_lim=c_lim, ...
        ttl=time_labels(i_t), ...
        show_colorbar=i_frame == numel(frame_ids));
    if ~isempty(radar_xy_km)
        scatter(radar_xy_km(:, 1), radar_xy_km(:, 2), ...
            p.Results.radar_marker_size, "r", "filled", MarkerEdgeColor="w")
    end
end
end
