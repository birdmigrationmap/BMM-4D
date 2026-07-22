function radar_xy_km = get_plot_radar_xy(g, radar)
%GET_PLOT_RADAR_XY Convert radar lon/lat metadata to centered plot coordinates.

if isempty(radar)
    radar_xy_km = zeros(0, 2);
    return
end

radar_xy_km = [ ...
    (g.f_X(radar.lon) - g.x(end) / 2) ./ 1000 ...
    (g.f_Y(radar.lat) - g.y(end) / 2) ./ 1000];
end
