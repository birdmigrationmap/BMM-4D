function radar_xy_km = get_plot_radar_xy(grid, radar)
%GET_PLOT_RADAR_XY Convert radar lon/lat metadata to centered plot coordinates.

if isempty(radar)
    radar_xy_km = [0 0];
    return
end

radar_xy_km = [ ...
    (grid.f_X(radar.lon) - grid.x(end) / 2) ./ 1000 ...
    (grid.f_Y(radar.lat) - grid.y(end) / 2) ./ 1000];
end
