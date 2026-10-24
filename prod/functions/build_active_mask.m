function active_mask_3d = build_active_mask(g, radar_scans, coverage_range_m, flight_layer_agl_m, column_mask)
%BUILD_ACTIVE_MASK Cells solved by the inversion.
%   Columns within COVERAGE_RANGE_M of at least one radar, or the columns of
%   COLUMN_MASK ([x y], e.g. build_domain_mask) when given, and in each column
%   the cells from the ground (cell top above g.dem) up to FLIGHT_LAYER_AGL_M
%   above it, capped by the grid top. Water columns are kept (birds over the
%   sea are part of the signal); g.dem is about 0 there.

if nargin < 5 || isempty(column_mask)
    column_mask = radar_columns(g, radar_scans, coverage_range_m);
end
ground = max(g.dem, 0);
active_mask_3d = column_mask ...
    & g.z + g.dz / 2 > ground ...
    & g.z - g.dz / 2 < ground + flight_layer_agl_m;
end

function column_mask = radar_columns(g, radar_scans, coverage_range_m)
% Columns within coverage_range_m of a radar.
radar_x = g.f_X(cellfun(@(s) s.lon, radar_scans(:, 1)));
radar_y = g.f_Y(cellfun(@(s) s.lat, radar_scans(:, 1)));
[X, Y] = ndgrid(g.x, g.y);
distance_to_radar = inf(g.sz(1:2));
for i_r = 1:numel(radar_x)
    distance_to_radar = min(distance_to_radar, hypot(X - radar_x(i_r), Y - radar_y(i_r)));
end
column_mask = distance_to_radar <= coverage_range_m;
end
