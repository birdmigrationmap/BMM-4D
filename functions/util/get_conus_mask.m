function conus_mask = get_conus_mask(g)
%GET_CONUS_MASK Return a CONUS state-polygon mask on g.x-by-g.y.

[Lon, Lat] = ndgrid(g.lon(:), g.lat(:));
conus_mask = false(g.sz(1:2));
states = shaperead(which("usastatehi.shp"), "UseGeoCoords", true);
for i_state = 1:numel(states)
    if max(states(i_state).Lat, [], "omitnan") < 24 || min(states(i_state).Lat, [], "omitnan") > 50 || ...
            max(states(i_state).Lon, [], "omitnan") < -125 || min(states(i_state).Lon, [], "omitnan") > -66
        continue
    end
    conus_mask = conus_mask | inpolygon(Lon, Lat, states(i_state).Lon, states(i_state).Lat);
end
end
