function column_mask = build_domain_mask(g, domain)
%BUILD_DOMAIN_MASK Geographic solve domain [x y]: US land plus a coastal strip.
%   Land of the contiguous US (state polygons) east of DOMAIN.west_lon, and
%   water (sea and Great Lakes) within DOMAIN.offshore_km of that land.
%   Canadian and Mexican land are outside.
%   Using a fixed geographic domain instead of the radar ranges avoids maps
%   shaped by the radar positions.

[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
states = states(~ismember({states.Name}, ["Alaska" "Hawaii"]));
in_us = false(g.sz(1:2));
for i_s = 1:numel(states)
    box = states(i_s).BoundingBox;
    near = lon_grid >= box(1, 1) & lon_grid <= box(2, 1) & lat_grid >= box(1, 2) & lat_grid <= box(2, 2);
    in_us(near) = in_us(near) | inpolygon(lon_grid(near), lat_grid(near), states(i_s).Lon, states(i_s).Lat);
end
us_land = in_us & ~g.mask_water & lon_grid >= domain.west_lon;

% Distance (km) from each cell to the nearest US land cell
distance_km = bwdist(us_land) * mean([g.dx g.dy]) / 1000;
us_water = g.mask_water & distance_km <= domain.offshore_km & lon_grid >= domain.west_lon;

column_mask = us_land | us_water;
end
