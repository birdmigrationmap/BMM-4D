function is_transition = find_land_water_transitions(g, source_id, dest_id)
%FIND_LAND_WATER_TRANSITIONS Identify horizontal connections crossing shorelines.

water_mask = bmm4d.grid.get_water_mask(g);
is_transition = water_mask(source_id) ~= water_mask(dest_id);
end
