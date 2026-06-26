function water_mask = get_water_mask(grid)
%GET_WATER_MASK Return the horizontal water mask for a reconstruction grid.

if isfield(grid, "mask_water")
    water_mask = grid.mask_water;
    return
end

water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
water_mask = logical(water_mask_data.Fonwater({grid.f_lat(grid.y), grid.f_lon(grid.x)}))';
end
