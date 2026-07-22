function water_mask = get_water_mask(g)
%GET_WATER_MASK Return the horizontal water mask for a reconstruction g.

if isfield(g, "mask_water")
    water_mask = g.mask_water;
    return
end

water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
water_mask = logical(water_mask_data.Fonwater({g.f_lat(g.y), g.f_lon(g.x)}))';
end
