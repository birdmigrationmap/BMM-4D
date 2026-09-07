function add_velocity_arrow(ax, uv_profile_time)

uv_mean = squeeze(mean(uv_profile_time, 1));
if norm(uv_mean) == 0
    return
end

arrow_dx_km = uv_mean(1) * 3.6;
arrow_dy_km = uv_mean(2) * 3.6;

hold(ax, "on")
quiver(ax, 0, 0, arrow_dx_km, arrow_dy_km, 0, ...
    "Color", "w", "LineWidth", 2.2, "MaxHeadSize", 0.55)
quiver(ax, 0, 0, arrow_dx_km, arrow_dy_km, 0, ...
    "Color", "k", "LineWidth", 1.2, "MaxHeadSize", 0.45)
end
