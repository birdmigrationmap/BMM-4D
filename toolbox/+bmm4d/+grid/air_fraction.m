function fraction = air_fraction(g, active_id)
%AIR_FRACTION Share of each active cell that is above the ground.
%   FRACTION = BMM4D.GRID.AIR_FRACTION(G, ACTIVE_ID) is [n_active 1], 1 for
%   cells fully above the ground (max(g.dem, 0)) and the part above it for
%   the lowest cell of a column, (cell top - ground) / dz. Density is per
%   km^3 of air, so this share enters the forward operator (birds below the
%   ground do not echo) and every column sum; without it a column jumps by
%   a whole layer where the ground crosses a layer boundary.

[ix, iy, iz] = ind2sub(g.sz, active_id(:));
ground = max(g.dem, 0);
z = squeeze(g.z);
fraction = min(max((z(iz) + g.dz / 2 - ground(sub2ind(g.sz(1:2), ix, iy))) / g.dz, 0), 1);
end
