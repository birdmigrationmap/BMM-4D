function gc = coarsen_grid(g, aggregate)
%COARSEN_GRID Coarse solve grid made of blocks of fine grid cells.
%   GC = COARSEN_GRID(G, [fx fy fz]) groups fx x fy x fz cells of G into one
%   cell, matching construct_F(..., aggregate=[fx fy fz]): coarse cell i
%   holds fine cells (i-1)*fx+1:i*fx (the last block may be partial). Cell
%   centres are block centres; the water mask is the block majority and the
%   DEM the block mean. [1 1 1] returns G unchanged.

aggregate = aggregate(:)';
if all(aggregate == 1)
    gc = g;
    return
end
gc = g;
gc.sz = ceil(g.sz ./ aggregate);
gc.dx = g.dx * aggregate(1);
gc.dy = g.dy * aggregate(2);
gc.dz = g.dz * aggregate(3);
gc.x = g.x(1) + ((0:gc.sz(1) - 1) + (aggregate(1) - 1) / (2 * aggregate(1))) * gc.dx;
gc.y = g.y(1) + ((0:gc.sz(2) - 1)' + (aggregate(2) - 1) / (2 * aggregate(2))) * gc.dy;
gc.z = reshape(g.z(1) + ((0:gc.sz(3) - 1) + (aggregate(3) - 1) / (2 * aggregate(3))) * gc.dz, 1, 1, []);
gc.lon = g.f_lon(gc.x);
gc.lat = g.f_lat(gc.y);

% Block statistics of the 2D maps
block_x = ceil((1:g.sz(1))' / aggregate(1));
block_y = ceil((1:g.sz(2)) / aggregate(2));
[bx, by] = ndgrid(block_x, block_y);
block_id = sub2ind(gc.sz(1:2), bx(:), by(:));
gc.mask_water = reshape(accumarray(block_id, double(g.mask_water(:)), [prod(gc.sz(1:2)) 1], @mean) > 0.5, gc.sz(1:2));
gc.dem = reshape(accumarray(block_id, g.dem(:), [prod(gc.sz(1:2)) 1], @mean), gc.sz(1:2));
gc.aggregate = aggregate;
end
