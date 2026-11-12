function [tiles, tile_id] = gd_estimate_tiles(log10_lo, log10_hi, range_m, az_step_deg, beam, table_s, density_grid, tile_m)
%GD_ESTIMATE_TILES Homogeneous-density estimate on non-overlapping tiles.
%   Cuts the sweep into tiles of about TILE_M along range and the same arc
%   along azimuth (the arc set at each tile's centre range), assumes one
%   density per tile and maximises its summed log likelihood
%   (BMM4D.RADAR.GATE.LOGLIK) over DENSITY_GRID. TILES is a table with the tile
%   centre range (km), the estimate, the number of informative gates and
%   the maximum log likelihood (summed over tiles for calibration).
%   TILE_ID [n_range n_az] gives the tile (row of TILES) of every gate.

range_m = range_m(:);
[n_range, n_az, ~] = size(log10_lo);
[log_p, known] = bmm4d.radar.gate.loglik(log10_lo, log10_hi, range_m, beam, table_s, density_grid);
% Tile index of every gate
rows_per_tile = max(1, round(tile_m / median(diff(range_m))));
row_block = ceil((1:n_range)' / rows_per_tile);
tile_id = zeros(n_range, n_az);
next_id = 0;
for b = unique(row_block)'
    rows = row_block == b;
    az_per_tile = max(1, round(tile_m / (mean(range_m(rows)) * deg2rad(az_step_deg))));
    az_block = ceil((1:n_az) / az_per_tile);
    tile_id(rows, :) = next_id + repmat(az_block, nnz(rows), 1);
    next_id = next_id + max(az_block);
end
% Summed log likelihood per tile and density
n_d = numel(density_grid);
member = sparse(tile_id(:), (1:n_range * n_az)', 1, next_id, n_range * n_az);
tile_ll = member * reshape(log_p, [], n_d);                       % [n_tile, n_d]
n_known = member * double(known(:));
[log_likelihood, i_best] = max(tile_ll, [], 2);
density = density_grid(i_best)';
center_km = (member * repmat(range_m, n_az, 1)) ./ full(sum(member, 2)) / 1000;
tiles = table(full(center_km), density(:), full(n_known), full(log_likelihood), ...
    VariableNames=["center_km" "density" "n_known" "log_likelihood"]);
end
