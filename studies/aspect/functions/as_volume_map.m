function [map, x_km, y_km] = as_volume_map(volume, sweep_ids, max_range_m, cell_m)
%AS_VOLUME_MAP Mean density of low sweeps on a Cartesian grid around the radar (dev).
%   [MAP, X_KM, Y_KM] = AS_VOLUME_MAP(VOLUME, SWEEP_IDS, MAX_RANGE_M, CELL_M)
%   averages the finite thinned gates of the sweeps SWEEP_IDS of VOLUME
%   (ground range within MAX_RANGE_M) on square cells of CELL_M metres
%   centred on the radar. MAP is [numel(X_KM) numel(Y_KM)] ([x y], NaN
%   where no gate); X_KM, Y_KM are the cell centres (km east and north).

edges = -max_range_m:cell_m:max_range_m;
n = numel(edges) - 1;
total = zeros(n, n);
count = zeros(n, n);
for i_elev = sweep_ids
    sweep = volume.dens.sweeps(i_elev);
    [i_range, i_az] = find(isfinite(sweep.data));
    ground_m = slant2ground(sweep.range(i_range), sweep.elev);
    ground_m = ground_m(:);
    az = sweep.az(i_az);
    x = ground_m .* sind(az(:));
    y = ground_m .* cosd(az(:));
    ix = discretize(x, edges);
    iy = discretize(y, edges);
    ok = isfinite(ix) & isfinite(iy) & ground_m <= max_range_m;
    values = sweep.data(sub2ind(size(sweep.data), i_range(ok), i_az(ok)));
    total = total + accumarray([ix(ok) iy(ok)], values, [n n]);
    count = count + accumarray([ix(ok) iy(ok)], 1, [n n]);
end
map = total ./ count;
x_km = (edges(1:end - 1) + cell_m / 2) / 1000;
y_km = x_km;
end
