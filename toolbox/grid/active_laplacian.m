function L = active_laplacian(active_mask_3d, water_mask, coeff_horver, store_class)
%ACTIVE_LAPLACIAN 3D smoothness operator restricted to active cells.
%   Same operator as toolbox/grid/LaplacianMatrix.m with water=true, but
%   built directly on the active cells (the full-grid operator does not fit
%   in memory for a continental grid). Neighbouring active cells are linked
%   with weight (1 + coeff_horver)/6 horizontally and (1 - 2 coeff_horver)/6
%   vertically; horizontal land-water pairs are cut; the diagonal is the sum
%   of the remaining link weights (inactive neighbours act as the grid
%   boundary). With every cell active, L equals LaplacianMatrix(g,
%   water=true, coeff_horver=coeff_horver).

if nargin < 4
    store_class = "double";
end
sz = size(active_mask_3d, 1:3);
active_id = find(active_mask_3d(:));
n_active = numel(active_id);
position = zeros(prod(sz), 1, "uint32");
position(active_id) = 1:n_active;
[ix, iy, iz] = ind2sub(sz, active_id);
water = water_mask(:);
xy_id = sub2ind(sz(1:2), ix, iy);

strides = [1, sz(1), sz(1) * sz(2)];
limits = {ix < sz(1), iy < sz(2), iz < sz(3)};
weights = [(1 + coeff_horver) / 6, (1 + coeff_horver) / 6, (1 - 2 * coeff_horver) / 6];
rows = cell(3, 1);
cols = cell(3, 1);
vals = cell(3, 1);
diagonal = zeros(n_active, 1);
for d = 1:3
    from = find(limits{d});
    neighbor = active_id(from) + strides(d);
    to = double(position(neighbor));
    keep = to > 0;
    if d < 3
        % Cut horizontal links between land and water columns.
        keep = keep & water(xy_id(from)) == water(xy_id(from) + strides(d));
    end
    from = from(keep);
    to = to(keep);
    rows{d} = [from; to];
    cols{d} = [to; from];
    vals{d} = -weights(d) .* ones(2 * numel(from), 1);
    diagonal = diagonal + weights(d) .* (accumarray(from, 1, [n_active 1]) + accumarray(to, 1, [n_active 1]));
end
L = sparse([vertcat(rows{:}); (1:n_active)'], [vertcat(cols{:}); (1:n_active)'], [vertcat(vals{:}); diagonal], n_active, n_active);
if store_class == "single"
    L = single(L);
end
end
