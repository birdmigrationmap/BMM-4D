function [interpolation_weights, radar_xy] = build_horizontal_idw_weights(grid, radar_scans, varargin)
%BUILD_HORIZONTAL_IDW_WEIGHTS Build shared radar-to-grid IDW interpolation weights.
% The returned sparse matrix has size [numel(grid.x)*numel(grid.y) n_radar]
% and maps one scalar value per radar to one scalar field on the full
% horizontal grid. Duplicate radar XY positions are handled internally by
% computing site-level weights and splitting each site's weight equally
% across all radars at that same location.

p = inputParser;
addParameter(p, "n_neighbors", 4);
addParameter(p, "power", 2);
parse(p, varargin{:});

n_radar = size(radar_scans, 1);
nxy = prod(grid.sz(1:2));
radar_xy = nan(n_radar, 2);

for i_r = 1:n_radar
    radar_xy(i_r, :) = [grid.f_X(radar_scans{i_r, 1}.lon) grid.f_Y(radar_scans{i_r, 1}.lat)];
end

[radar_xy_unique, ~, site_id] = unique(radar_xy, "rows", "stable");
n_site = size(radar_xy_unique, 1);
site_count = accumarray(site_id, 1, [n_site 1]);
n_neighbors = min(p.Results.n_neighbors, n_site);
[X, Y] = ndgrid(grid.x(:), grid.y(:));
grid_xy = [X(:) Y(:)];
chunk_size = 10000;
row_ids = zeros(nxy * n_neighbors, 1);
col_ids = zeros(nxy * n_neighbors, 1);
values = zeros(nxy * n_neighbors, 1);

for row_start = 1:chunk_size:nxy
    row_end = min(row_start + chunk_size - 1, nxy);
    chunk_ids = row_start:row_end;
    dx = grid_xy(chunk_ids, 1) - radar_xy_unique(:, 1)';
    dy = grid_xy(chunk_ids, 2) - radar_xy_unique(:, 2)';
    [distance_squared, site_ids] = mink(dx.^2 + dy.^2, n_neighbors, 2);
    weights = distance_squared.^(-p.Results.power / 2);

    zero_distance = distance_squared == 0;
    rows_with_zero = any(zero_distance, 2);
    weights(rows_with_zero, :) = zero_distance(rows_with_zero, :);
    weights = weights ./ sum(weights, 2);

    output_ids = (row_start - 1) * n_neighbors + (1:numel(chunk_ids) * n_neighbors);
    row_ids(output_ids) = repelem(chunk_ids(:), n_neighbors);
    col_ids(output_ids) = reshape(site_ids', [], 1);
    values(output_ids) = reshape(weights', [], 1);
end

site_weights = sparse(row_ids, col_ids, values, nxy, n_site);
site_membership = sparse(site_id, 1:n_radar, 1 ./ site_count(site_id), n_site, n_radar);
interpolation_weights = site_weights * site_membership;
end
