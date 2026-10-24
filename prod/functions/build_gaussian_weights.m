function interpolation_weights = build_gaussian_weights(g, radar_scans, length_scale_m)
%BUILD_GAUSSIAN_WEIGHTS Smooth radar-to-grid weights [nxy n_radar].
%   Each grid column averages the radars within 3 length scales with
%   weights exp(-d^2 / (2 L^2)), normalized to sum to 1. Unlike the
%   nearest-neighbour IDW of the paper prior, the weights vary smoothly
%   everywhere (no seams where the neighbour set changes, no peak at each
%   radar). Columns with no radar within 3 L take the nearest radar.

n_radar = size(radar_scans, 1);
nxy = prod(g.sz(1:2));
radar_xy = zeros(n_radar, 2);
for i_r = 1:n_radar
    radar_xy(i_r, :) = [g.f_X(radar_scans{i_r, 1}.lon) g.f_Y(radar_scans{i_r, 1}.lat)];
end
[X, Y] = ndgrid(g.x(:), g.y(:));
grid_xy = [X(:) Y(:)];
cutoff_squared = (3 * length_scale_m) ^ 2;
chunk_size = 20000;
parts = cell(ceil(nxy / chunk_size), 1);
for i_chunk = 1:numel(parts)
    rows = ((i_chunk - 1) * chunk_size + 1):min(i_chunk * chunk_size, nxy);
    distance_squared = (grid_xy(rows, 1) - radar_xy(:, 1)') .^ 2 + (grid_xy(rows, 2) - radar_xy(:, 2)') .^ 2;
    weights = exp(-distance_squared / (2 * length_scale_m ^ 2));
    weights(distance_squared > cutoff_squared) = 0;
    [~, nearest] = min(distance_squared, [], 2);
    isolated = sum(weights, 2) == 0;
    weights(sub2ind(size(weights), find(isolated), nearest(isolated))) = 1;
    weights = weights ./ sum(weights, 2);
    [r, k, w] = find(weights);
    parts{i_chunk} = [rows(r)' k w];
end
triplets = vertcat(parts{:});
interpolation_weights = sparse(triplets(:, 1), triplets(:, 2), triplets(:, 3), nxy, n_radar);
end
