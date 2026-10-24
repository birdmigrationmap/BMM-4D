function [ratio, edges_km, bias_index] = radar_distance_bias(field, active_id, g, radar_selected, edges_km)
%RADAR_DISTANCE_BIAS Column density as a function of distance to the nearest radar.
%   FIELD is [n_active n_time] (birds km^-3) on ACTIVE_ID of grid G. For each
%   time, each column density is divided by its regional mean (Gaussian
%   blur with sigma 150 km over the active columns), so large-scale
%   patterns cancel; the mean of this ratio in each distance bin
%   (EDGES_KM, distance to the nearest radar of RADAR_SELECTED) is RATIO
%   [n_bin n_time] after normalizing each curve to mean 1; it is 1
%   everywhere when the map has no radar imprint. Times with almost no
%   birds (domain mean < 5 birds km^-2) are NaN. BIAS_INDEX [1 n_time] is the RMS of log(RATIO): 0 is unbiased,
%   0.1 is about a 10% systematic difference between distance bins.

if nargin < 5
    edges_km = [0 25 50 75 100 125 150 200 300];
end
[ix, iy, ~] = ind2sub(g.sz, active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for i_r = 1:height(radar_selected)
    d = 6371 * acos(min(1, sind(lat_grid) .* sind(radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(radar_selected.lat(i_r)) .* cosd(lon_grid - radar_selected.lon(i_r))));
    distance_km = min(distance_km, d);
end
active_column = unique(column_id);
is_active = false(g.sz(1:2));
is_active(active_column) = true;
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
bin = discretize(distance_km(active_column), edges_km);
n_bin = numel(edges_km) - 1;
n_time = size(field, 2);
ratio = nan(n_bin, n_time);
for i_t = 1:n_time
    column = reshape(accumarray(column_id, double(field(:, i_t)), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    column = max(column, 0) .* is_active;
    if mean(column(active_column)) < 5
        continue
    end
    % Regional mean by normalized convolution over the active columns
    regional = imgaussfilt(column, sigma_cells) ./ active_weight;
    relative = column(active_column) ./ max(regional(active_column), 1);
    ok = ~isnan(bin);
    bin_mean = accumarray(bin(ok), relative(ok), [n_bin 1], @mean, nan);
    % Normalize the curve to mean 1: only its shape with distance matters
    ratio(:, i_t) = bin_mean / mean(bin_mean, "omitnan");
end
bias_index = sqrt(mean(log(ratio) .^ 2, 1, "omitnan"));
end
