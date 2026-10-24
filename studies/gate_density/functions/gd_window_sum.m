function summed = gd_window_sum(values, range_m, az_step_deg, radius_m)
%GD_WINDOW_SUM Sum over the gates within a radius of each gate of a sweep.
%   VALUES is [n_range n_az (n_layer)]. The window is RADIUS_M along range
%   and the same arc length along azimuth on each side (a range x arc box,
%   the shape of the read_rd smoothing window, approximating an x-y
%   radius); azimuth wraps around north. Uses moving sums, so it is fast.

range_m = range_m(:);
n_az = size(values, 2);
half_range = round(radius_m / median(diff(range_m)));
summed = movsum(values, [half_range half_range], 1, Endpoints=0);
half_az = min(floor(n_az / 2), round(radius_m ./ (range_m * deg2rad(az_step_deg))));
for h = unique(half_az(half_az > 0))'
    rows = half_az == h;
    padded = cat(2, summed(rows, end - h + 1:end, :), summed(rows, :, :), summed(rows, 1:h, :));
    padded = movsum(padded, [h h], 2, Endpoints=0);
    summed(rows, :, :) = padded(:, h + 1:h + n_az, :);
end
end
