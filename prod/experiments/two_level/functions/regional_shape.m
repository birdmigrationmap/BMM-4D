function coefficients = regional_shape(g, radar_scans, radar_profiles, shape_patterns, z_agl, n_pattern, length_scale_m, ground_smoothing_m)
%REGIONAL_SHAPE Shape coefficients of the radar profiles, interpolated to every column.
%   Each radar profile (RADAR_PROFILES [nz n_radar n_time], layers above sea
%   level, from build_density_prior_fast) is shifted to height above the
%   smoothed ground at the radar and projected on the first N_PATTERN shape
%   patterns; its coefficients c_k = a_k / a_1 (k = 2..K) are interpolated
%   to every column with Gaussian weights over LENGTH_SCALE_M, weighted by
%   the radar's column amount so radars with few birds do not set the
%   shape. Radars whose grid top leaves less than the pattern height range
%   above ground get no weight. Returns [nxy K-1 n_time]; c = 0 everywhere
%   is the night's mean shape.

[nz, n_radar, n_time] = size(radar_profiles);
ground_smooth = imgaussfilt(max(g.dem, 0), ground_smoothing_m / mean([g.dx g.dy]));
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
radar_lon = cellfun(@(s) s.lon, radar_scans(:, 1));
radar_lat = cellfun(@(s) s.lat, radar_scans(:, 1));
radar_ground = interp2(lon_grid', lat_grid', ground_smooth', radar_lon, radar_lat);
z_asl = squeeze(g.z);
basis = shape_patterns(:, 1:n_pattern);

radar_coefficients = zeros(n_radar, n_pattern - 1, n_time);
radar_amount = zeros(n_radar, n_time);
for i_r = 1:n_radar
    if z_asl(end) - radar_ground(i_r) < z_agl(end)
        continue
    end
    profile_agl = max(interp1(z_asl - radar_ground(i_r), reshape(radar_profiles(:, i_r, :), nz, n_time), z_agl, "linear", 0), 0);
    amounts = basis \ profile_agl;                 % [K n_time]
    radar_coefficients(i_r, :, :) = reshape(amounts(2:end, :) ./ max(amounts(1, :), eps), 1, [], n_time);
    radar_amount(i_r, :) = sum(profile_agl, 1) * g.dz / 1000;
end

weights = build_gaussian_weights(g, radar_scans, length_scale_m);
coefficients = zeros(prod(g.sz(1:2)), n_pattern - 1, n_time);
for i_t = 1:n_time
    amount = radar_amount(:, i_t);
    total = weights * amount;
    coefficients(:, :, i_t) = (weights * (amount .* radar_coefficients(:, :, i_t))) ./ max(total, eps);
end
end
