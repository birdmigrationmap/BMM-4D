function prior_active = ar_profile_prior(g, radar_scans, radar_profiles, active_mask_3d, length_scale_m, vertical, options)
%AR_PROFILE_PRIOR Radar-profile prior, with a choice of vertical rule over terrain (dev).
%   PRIOR_ACTIVE = AR_PROFILE_PRIOR(G, RADAR_SCANS, RADAR_PROFILES,
%   ACTIVE_MASK_3D, LENGTH_SCALE_M, VERTICAL) interpolates the radar
%   profiles of the production prior (build_density_prior_fast,
%   info.radar_profiles [nz n_radar n_time], layers above sea level)
%   between radars with Gaussian weights of LENGTH_SCALE_M
%   (build_gaussian_weights) and returns the active cells [n_active n_time].
%   - VERTICAL "asl": layer by layer, as the production prior (equal to it
%     at the production length scale);
%   - VERTICAL "agl": each radar's profile is read from the lowest active
%     cell of the radar's column upward (height above that radar's
%     ground), and laid in every column from its own lowest active cell
%     upward. Over flat ground this is the "asl" prior; over a ridge the
%     dense low layer rises with the ground instead of being cut by it.
%   - VERTICAL "sigma": terrain-following and count-conserving, as the
%     "sigma" links of ar_steps. Each radar's profile, from its column's
%     lowest to highest active cell, becomes birds per unit of relative
%     height (density x air thickness of the radar's column) on N_SIGMA
%     bins; these are interpolated between radars and spread over each
%     column's own active layers (divided by its air thickness). A column
%     over a ridge then holds the same birds as the radars' columns,
%     squeezed into less air.
%   - VERTICAL "smooth": as the "smooth" links of ar_steps. Each radar's
%     profile is read by height above the ground smoothed over
%     OPTIONS.smooth_length_m (G_L, default 100 km) at the radar, and laid
%     at the same height above G_L in every column: a cell at z takes the
%     radar's profile at z + G_L(radar) - G_L(column). Infinite length: the
%     "asl" prior.
%   - VERTICAL "hybrid": as the "hybrid" links of ar_steps, from the
%     radar's column to each column (ground rising by dG = G(column) -
%     G(radar), grid ground): a cell h above its ground takes the radar's
%     profile at h / (1 - dG / H) above the radar's ground where dG < H and
%     h < H - dG, and at the same height above sea level elsewhere (H =
%     OPTIONS.fade_height_m, default 1000 m).
%   The "smooth" and "hybrid" profiles are interpolated linearly in height
%   and held constant beyond the grid's first and last layer.

arguments
    g
    radar_scans
    radar_profiles
    active_mask_3d
    length_scale_m
    vertical (1, 1) string
    options.smooth_length_m (1, 1) double = 100e3
    options.fade_height_m (1, 1) double = 1000
end

nz = g.sz(3);
nxy = prod(g.sz(1:2));
n_radar = size(radar_profiles, 2);
n_time = size(radar_profiles, 3);
weights = bmm4d.prior.build_gaussian_weights(g, radar_scans, length_scale_m);
active_id = find(active_mask_3d(:));
[ax, ay, az] = ind2sub(g.sz, active_id);
active_column = sub2ind(g.sz(1:2), ax, ay);
if vertical == "sigma"
    prior_active = sigma_prior(g, radar_scans, radar_profiles, active_mask_3d, weights);
    return
elseif vertical == "smooth" || vertical == "hybrid"
    prior_active = shifted_prior(g, radar_scans, radar_profiles, active_mask_3d, weights, vertical, options);
    return
end
if vertical == "asl"
    layer = az;                                       % sea-level layer of each active cell
    radar_lowest = ones(n_radar, 1);
else
    [has_active, lowest] = max(active_mask_3d, [], 3);
    lowest(~has_active) = nz + 1;
    layer = az - lowest(active_column) + 1;           % 1 = lowest active cell
    radar_column = sub2ind(g.sz(1:2), ...
        min(max(round(interp1(g.x, 1:g.sz(1), g.f_X(cellfun(@(s) s.lon, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(1)), ...
        min(max(round(interp1(g.y, 1:g.sz(2), g.f_Y(cellfun(@(s) s.lat, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(2)));
    radar_lowest = min(lowest(radar_column), nz);
end

prior_active = zeros(numel(active_id), n_time, "single");
for i_t = 1:n_time
    % Profiles from each radar's lowest layer up (zero above the top of the grid)
    profiles = zeros(nz, n_radar);
    for i_r = 1:n_radar
        above = radar_lowest(i_r):nz;
        profiles(1:numel(above), i_r) = radar_profiles(above, i_r, i_t);
    end
    field = full(weights * profiles');                % nxy x nz
    prior_active(:, i_t) = single(field(active_column + (layer - 1) * nxy));
end
end

function prior_active = sigma_prior(g, radar_scans, radar_profiles, active_mask_3d, weights)
% Terrain-following, count-conserving prior (VERTICAL "sigma")
n_sigma = 60;
nz = g.sz(3);
n_radar = size(radar_profiles, 2);
n_time = size(radar_profiles, 3);
active_id = find(active_mask_3d(:));
[ax, ay, az] = ind2sub(g.sz, active_id);
active_column = sub2ind(g.sz(1:2), ax, ay);
[~, lowest] = max(active_mask_3d, [], 3);
n_layers = sum(active_mask_3d, 3);
air = zeros(g.sz);
air(active_id) = bmm4d.grid.air_fraction(g, active_id);
thickness = sum(air, 3);
radar_column = sub2ind(g.sz(1:2), ...
    min(max(round(interp1(g.x, 1:g.sz(1), g.f_X(cellfun(@(s) s.lon, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(1)), ...
    min(max(round(interp1(g.y, 1:g.sz(2), g.f_Y(cellfun(@(s) s.lat, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(2)));
% Radar columns without active cells (radar outside the domain): full grid
radar_lowest = lowest(radar_column);
radar_layers = n_layers(radar_column);
radar_thickness = thickness(radar_column);
outside = radar_layers == 0;
radar_lowest(outside) = 1;
radar_layers(outside) = nz;
radar_thickness(outside) = nz;
% Bin of each active cell: centre of its share of the column
sigma_bin = min(floor((az - lowest(active_column) + 0.5) ./ n_layers(active_column) * n_sigma) + 1, n_sigma);
sigma_centre = ((1:n_sigma) - 0.5) / n_sigma;
prior_active = zeros(numel(active_id), n_time, "single");
for i_t = 1:n_time
    per_sigma = zeros(n_sigma, n_radar);
    for i_r = 1:n_radar
        layer = radar_lowest(i_r) + min(floor(sigma_centre * radar_layers(i_r)), radar_layers(i_r) - 1);
        per_sigma(:, i_r) = radar_profiles(layer, i_r, i_t) * radar_thickness(i_r);
    end
    field = full(weights * per_sigma');               % nxy x n_sigma, birds per unit relative height
    prior_active(:, i_t) = single(field(active_column + (sigma_bin - 1) * prod(g.sz(1:2))) ./ thickness(active_column));
end
end


function prior_active = shifted_prior(g, radar_scans, radar_profiles, active_mask_3d, weights, vertical, options)
% Each cell reads every radar's profile at the height that the "smooth" or
% "hybrid" rule maps to the radar's column (VERTICAL "smooth" / "hybrid")
z = squeeze(g.z);
n_radar = size(radar_profiles, 2);
n_time = size(radar_profiles, 3);
active_id = find(active_mask_3d(:));
[ax, ay, az] = ind2sub(g.sz, active_id);
active_column = sub2ind(g.sz(1:2), ax, ay);
ground = max(g.dem, 0);
radar_column = sub2ind(g.sz(1:2), ...
    min(max(round(interp1(g.x, 1:g.sz(1), g.f_X(cellfun(@(s) s.lon, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(1)), ...
    min(max(round(interp1(g.y, 1:g.sz(2), g.f_Y(cellfun(@(s) s.lat, radar_scans(:, 1))), "linear", "extrap")), 1), g.sz(2)));
if vertical == "smooth"
    if isinf(options.smooth_length_m)
        reference = zeros(g.sz(1:2));
    else
        reference = imgaussfilt(ground, options.smooth_length_m / g.dx, Padding="replicate");
    end
else
    reference = ground;
    height = max(z(az) - ground(active_column), 0);     % height above the cell's ground
end
% Height read in each radar's column, as a fractional layer index
layer_in_radar = zeros(numel(active_id), n_radar);
for i_r = 1:n_radar
    rise = reference(active_column) - reference(radar_column(i_r));
    if vertical == "smooth"
        z_radar = z(az) - rise;
    else
        z_radar = z(az);                                % level: same height above sea level
        H = options.fade_height_m;
        follow = rise < H & height < H - rise;
        z_radar(follow) = reference(radar_column(i_r)) + height(follow) ./ (1 - rise(follow) / H);
    end
    layer_in_radar(:, i_r) = min(max((z_radar - z(1)) / g.dz + 1, 1), g.sz(3));
end
lower = floor(layer_in_radar);
upper = min(lower + 1, g.sz(3));
fraction = layer_in_radar - lower;
cell_weights = full(weights(active_column, :));       % n_active x n_radar
prior_active = zeros(numel(active_id), n_time, "single");
for i_t = 1:n_time
    value = zeros(numel(active_id), 1);
    for i_r = 1:n_radar
        profile = radar_profiles(:, i_r, i_t);
        value = value + cell_weights(:, i_r) .* ((1 - fraction(:, i_r)) .* profile(lower(:, i_r)) + fraction(:, i_r) .* profile(upper(:, i_r)));
    end
    prior_active(:, i_t) = single(value);
end
end
