function [eta, bird_count] = gd_simulate_sweep(density, range_m, az_step_deg, beam)
%GD_SIMULATE_SWEEP Reflectivity of a sweep through randomly placed birds.
%   DENSITY (birds km^-3) is a scalar, one value per range gate, or a map
%   [n_range n_az] (birds are then thinned by the density of their
%   nearest radial). For each gate centre in RANGE_M, birds fill the shell of one
%   gate length around the radar within +-sample_radius beam widths of the
%   sweep elevation (Poisson count), at uniform range, azimuth and
%   elevation offset, with lognormal cross-sections (BMM4D.RADAR.GATE.BEAM). Each bird
%   adds to every gate within the range reach of its position (two-way
%   weight f^4 W^2), so neighbouring range gates share birds; shells within
%   that reach are also simulated beyond the first and last gate. Radials every
%   AZ_STEP_DEG overlap when the step is below the beam width, as on the
%   WSR-88D (0.5 deg super-resolution radials, 0.95 deg beam), so one bird
%   adds to several radials. ETA [n_range n_az] is the reflectivity
%   (cm^2 km^-3) = sum(sigma w) / V_eff; BIRD_COUNT is the number of birds
%   inside each gate's counting volume.

range_m = range_m(:);
n_range = numel(range_m);
az = (0:az_step_deg:360 - az_step_deg);
n_az = numel(az);
if isvector(density) || isscalar(density)
    density = density(:) .* ones(n_range, n_az);
end
reach_deg = rad2deg(beam.sample_radius * beam.width_rad);                % in elevation
reach_az_deg = rad2deg(beam.sample_radius * beam.width_az_rad);          % in azimuth (rotation broadening)
n_reach = ceil(reach_az_deg / az_step_deg);
reach_m = beam.range_reach * beam.range_sd_m;
spacing_m = median(diff(range_m));
n_reach_range = ceil(reach_m / spacing_m);
mu = log(beam.sigma_mean_cm2) - beam.sigma_log_sd ^ 2 / 2;
eta = zeros(n_range, n_az);
bird_count = zeros(n_range, n_az);
% Shells of one gate spacing around each gate centre, padded by n_reach_range shells on both sides
for i_s = 1 - n_reach_range:n_range + n_reach_range
    i_d = min(max(i_s, 1), n_range);                           % density row of this shell
    shell_centre_m = range_m(1) + (i_s - 1) * spacing_m;
    shell_km3 = spacing_m * shell_centre_m ^ 2 * 2 * reach_deg * pi / 180 * 2 * pi / 1e9;
    row_max = max(density(i_d, :));
    n_bird = poissrnd(row_max * shell_km3);
    if n_bird == 0
        continue
    end
    bird_az = 360 * rand(n_bird, 1);
    % Thinning: keep a bird with probability density / row maximum
    keep = rand(n_bird, 1) < density(i_d, mod(round(bird_az / az_step_deg), n_az) + 1)' / row_max;
    bird_az = bird_az(keep);
    n_bird = numel(bird_az);
    if n_bird == 0
        continue
    end
    bird_range_m = shell_centre_m + spacing_m * (rand(n_bird, 1) - 0.5);
    bird_el = reach_deg * (2 * rand(n_bird, 1) - 1);
    sigma = exp(mu + beam.sigma_log_sd * randn(n_bird, 1));
    nearest = round(bird_az / az_step_deg);
    for i_r = max(1, i_s - n_reach_range):min(n_range, i_s + n_reach_range)
        offset_m = bird_range_m - range_m(i_r);
        in_range = abs(offset_m) <= reach_m;
        if ~any(in_range)
            continue
        end
        w_range = exp(-offset_m .^ 2 / (2 * beam.range_sd_m ^ 2)) .* in_range;
        for offset = -n_reach:n_reach
            radial = mod(nearest + offset, n_az) + 1;
            d_az = mod(bird_az - az(radial)' + 180, 360) - 180;
            a2 = (deg2rad(d_az) / beam.width_az_rad) .^ 2 + (deg2rad(bird_el) / beam.width_rad) .^ 2;   % in beam widths
            inside = a2 <= beam.sample_radius ^ 2;
            w = exp(-8 * log(2) * a2) .* inside .* w_range;
            eta(i_r, :) = eta(i_r, :) + accumarray(radial, sigma .* w, [n_az 1])';
            bird_count(i_r, :) = bird_count(i_r, :) + accumarray(radial, inside & in_range, [n_az 1])';
        end
    end
end
eta = eta ./ beam.volume_eff_km3(range_m);
end
