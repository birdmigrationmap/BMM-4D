function [eta, bird_count] = gd_simulate_sweep(density, range_m, az_step_deg, beam)
%GD_SIMULATE_SWEEP Reflectivity of a sweep through randomly placed birds.
%   DENSITY (birds km^-3) is a scalar, one value per range gate, or a map
%   [n_range n_az] (birds are then thinned by the density of their
%   nearest radial). For each gate centre in RANGE_M, birds fill the shell of one
%   gate length around the radar within +-sample_radius beam widths of the
%   sweep elevation (Poisson count), at uniform azimuth and elevation
%   offset, with lognormal cross-sections (GD_BEAM). Radials every
%   AZ_STEP_DEG overlap when the step is below the beam width, as on the
%   WSR-88D (0.5 deg super-resolution radials, 0.95 deg beam), so one bird
%   adds to several radials. ETA [n_range n_az] is the reflectivity
%   (cm^2 km^-3) = sum(sigma w) / V_eff; BIRD_COUNT is the number of birds
%   inside each gate's sampled disc.

range_m = range_m(:);
n_range = numel(range_m);
az = (0:az_step_deg:360 - az_step_deg);
n_az = numel(az);
if isvector(density) || isscalar(density)
    density = density(:) .* ones(n_range, n_az);
end
reach_deg = rad2deg(beam.sample_radius * beam.width_rad);
n_reach = ceil(reach_deg / az_step_deg);
mu = log(beam.sigma_mean_cm2) - beam.sigma_log_sd ^ 2 / 2;
eta = zeros(n_range, n_az);
bird_count = zeros(n_range, n_az);
for i_r = 1:n_range
    shell_km3 = beam.gate_length_m * range_m(i_r) ^ 2 * 2 * reach_deg * pi / 180 * 2 * pi / 1e9;
    row_max = max(density(i_r, :));
    n_bird = poissrnd(row_max * shell_km3);
    if n_bird == 0
        continue
    end
    bird_az = 360 * rand(n_bird, 1);
    % Thinning: keep a bird with probability density / row maximum
    keep = rand(n_bird, 1) < density(i_r, mod(round(bird_az / az_step_deg), n_az) + 1)' / row_max;
    bird_az = bird_az(keep);
    n_bird = numel(bird_az);
    if n_bird == 0
        continue
    end
    bird_el = reach_deg * (2 * rand(n_bird, 1) - 1);
    sigma = exp(mu + beam.sigma_log_sd * randn(n_bird, 1));
    nearest = round(bird_az / az_step_deg);
    for offset = -n_reach:n_reach
        radial = mod(nearest + offset, n_az) + 1;
        d_az = mod(bird_az - az(radial)' + 180, 360) - 180;
        a2 = deg2rad(d_az) .^ 2 + deg2rad(bird_el) .^ 2;
        inside = a2 <= (beam.sample_radius * beam.width_rad) ^ 2;
        w = exp(-8 * log(2) * a2 / beam.width_rad ^ 2) .* inside;
        eta(i_r, :) = eta(i_r, :) + accumarray(radial, sigma .* w, [n_az 1])';
        bird_count(i_r, :) = bird_count(i_r, :) + accumarray(radial, inside, [n_az 1])';
    end
    eta(i_r, :) = eta(i_r, :) / beam.volume_eff_km3(range_m(i_r));
end
end
