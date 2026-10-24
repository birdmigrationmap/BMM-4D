function [density, steps] = gd_read_rd_filters(dbz, range_m, az_step_deg)
%GD_READ_RD_FILTERS Our current cleaning steps of read_rd on one sweep.
%   Copies the steps of read_rd (functions/radar/read_rd.m, prod settings)
%   that act on a bird-only sweep, in order, so their effect can be
%   measured on simulated sweeps: density cap (thr_dens = 3000 birds
%   km^-3, set in clean_volume), sweep-wide outlier removal (median +- 3
%   MAD in dBZ), spatial continuity (5 x 5 moving median +- 3 MAD), no-echo
%   gates as 0, then moving medians over 2 km along range and along
%   azimuth (smoothing = 2000). DBZ is [n_range n_az] with NaN for no echo.
%   STEPS holds the density after each step (struct of [n_range n_az]).

density_of = @(x) z_to_refl(idb(x)) / 11;
no_echo = isnan(dbz);

% Density cap
dbz(density_of(dbz) > 3000) = nan;
steps.cap = density_of(dbz);

% Sweep-wide outliers: median +- 3 MAD (scaled) in dBZ
sweep_median = median(dbz(:), "omitnan");
mad_scaled = -1 / (sqrt(2) * erfcinv(3 / 2)) * median(abs(dbz(:) - sweep_median), "omitnan");
dbz(dbz < sweep_median - 3 * mad_scaled | dbz > sweep_median + 3 * mad_scaled) = nan;
steps.outliers = density_of(dbz);

% Spatial continuity: 5 x 5 moving median +- 3 times the mean local MAD
neighbours = nan([size(dbz) 25]);
k = 1;
for i = -2:2
    for j = -2:2
        neighbours(:, :, k) = circshift(dbz, [i j]);
        k = k + 1;
    end
end
local_median = median(neighbours, 3, "omitnan");
local_mad = -1 / (sqrt(2) * erfcinv(3 / 2)) * mean(median(abs(neighbours - local_median), 3, "omitnan"), "all", "omitnan");
dbz(dbz < local_median - 3 * local_mad | dbz > local_median + 3 * local_mad) = nan;
steps.continuity = density_of(dbz);

% No-echo gates are observations of zero density
dbz(no_echo) = -inf;
density = density_of(dbz);
steps.no_echo = density;

% Moving medians over 2 km along range, then along azimuth
window_range = 2000 / median(diff(range_m));
density = smoothdata(density, 1, "movmedian", window_range, "omitnan");
window_az = max(1, round(2000 ./ (deg2rad(az_step_deg) * range_m(:))));
for w = unique(window_az)'
    if w > 1
        rows = window_az == w;
        density(rows, :) = smoothdata(density(rows, :), 2, "movmedian", w, "omitnan");
    end
end
steps.smoothing = density;
end
