function [density, interval, log_likelihood] = gd_estimate_density(log10_lo, log10_hi, volume_sample_km3, table_s, density_grid)
%GD_ESTIMATE_DENSITY Bird density of a window of gates, assumed homogeneous.
%   Each gate is known only as an interval of its summed bird contribution
%   s = eta * V_eff (log10 cm^2): [LOG10_LO, LOG10_HI]. Use -Inf / Inf for
%   open ends (no echo: below the detection limit; filtered gates: between
%   the output and the output plus the power removed). The number of birds
%   per gate is Poisson with mean density * VOLUME_SAMPLE_KM3 (the gates'
%   mean sampled volume) and their sum follows TABLE_S
%   (GD_BIRD_SUM_TABLE), so each gate has probability
%   F(hi) - F(lo) with F = sum_N Poisson(N) cdf_N. DENSITY (birds km^-3)
%   maximizes the summed log probability over DENSITY_GRID; INTERVAL is
%   the approximate 95% likelihood interval (log likelihood within 1.92
%   of the maximum). Gates are treated as independent, so the interval is
%   too narrow where radials overlap.

lambda = density_grid(:)' * volume_sample_km3;
n = (0:table_s.n_max)';
poisson = exp(n .* log(lambda) - lambda - gammaln(n + 1));        % [n_max + 1, n_grid]
too_many = 1 - sum(poisson, 1) > 1e-3;
mixture_cdf = table_s.cdf * poisson;                              % [n_s, n_grid]
lo = max(min(log10_lo(:), table_s.log10_s(end)), table_s.log10_s(1));
hi = max(min(log10_hi(:), table_s.log10_s(end)), table_s.log10_s(1));
cdf_lo = interp1(table_s.log10_s, mixture_cdf, lo);
cdf_hi = interp1(table_s.log10_s, mixture_cdf, hi);
cdf_lo(isinf(log10_lo(:)) & log10_lo(:) < 0, :) = 0;
cdf_hi(isinf(log10_hi(:)) & log10_hi(:) > 0, :) = 1;
log_likelihood = sum(log(max(cdf_hi - cdf_lo, 1e-12)), 1);
log_likelihood(too_many) = -Inf;                                   % beyond the table: N_max too small
[best, i_best] = max(log_likelihood);
density = density_grid(i_best);
inside = find(log_likelihood >= best - 1.92);
interval = density_grid([inside(1) inside(end)]);
end
