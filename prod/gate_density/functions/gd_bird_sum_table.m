function table_s = gd_bird_sum_table(beam, n_max, n_sample)
%GD_BIRD_SUM_TABLE Distribution of the summed bird contribution in a gate.
%   For N = 0..N_MAX birds placed uniformly in the sampled disc of a gate
%   (GD_BEAM), the gate sum s = sum(sigma_i w_i) (cm^2) does not depend on
%   range. TABLE_S.cdf(:, N + 1) is its cumulative distribution on the grid
%   TABLE_S.log10_s, from N_SAMPLE Monte Carlo draws (s = 0 for N = 0, so
%   that column is 1 everywhere). A gate's reflectivity is s / V_eff(r).

rng(1)
weight_of = @(u) exp(-8 * log(2) * beam.sample_radius ^ 2 * u);      % u = (a / (R width))^2 ~ U(0, 1)
mu = log(beam.sigma_mean_cm2) - beam.sigma_log_sd ^ 2 / 2;
contribution = exp(mu + beam.sigma_log_sd * randn(n_sample, n_max)) .* weight_of(rand(n_sample, n_max));
sums = cumsum(contribution, 2);
log10_s = (-8:0.02:5)';
cdf = ones(numel(log10_s), n_max + 1);
for n = 1:n_max
    cdf(:, n + 1) = mean(log10(sums(:, n))' <= log10_s, 2);
end
table_s = struct("log10_s", log10_s, "cdf", cdf, "n_max", n_max, "beam", beam);
end
