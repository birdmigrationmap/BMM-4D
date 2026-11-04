function table_s = bird_sum_table(beam, n_max, n_sample)
%BIRD_SUM_TABLE Distribution of the summed bird contribution in a gate.
%   For N = 0..N_MAX birds placed uniformly in the counting volume of a
%   gate (BMM4D.RADAR.GATE.BEAM), the gate sum s = sum(sigma_i w_i) (cm^2),
%   with w the two-way weight f^4 W^2, depends neither on range nor on the
%   beam widths (birds are uniform in the scaled ellipse of the beam). Its
%   cumulative distribution on the grid TABLE_S.log10_s is drawn by Monte
%   Carlo (N_SAMPLE gates); beyond N_MAX the sum is taken as normal
%   (central limit), with the exact mean and variance of one bird.
%   The number of birds in the counting volume is Poisson with mean
%   lambda = density * V_sample(r), so a gate's distribution depends on
%   lambda only: TABLE_S.cdf(:, k) is the Poisson mixture at
%   lambda = 10^log10_lambda(k), on a grid fine enough (1% steps) to be
%   interpolated (LOGLIK). A gate's reflectivity is s / V_eff(r).
%   The draws use their own random stream (seed 1), so the table is the
%   same every time and the caller's random numbers are untouched.

stream = RandStream("mt19937ar", Seed=1);
% Angular weight f^4 with u = (a / (R width))^2 ~ U(0, 1), times range weight W^2 with x / sigma_r ~ U(-reach, reach)
weight_of = @(u, v) exp(-8 * log(2) * beam.sample_radius ^ 2 * u) .* exp(-(beam.range_reach * (2 * v - 1)) .^ 2 / 2);
mu = log(beam.sigma_mean_cm2) - beam.sigma_log_sd ^ 2 / 2;
log10_s = (-8:0.02:5)';
s_grid = 10 .^ log10_s;
edges = [-Inf; log10_s];
cdf_n = ones(numel(log10_s), n_max + 1);
% Add one bird at a time to each of the N_SAMPLE gates (a running sum keeps memory small for large N_MAX)
sums = zeros(n_sample, 1);
for n = 1:n_max
    sums = sums + exp(mu + beam.sigma_log_sd * randn(stream, n_sample, 1)) .* weight_of(rand(stream, n_sample, 1), rand(stream, n_sample, 1));
    cdf_n(:, n + 1) = cumsum(histcounts(log10(sums), edges))' / n_sample;
end
% Exact moments of one bird's contribution: sigma, f^4 and W^2 are independent
a = 8 * log(2) * beam.sample_radius ^ 2;
reach = beam.range_reach;
mean_one = beam.sigma_mean_cm2 * (1 - exp(-a)) / a * sqrt(2 * pi) / (2 * reach) * erf(reach / sqrt(2));
var_one = beam.sigma_mean_cm2 ^ 2 * exp(beam.sigma_log_sd ^ 2) * (1 - exp(-2 * a)) / (2 * a) ...
    * sqrt(pi) / (2 * reach) * erf(reach) - mean_one ^ 2;

% Poisson mixture over N at each lambda, keeping the counts within 10 sd of the mean.
% For large counts, neighbouring N give nearly the same distribution, so counts are
% grouped in blocks of about sqrt(lambda) / 4 and each block uses the cdf at its mean N.
log10_lambda = -5:0.005:5.5;
cdf = zeros(numel(log10_s), numel(log10_lambda));
for k = 1:numel(log10_lambda)
    lambda = 10 ^ log10_lambda(k);
    spread = 10 * sqrt(lambda) + 10;
    n = (max(0, floor(lambda - spread)):ceil(lambda + spread))';
    probability = exp(n .* log(lambda) - lambda - gammaln(n + 1));
    block = floor((n - n(1)) / max(1, floor(sqrt(lambda) / 4))) + 1;
    weight = accumarray(block, probability);
    center = round(accumarray(block, probability .* n) ./ max(weight, realmin));
    in_table = center <= n_max;
    cdf(:, k) = cdf_n(:, center(in_table) + 1) * weight(in_table);
    if any(~in_table)
        mean_sum = center(~in_table)' * mean_one;
        sd_sum = sqrt(center(~in_table)' * var_one);
        cdf(:, k) = cdf(:, k) + 0.5 * erfc(-(s_grid - mean_sum) ./ (sqrt(2) * sd_sum)) * weight(~in_table);
    end
end
table_s = struct("log10_s", log10_s, "log10_lambda", log10_lambda, "cdf", cdf, "beam", beam);
end
