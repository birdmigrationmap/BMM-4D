function result = gd_estimate_sweep(log10_lo, log10_hi, range_m, az_step_deg, beam, table_s, varargin)
%GD_ESTIMATE_SWEEP Homogeneous-density estimate around every gate of a sweep.
%   For each gate, all gates within RADIUS_M along range and within the
%   same arc length along azimuth (a range x arc box, the window shape of
%   the read_rd smoothing, approximating an x-y radius) are assumed to
%   share one bird density. Each gate is an interval [LOG10_LO, LOG10_HI]
%   of its summed bird contribution s (GD_GATE_BOUNDS); -Inf/Inf for both
%   ends means no information. The gate's log probability is computed for
%   every density of DENSITY_GRID from the Poisson mixture at its own range
%   (V ∝ r^2), then summed over the window with moving sums, and the
%   window maximum is refined by a parabola in log density.
%   RESULT has density, lower and upper (approximate 95% likelihood
%   interval, too narrow where radials overlap) [n_range n_az], and
%   n_known (gates with information in the window).
%   Options: radius_m (2000), density_grid (logspace(-1, 3.7, 48)).

p = inputParser;
addParameter(p, "radius_m", 2000);
addParameter(p, "density_grid", logspace(-1, 3.7, 48));
parse(p, varargin{:});
opt = p.Results;
range_m = range_m(:);
[n_range, n_az] = size(log10_lo);
grid_d = opt.density_grid(:)';
n_d = numel(grid_d);
n_s = numel(table_s.log10_s);
n = (0:table_s.n_max)';

% Mixture CDF of s for every range row and density: [n_s, n_d] per row
known = ~(isinf(log10_lo) & isinf(log10_hi));
index_of = @(x) min(max(round((x - table_s.log10_s(1)) / (table_s.log10_s(2) - table_s.log10_s(1))) + 1, 1), n_s);
lo_index = index_of(log10_lo);
hi_index = index_of(log10_hi);
log_p = zeros(n_range, n_az, n_d);
volume = beam.volume_sample_km3(range_m);
for i_r = 1:n_range
    if ~any(known(i_r, :))
        continue
    end
    lambda = grid_d * volume(i_r);
    poisson = exp(n .* log(lambda) - lambda - gammaln(n + 1));
    cdf = table_s.cdf * poisson;                                  % [n_s, n_d]
    cdf_lo = cdf(lo_index(i_r, :), :);
    cdf_hi = cdf(hi_index(i_r, :), :);
    cdf_lo(isinf(log10_lo(i_r, :)) & log10_lo(i_r, :) < 0, :) = 0;
    cdf_hi(isinf(log10_hi(i_r, :)) & log10_hi(i_r, :) > 0, :) = 1;
    row = log(max(cdf_hi - cdf_lo, 1e-12));
    row(~known(i_r, :), :) = 0;
    row(:, 1 - sum(poisson, 1) > 1e-3) = -1e6;                    % beyond the table
    log_p(i_r, :, :) = reshape(row, 1, n_az, n_d);
end

% Window sums (GD_WINDOW_SUM)
window_ll = gd_window_sum(log_p, range_m, az_step_deg, opt.radius_m);
n_known = gd_window_sum(double(known), range_m, az_step_deg, opt.radius_m);

% Maximum, parabolic refinement in log density, 1.92 interval
[best, i_best] = max(window_ll, [], 3);
log_grid = log(grid_d);
step = log_grid(2) - log_grid(1);
i_mid = min(max(i_best, 2), n_d - 1);
pick = @(k) window_ll(sub2ind(size(window_ll), repmat((1:n_range)', 1, n_az), repmat(1:n_az, n_range, 1), k));
y_minus = pick(i_mid - 1);
y_mid = pick(i_mid);
y_plus = pick(i_mid + 1);
curvature = y_plus - 2 * y_mid + y_minus;
shift = 0.5 * (y_minus - y_plus) ./ curvature;
shift(~(curvature < 0) | abs(shift) > 1) = 0;
density = exp(log_grid(i_mid) + shift * step);
inside = window_ll >= best - 1.92;
first = nan(n_range, n_az);
last = nan(n_range, n_az);
[~, first(:)] = max(reshape(inside, [], n_d), [], 2);
[~, last_flip] = max(fliplr(reshape(inside, [], n_d)), [], 2);
last(:) = n_d + 1 - last_flip;
lower = grid_d(first);
upper = grid_d(last);
no_data = n_known == 0;
density(no_data) = nan;
lower(no_data) = nan;
upper(no_data) = nan;
result = struct("density", density, "lower", lower, "upper", upper, "n_known", n_known);
end
