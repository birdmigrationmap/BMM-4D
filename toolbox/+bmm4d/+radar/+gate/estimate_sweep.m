function result = estimate_sweep(log10_lo, log10_hi, range_m, az_step_deg, beam, table_s, varargin)
%ESTIMATE_SWEEP Homogeneous-density estimate around every gate of a sweep.
%   For each gate, all gates within RADIUS_M along range and within the
%   same arc length along azimuth (BMM4D.RADAR.GATE.WINDOW_SUM) are assumed
%   to share one bird density. Each gate's log likelihood at every density
%   of DENSITY_GRID (BMM4D.RADAR.GATE.LOGLIK) is summed over the window, and
%   the window maximum is refined by a parabola in log density.
%   RESULT has density, lower and upper (approximate 95% likelihood
%   interval, too narrow where radials overlap) [n_range n_az], log_sd (the
%   standard deviation of log density from the curvature of the window log
%   likelihood at its maximum), sd (standard deviation of the density in
%   birds km^-3: the window likelihood normalized over DENSITY_GRID, uniform
%   prior in density; its spread, also where the likelihood only bounds the
%   density from above; at least density x log_sd, where the likelihood is
%   narrower than the grid step) and n_known (gates with information in the
%   window). Windows without information are NaN.
%   Options: radius_m (1000, a 2 km window), density_grid (logspace(-1, 3.7, 48)).

p = inputParser;
addParameter(p, "radius_m", 1000);
addParameter(p, "density_grid", logspace(-1, 3.7, 48));
parse(p, varargin{:});
opt = p.Results;
range_m = range_m(:);
[n_range, n_az] = size(log10_lo);
grid_d = opt.density_grid(:)';
n_d = numel(grid_d);

[log_p, known] = bmm4d.radar.gate.loglik(log10_lo, log10_hi, range_m, beam, table_s, grid_d);
window_ll = bmm4d.radar.gate.window_sum(log_p, range_m, az_step_deg, opt.radius_m);
n_known = bmm4d.radar.gate.window_sum(double(known), range_m, az_step_deg, opt.radius_m);

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
log_sd = step ./ sqrt(max(-curvature, 1e-6));
inside = window_ll >= best - 1.92;
first = nan(n_range, n_az);
last = nan(n_range, n_az);
[~, first(:)] = max(reshape(inside, [], n_d), [], 2);
[~, last_flip] = max(fliplr(reshape(inside, [], n_d)), [], 2);
last(:) = n_d + 1 - last_flip;
lower = grid_d(first);
upper = grid_d(last);
% Spread of the normalized likelihood in density (grid cells of width ~ density)
weight = exp(window_ll - best) .* reshape(grid_d, 1, 1, []);
weight = weight ./ sum(weight, 3);
mean_density = sum(weight .* reshape(grid_d, 1, 1, []), 3);
sd = sqrt(max(sum(weight .* reshape(grid_d, 1, 1, []) .^ 2, 3) - mean_density .^ 2, 0));
clear weight
% A likelihood narrower than the grid step collapses on one grid density:
% there the curvature (log_sd) gives the spread
sd = max(sd, density .* log_sd);
no_data = n_known == 0;
density(no_data) = nan;
lower(no_data) = nan;
upper(no_data) = nan;
log_sd(no_data) = nan;
sd(no_data) = nan;
result = struct("density", density, "lower", lower, "upper", upper, "log_sd", log_sd, "sd", sd, "n_known", n_known);
end
