function fit = as_overlap_fit_radar(D, rows, radar_height, n_radar, radar_terms, prior_sd)
%AS_OVERLAP_FIT_RADAR Shared aspect function with fixed azimuthal patterns per radar (dev).
%   FIT = AS_OVERLAP_FIT_RADAR(D, ROWS, RADAR_HEIGHT, N_RADAR, RADAR_TERMS)
%   fits on the pairs ROWS of as_overlap_design the "height" form of
%   as_fit_aspect_model (side-on amplitude b2 = b2_per_km h, front/back b1,
%   beam slopes, pair-frame-layer fixed effects) and, when RADAR_TERMS is
%   true, a fixed azimuthal pattern of every radar in the geographic frame:
%       g_r(az) = p1 cos az + q1 sin az + p2 cos 2 az + q2 sin 2 az,
%   az the compass azimuth from the radar. Partial blockage, calibration by
%   sector and lasting gradients stay in the geographic frame, while the
%   aspect turns with the birds' track: on nights with different tracks
%   the two can be told apart.
%   The radar patterns sum to zero over the radars (each harmonic): they
%   are how a radar differs from the others. A pattern shared by every
%   radar would equal the aspect wherever the track is the same, and in the
%   overlap the density field cancels, so no real shared pattern is
%   expected there.
%   PRIOR_SD (default Inf: none) is a Gaussian prior on each radar harmonic
%   (log units), a ridge penalty 1 / PRIOR_SD^2 on the radar columns. A
%   radar sees each neighbour over a narrow azimuth window, and the
%   pair-frame-layer fixed effect removes what is constant in it, so
%   without a prior the patterns of weakly constrained radars are almost
%   free (amplitudes of 10-70 on the 21-22 Oct night) and, through the
%   sum-to-zero constraint, move the shared terms.
%   Sparse design (4 columns per radar), weights and Huber passes as
%   as_overlap_fit, standard errors clustered by radar pair.
%   FIT: names, beta, se of [b2_per_km b1 beam slopes]; radar [N_RADAR x 4]
%   (p1 q1 p2 q2, NaN for radars without pairs) and radar_se;
%   sigma_mismatch, n, n_pair.

if nargin < 5
    radar_terms = true;
end
if nargin < 6
    prior_sd = inf;
end
height_km = (D.height_asl_m(rows) - (radar_height(D.radar_a(rows)) + radar_height(D.radar_b(rows))) / 2) / 1000;
beam = D.X(rows, D.beam_columns);
beam_used = any(beam ~= 0, 1);
main = [D.X(rows, 1) .* height_km, D.X(rows, 3), beam(:, beam_used)];
names = ["b2_per_km" "b1" D.names(D.beam_columns(beam_used))];
n_main = numel(names);
n = nnz(rows);
X = sparse(main);
C = zeros(0, n_main);
radar_column = [];
if radar_terms
    % Radar a enters with +, radar b with - (log ratio a - b)
    az_a = deg2rad(mod(D.psi_a(rows) + D.track(rows), 360));
    az_b = deg2rad(mod(D.psi_b(rows) + D.track(rows), 360));
    ra = D.radar_a(rows);
    rb = D.radar_b(rows);
    basis = @(az) [cos(az) sin(az) cos(2 * az) sin(2 * az)];
    row_id = repmat((1:n)', 1, 4);
    col_a = 4 * (ra - 1) + (1:4);
    col_b = 4 * (rb - 1) + (1:4);
    R = sparse([row_id(:); row_id(:)], [col_a(:); col_b(:)], [reshape(basis(az_a), [], 1); ...
        reshape(-basis(az_b), [], 1)], n, 4 * n_radar);
    % Radars without pairs drop out
    radar_column = find(full(any(R ~= 0, 1)));
    X = [X R(:, radar_column)];
    % Sum-to-zero constraint of each harmonic over the radars
    harmonic = mod(radar_column - 1, 4) + 1;
    C = [zeros(4, n_main), double((1:4)' == harmonic)];
end
y = D.y(rows);
[~, ~, group] = unique(D.group(rows));
[~, ~, cluster] = unique(D.pair(rows));
var_est = D.var_est(rows);
var_est(~isfinite(var_est)) = median(var_est, "omitnan");

sigma2 = var(y);
huber = ones(n, 1);
for pass = 1:5
    w = huber ./ (var_est + sigma2);
    [y_t, X_t] = demean(y, X, w, group);
    A = full(X_t' * (w .* X_t));
    if radar_terms && isfinite(prior_sd)
        penalty = [zeros(n_main, 1); ones(numel(radar_column), 1) / prior_sd ^ 2];
        A = A + diag(penalty);
    end
    [beta, M] = solve_constrained(A, full(X_t' * (w .* y_t)), C);
    e = y_t - X_t * beta;
    if pass <= 2
        sigma2 = max(mean(e .^ 2 - var_est), 0.01 * var(y));
    elseif pass < 5
        z = e ./ sqrt(var_est + sigma2);
        k = 1.345 * 1.4826 * median(abs(z - median(z)));
        huber = min(1, k ./ abs(z));
    end
end

% Cluster-robust covariance (by radar pair) of the constrained estimate
score = X_t .* (w .* e);
S_cluster = sparse(cluster, 1:n, 1) * score;
V = M * full(S_cluster' * S_cluster) * M;
se = sqrt(max(diag(V), 0))';
fit = struct("names", names, "beta", beta(1:n_main)', "se", se(1:n_main), ...
    "sigma_mismatch", sqrt(sigma2), "n", n, "n_pair", max(cluster), "radar", nan(n_radar, 4), ...
    "radar_se", nan(n_radar, 4));
if radar_terms
    radar = nan(4 * n_radar, 1);
    radar_se = nan(4 * n_radar, 1);
    radar(radar_column) = beta(n_main + 1:end);
    radar_se(radar_column) = se(n_main + 1:end);
    fit.radar = reshape(radar, 4, n_radar)';
    fit.radar_se = reshape(radar_se, 4, n_radar)';
end
end

function [y_t, X_t] = demean(y, X, w, group)
% Weighted demeaning within each pair-frame-layer group (sparse)
n_group = max(group);
G = sparse(group, (1:numel(y))', w, n_group, numel(y));
sw = full(sum(G, 2));
scale = spdiags(1 ./ sw, 0, n_group, n_group);
y_mean = scale * (G * y);
X_mean = scale * (G * X);
y_t = y - y_mean(group);
X_t = X - X_mean(group, :);
end

function [beta, M] = solve_constrained(A, b, C)
% Least squares with the linear constraints C beta = 0 (KKT system); M is
% the constrained inverse of A (top-left block of the KKT inverse), for
% the sandwich covariance
p = size(A, 1);
m = size(C, 1);
K = [A C'; C zeros(m)];
K_inv = pinv(K);
solution = K_inv * [b; zeros(m, 1)];
beta = solution(1:p);
M = K_inv(1:p, 1:p);
end
