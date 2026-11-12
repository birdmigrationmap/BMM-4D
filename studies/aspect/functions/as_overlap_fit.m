function fit = as_overlap_fit(D, columns, rows)
%AS_OVERLAP_FIT Aspect harmonics from radar-pair log ratios (dev).
%   FIT = AS_OVERLAP_FIT(D, COLUMNS, ROWS) regresses the log ratios D.y of
%   as_overlap_design on the regressors D.X(:, COLUMNS) (names or indices)
%   with a fixed effect per pair and frame (removed by weighted demeaning
%   within each group), over the pairs ROWS (logical; default all).
%   Weights 1 / (estimator variance + mismatch variance); the mismatch
%   variance (beams not seeing exactly the same birds) is the mean excess
%   of the squared residuals over the estimator variance (method of
%   moments, two passes), then Huber weights (k = 1.345 robust sd) for the
%   heavy tails. Standard errors are clustered by radar pair (the residuals
%   of one pair are correlated in space and time).
%   FIT: names, beta, se, sigma_mismatch, n, n_pair, and the residuals
%   (resid, demeaned y_tilde) of ROWS for scoring.

if nargin < 3 || isempty(rows)
    rows = true(size(D.y));
end
if isstring(columns) || iscellstr(columns)
    [~, columns] = ismember(string(columns), D.names);
end
% Regressors without any value in ROWS (e.g. beam slope of an empty layer) are dropped
columns = columns(any(D.X(rows, columns) ~= 0, 1));
y = D.y(rows);
X = D.X(rows, columns);
group = D.group(rows);
[~, ~, group] = unique(group);
cluster = D.pair(rows);
var_est = D.var_est(rows);
var_est(~isfinite(var_est)) = median(var_est, "omitnan");

sigma2 = var(y);
huber = ones(size(y));
for pass = 1:4
    w = huber ./ (var_est + sigma2);
    [y_t, X_t] = demean(y, X, w, group);
    beta = pinv(X_t' * (w .* X_t)) * (X_t' * (w .* y_t));
    e = y_t - X_t * beta;
    if pass <= 2
        sigma2 = max(mean(e .^ 2 - var_est), 0.01 * var(y));
    else
        z = e ./ sqrt(var_est + sigma2);
        k = 1.345 * 1.4826 * median(abs(z - median(z)));
        huber = min(1, k ./ abs(z));
    end
end
w = huber ./ (var_est + sigma2);
[y_t, X_t] = demean(y, X, w, group);
A = X_t' * (w .* X_t);
% pinv: nuisance columns can be collinear in small subsets (one radar)
A_inv = pinv(A);
beta = A_inv * (X_t' * (w .* y_t));
e = y_t - X_t * beta;

% Cluster-robust covariance (by radar pair)
score = (w .* e) .* X_t;
[~, ~, c_id] = unique(cluster);
S = zeros(size(X_t, 2));
for j = 1:size(X_t, 2)
    sj = accumarray(c_id, score(:, j));
    for k2 = 1:size(X_t, 2)
        S(j, k2) = sum(sj .* accumarray(c_id, score(:, k2)));
    end
end
V = A_inv * S * A_inv;
fit = struct("names", D.names(columns), "beta", beta', "se", sqrt(diag(V))', ...
    "sigma_mismatch", sqrt(sigma2), "n", numel(y), "n_pair", max(c_id), "resid", e, "y_tilde", y_t, ...
    "weights", w);
end

function [y_t, X_t] = demean(y, X, w, group)
% Weighted demeaning within each group (pair and frame fixed effect)
sw = accumarray(group, w);
y_t = y - mean_by(y, w, group, sw);
X_t = X;
for j = 1:size(X, 2)
    X_t(:, j) = X(:, j) - mean_by(X(:, j), w, group, sw);
end
end

function m = mean_by(values, w, group, sw)
group_mean = accumarray(group, w .* values) ./ sw;
m = group_mean(group);
end
