function score = as_overlap_score(D, y, rows, weights)
%AS_OVERLAP_SCORE Disagreement of overlapping radars after a correction (dev).
%   SCORE = AS_OVERLAP_SCORE(D, Y, ROWS, WEIGHTS) for corrected log ratios
%   Y (radar a minus radar b, same cell and time; D from
%   as_overlap_design) on the pairs ROWS, with fixed WEIGHTS (the same for
%   every correction, so the scores compare):
%   - disagreement: weighted RMS of Y after the pair-frame-layer fixed
%     effects and the beam slopes (what the two radars still disagree on,
%     the aspect included);
%   - c2, c1 (and se): aspect harmonics left in Y (as as_overlap_fit, but
%     with the fixed weights): 0 when the correction removes the aspect.

y = y(rows);
w = weights(rows);
[~, ~, group] = unique(D.group(rows));
sw = accumarray(group, w);
demean = @(values) values - index(accumarray(group, w .* values) ./ sw, group);
y_t = demean(y);
B = D.X(rows, D.beam_columns);
B = B(:, any(B ~= 0, 1));
for j = 1:size(B, 2)
    B(:, j) = demean(B(:, j));
end
H = D.X(rows, [1 3]);
for j = 1:2
    H(:, j) = demean(H(:, j));
end
% Disagreement after the beam slopes only
beta_b = (B' * (w .* B)) \ (B' * (w .* y_t));
e = y_t - B * beta_b;
score.disagreement = sqrt(sum(w .* e .^ 2) / sum(w));
% Aspect harmonics left, standard errors clustered by radar pair
X = [H B];
A = X' * (w .* X);
beta = A \ (X' * (w .* y_t));
r = y_t - X * beta;
[~, ~, cluster] = unique(D.pair(rows));
S = zeros(size(X, 2));
G = zeros(max(cluster), size(X, 2));
for j = 1:size(X, 2)
    G(:, j) = accumarray(cluster, w .* r .* X(:, j));
end
S = G' * G + S;
V = A \ S / A;
score.c2 = beta(1);
score.c1 = beta(2);
score.c2_se = sqrt(V(1, 1));
score.c1_se = sqrt(V(2, 2));
score.n = numel(y);
end

function values = index(values, id)
values = values(id);
end
