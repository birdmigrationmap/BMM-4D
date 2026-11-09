function weights = observation_weights(observations, prediction, error_floor)
%OBSERVATION_WEIGHTS Weight of each observation from the estimator's sd.
%   WEIGHTS = OBSERVATION_WEIGHTS(OBSERVATIONS, PREDICTION, ERROR_FLOOR)
%   returns [n_rows n_time] weights 1 / sigma^2 for the data term of
%   solve_advection_lean (0 where there is no observation).
%   OBSERVATIONS.sd is the spread of each window's likelihood (birds
%   km^-3, bmm4d.radar.gate.estimate_sweep). Taken at the window's own
%   reading it would give low readings more weight because they are low,
%   so sigma is the sd the estimator has at the PREDICTION ([n_rows
%   n_time], e.g. F x the density prior): the median sd / density of all
%   observations, by expected birds per gate lambda = density x counting
%   volume (bmm4d.radar.gate.beam) and range, read at the prediction and
%   times it. ERROR_FLOOR (birds km^-3) is added in quadrature for what the
%   estimator does not see (forward model, density varying within a cell);
%   it also caps the largest weights, which set the number of PCG
%   iterations.
%   Weights are scaled to the mean variance, mean(sigma^2) / sigma^2, so
%   that weight 1 stands for an observation of average error, as the
%   equal weights the regularization strengths were tuned with. (Scaled to
%   mean weight 1 instead, weight 1 would stand for the most precise,
%   nearly empty windows and the field would be over-smoothed.)

[row, col, observed] = find(observations.values);
sd = full(observations.sd(sub2ind(size(observations.sd), row, col)));
range_km = double(observations.row_range_m(row)) / 1000;
beam = bmm4d.radar.gate.beam();
volume_sample = beam.volume_sample_km3(range_km * 1000);

% Relative sd of the estimator by lambda and range (empty cells: nearest lambda, then range)
lambda_edges = 10 .^ (-2:0.25:3.5);
range_edges_km = [0 5 10 15 20 30 40 60 80 120 160 250];
n_lambda = numel(lambda_edges) - 1;
n_range = numel(range_edges_km) - 1;
range_bin = min(max(discretize(min(range_km, range_edges_km(end)), range_edges_km), 1), n_range);
lambda_bin_of = @(lambda) min(max(discretize(min(max(lambda, lambda_edges(1)), lambda_edges(end)), lambda_edges), 1), n_lambda);
known = sd > 0 & observed > 1e-6;
relative_sd = accumarray([lambda_bin_of(observed(known) .* volume_sample(known)) range_bin(known)], ...
    sd(known) ./ observed(known), [n_lambda n_range], @median, NaN);
relative_sd = fillmissing(fillmissing(relative_sd, "nearest", 1), "nearest", 2);

% sigma at the prediction, the floor added in quadrature
predicted = max(full(prediction(sub2ind(size(prediction), row, col))), 0.1);
relative = relative_sd(sub2ind([n_lambda n_range], lambda_bin_of(predicted .* volume_sample), range_bin));
sigma2 = (predicted .* relative) .^ 2 + error_floor ^ 2;

weights = zeros(size(observations.values));
weights(sub2ind(size(weights), row, col)) = mean(sigma2) ./ sigma2;
end
