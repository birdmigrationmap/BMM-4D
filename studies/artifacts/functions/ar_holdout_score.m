function score = ar_holdout_score(S, radar_name, solve)
%AR_HOLDOUT_SCORE Prediction error of one radar left out of the solve (dev).
%   SCORE = AR_HOLDOUT_SCORE(S, RADAR_NAME, SOLVE) gives the rows of
%   RADAR_NAME weight 0 and calls SOLVE(S_OUT, KEEP), KEEP marking the
%   other radars: SOLVE must build its prior from the profiles of the KEEP
%   radars only (e.g. ar_profile_prior with S.radar_profiles(:, KEEP, :)),
%   otherwise the prior still carries the left-out radar's own profile.
%   Then compares the prediction of that radar's observed rows at 30-120 km
%   with its data:
%   - relative_error: relative L2 error;
%   - median_abs_log_ratio: median |log(predicted / observed)|, both
%     floored at 0.1 birds km^-3.
% The data of the other radars then decide the field where this radar
% looks, which the misfit of a solve that fits all radars cannot test.
% The winds of the advection still come from all radars.

i_radar = find(string(S.radar_selected.name) == radar_name);
rows = S.observations.row_radar_id == i_radar;
test_rows = rows & S.observations.row_range_m >= 30000 & S.observations.row_range_m <= 120000;
keep = true(size(S.radar_scans, 1), 1);
keep(i_radar) = false;
S_out = S;
S_out.row_weights(rows, :) = 0;
result = solve(S_out, keep);
predicted = S.observations.F(test_rows, S.active_id) * double(result.field);
truth = full(S.observations.values(test_rows, :));
mask = truth ~= 0;
score = struct("relative_error", norm(predicted(mask) - truth(mask)) / norm(truth(mask)), ...
    "median_abs_log_ratio", median(abs(log(max(predicted(mask), 0.1) ./ max(truth(mask), 0.1)))));
end
