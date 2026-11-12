function score = nr_score(S, field, cut_code, drop, truth)
%NR_SCORE Near-radar scores of a reconstruction (dev).
%   SCORE = NR_SCORE(S, FIELD, CUT_CODE, DROP) scores FIELD [n_active
%   n_time] against the case S (ar_load, all observations) when the entries
%   DROP (nr_drop_rows) were left out of the solve:
%   - distance: column against the distance to the nearest radar, each
%     radar and frame relative to its own 60-100 km ring, median over
%     radar-frames (bins 0-10, 10-20, 20-30, 30-40, 40-60, 60-100 km);
%   - radars: per radar, median over frames of its 0-20 km column over its
%     60-100 km column;
%   - residuals: observed and predicted sums of the observed entries (both
%     above 0.5 birds km^-3) by cut class (split / batch / other non-split),
%     range ring and whether the entry was dropped; observed / predicted
%     of the dropped rows is the deficit the rule removed, that of the kept
%     batch rows beyond the cut what it left;
%   - misfit_trusted: relative misfit of the kept entries at 30-120 km.
%   SCORE = NR_SCORE(..., TRUTH) also compares the column with a known
%   truth field (synthetic tests): truth_distance is the median over
%   radar-frames of reconstructed / true column per distance bin, and
%   truth_error the median absolute log ratio per bin.

edges = [0 10 20 30 40 60 100];
reference_bin = 6;
n_time = size(field, 2);

% Distance and nearest radar of every grid column
[lon, lat] = ndgrid(double(S.g.lon(:)), double(S.g.lat(:)));
dist_km = inf(numel(lon), 1);
nearest = zeros(numel(lon), 1);
for k = 1:height(S.radar_selected)
    d = 6371 * hypot(deg2rad(lon(:) - S.radar_selected.lon(k)) .* cosd(S.radar_selected.lat(k)), ...
        deg2rad(lat(:) - S.radar_selected.lat(k)));
    closer = d < dist_km;
    dist_km(closer) = d(closer);
    nearest(closer) = k;
end
bin = discretize(dist_km, edges);
n_radar = height(S.radar_selected);
n_bin = numel(edges) - 1;

% Column per radar, distance bin and frame
column = reshape(ar_columns(S, field), [], n_time);
use = ~isnan(bin) & all(isfinite(column), 2);
bin_mean = @(col) accumarray([nearest(use) bin(use)], col(use), [n_radar n_bin], @mean, NaN);
relative = nan(n_radar, n_time, n_bin);
for t = 1:n_time
    m = bin_mean(column(:, t));
    relative(:, t, :) = m ./ m(:, reference_bin);
end
flat = reshape(relative, [], n_bin);
complete = all(isfinite(flat) & flat > 0, 2);
score.distance = table(edges(1:end-1)', edges(2:end)', median(flat(complete, :), 1)', ...
    VariableNames=["lo_km" "hi_km" "column_ratio"]);
near_column = (relative(:, :, 1) + relative(:, :, 2)) / 2;
score.radars = table(S.radar_selected.name, median(near_column, 2, "omitnan"), ...
    VariableNames=["radar" "column_0_20_ratio"]);

% Observed against predicted by cut class, ring and drop
observed = full(S.observations.values);
predicted = S.observations.F(:, S.active_id) * double(field);
range_km = repmat(double(S.observations.row_range_m) / 1000, 1, n_time);
ring = discretize(range_km, edges);
class = ones(size(observed));                  % 1 split
class(cut_code == 4) = 2;                      % 2 batch
class(cut_code == 2 | cut_code == 3 | cut_code == 5) = 3;   % 3 other non-split
class(cut_code == 0) = 4;                      % 4 unknown
ok = observed > 0.5 & predicted > 0.5 & ~isnan(ring);
key = [class(ok) ring(ok) 1 + drop(ok)];
sz = [4 n_bin 2];
n = accumarray(key, 1, sz);
observed_sum = accumarray(key, observed(ok), sz);
predicted_sum = accumarray(key, predicted(ok), sz);
[i_class, i_ring, i_drop] = ndgrid(1:4, 1:n_bin, 1:2);
class_names = ["split" "batch" "other" "unknown"];
has = n(:) > 0;
score.residuals = table(class_names(i_class(has))', edges(i_ring(has))', edges(i_ring(has) + 1)', i_drop(has) == 2, ...
    n(has), observed_sum(has), predicted_sum(has), observed_sum(has) ./ predicted_sum(has), ...
    VariableNames=["cuts" "lo_km" "hi_km" "dropped" "n" "observed" "predicted" "ratio"]);
trusted = ok & ~drop & range_km >= 30 & range_km <= 120;
score.misfit_trusted = norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted));

% Against a known truth
if nargin > 4
    truth_column = reshape(ar_columns(S, truth), [], n_time);
    ratio = nan(n_radar, n_time, n_bin);
    for t = 1:n_time
        ratio(:, t, :) = bin_mean(column(:, t)) ./ bin_mean(truth_column(:, t));
    end
    flat = reshape(ratio, [], n_bin);
    good = all(isfinite(flat) & flat > 0, 2);
    score.distance.truth_ratio = median(flat(good, :), 1)';
    score.distance.truth_error = median(abs(log(flat(good, :))), 1)';
end
end
