function [range_height_metrics, case_metrics] = aggregate_static_recovery_metrics(truth_fields, inverse_fields, scans, bins)
%AGGREGATE_STATIC_RECOVERY_METRICS Aggregate static recovery by range-height.

n_bin = bins.n_bin;
sqerr_sum = zeros(n_bin, 1);
normalized_sqerr_sum = zeros(n_bin, 1);
bias_sum = zeros(n_bin, 1);
truth_sum = zeros(n_bin, 1);
inverse_sum = zeros(n_bin, 1);
truth_sq_sum = zeros(n_bin, 1);
inverse_sq_sum = zeros(n_bin, 1);
truth_inverse_sum = zeros(n_bin, 1);
n_value = zeros(n_bin, 1);

case_metrics = table('Size', [0 8], ...
    'VariableTypes', ["double" "double" "string" "double" "double" "double" "double" "double"], ...
    'VariableNames', ["i_case" "i_sim" "station" "rmse" "nrmse" "bias" "correlation" "truth_rms"]);

[n_sim, n_case] = size(truth_fields);
for i_case = 1:n_case
    for i_sim = 1:n_sim
        truth_field = truth_fields{i_sim, i_case};
        inverse_field = inverse_fields{i_sim, i_case};
        truth_value = double(truth_field(bins.valid_mask));
        inverse_value = double(inverse_field(bins.valid_mask));
        residual = inverse_value - truth_value;
        truth_rms = sqrt(mean(truth_value .^ 2));
        normalized_residual = residual ./ truth_rms;

        sqerr_sum = sqerr_sum + accumarray(bins.bin_id, residual .^ 2, [n_bin 1], @sum, 0);
        normalized_sqerr_sum = normalized_sqerr_sum + accumarray(bins.bin_id, normalized_residual .^ 2, [n_bin 1], @sum, 0);
        bias_sum = bias_sum + accumarray(bins.bin_id, residual, [n_bin 1], @sum, 0);
        truth_sum = truth_sum + accumarray(bins.bin_id, truth_value, [n_bin 1], @sum, 0);
        inverse_sum = inverse_sum + accumarray(bins.bin_id, inverse_value, [n_bin 1], @sum, 0);
        truth_sq_sum = truth_sq_sum + accumarray(bins.bin_id, truth_value .^ 2, [n_bin 1], @sum, 0);
        inverse_sq_sum = inverse_sq_sum + accumarray(bins.bin_id, inverse_value .^ 2, [n_bin 1], @sum, 0);
        truth_inverse_sum = truth_inverse_sum + accumarray(bins.bin_id, truth_value .* inverse_value, [n_bin 1], @sum, 0);
        n_value = n_value + accumarray(bins.bin_id, ones(size(residual)), [n_bin 1], @sum, 0);

        pattern_correlation = corr(truth_value, inverse_value);
        if ~isfinite(pattern_correlation)
            pattern_correlation = nan;
        end

        case_metrics = [case_metrics; {i_case, i_sim, scans.station(i_case), ...
            sqrt(mean(residual .^ 2)), sqrt(mean(normalized_residual .^ 2)), mean(residual), pattern_correlation, truth_rms}]; %#ok<AGROW>
    end
end

rmse = sqrt(sqerr_sum ./ n_value);
normalized_rmse = sqrt(normalized_sqerr_sum ./ n_value);
bias = bias_sum ./ n_value;
correlation = (truth_inverse_sum - truth_sum .* inverse_sum ./ n_value) ./ ...
    sqrt((truth_sq_sum - truth_sum .^ 2 ./ n_value) .* (inverse_sq_sum - inverse_sum .^ 2 ./ n_value));
correlation(~isfinite(correlation)) = nan;

range_height_metrics = struct( ...
    "rmse", reshape(rmse, bins.bin_size), ...
    "normalized_rmse", reshape(normalized_rmse, bins.bin_size), ...
    "bias", reshape(bias, bins.bin_size), ...
    "correlation", reshape(correlation, bins.bin_size), ...
    "n_value", reshape(n_value, bins.bin_size));
end
