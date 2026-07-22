function metrics = summarize_static_recovery_metrics(truth_field, inverse_field, g, bins)
%SUMMARIZE_STATIC_RECOVERY_METRICS Compute static synthetic recovery metrics.

truth_value = double(truth_field(bins.valid_mask));
inverse_value = double(inverse_field(bins.valid_mask));
residual = inverse_value - truth_value;
truth_rms = sqrt(mean(truth_value .^ 2));
normalized_residual = residual ./ truth_rms;

metrics.rmse = sqrt(mean(residual .^ 2));
metrics.nrmse = sqrt(mean(normalized_residual .^ 2));
metrics.bias = mean(residual);
metrics.correlation = corr(truth_value, inverse_value);
if ~isfinite(metrics.correlation)
    metrics.correlation = nan;
end
descriptor = summarize_static_truth_descriptors({truth_field}, g, table("case", VariableNames="station"));
metrics.truth_rms = descriptor.truth_rms;
metrics.mean_altitude_m = descriptor.mean_altitude_m;
metrics.low_altitude_fraction = descriptor.low_altitude_fraction;
metrics.horizontal_roughness = descriptor.horizontal_roughness;

truth_sum = accumarray(bins.bin_id, truth_value, [bins.n_bin 1], @sum, 0);
inverse_sum = accumarray(bins.bin_id, inverse_value, [bins.n_bin 1], @sum, 0);
truth_sq_sum = accumarray(bins.bin_id, truth_value .^ 2, [bins.n_bin 1], @sum, 0);
inverse_sq_sum = accumarray(bins.bin_id, inverse_value .^ 2, [bins.n_bin 1], @sum, 0);
truth_inverse_sum = accumarray(bins.bin_id, truth_value .* inverse_value, [bins.n_bin 1], @sum, 0);
n_value = accumarray(bins.bin_id, ones(size(truth_value)), [bins.n_bin 1], @sum, 0);

metrics.range_height_rmse = reshape(sqrt(accumarray(bins.bin_id, residual .^ 2, [bins.n_bin 1], @sum, 0) ./ n_value), bins.bin_size);
metrics.range_height_nrmse = reshape(sqrt(accumarray(bins.bin_id, normalized_residual .^ 2, [bins.n_bin 1], @sum, 0) ./ n_value), bins.bin_size);
metrics.range_height_correlation = reshape( ...
    (truth_inverse_sum - truth_sum .* inverse_sum ./ n_value) ./ ...
    sqrt((truth_sq_sum - truth_sum .^ 2 ./ n_value) .* (inverse_sq_sum - inverse_sum .^ 2 ./ n_value)), ...
    bins.bin_size);
metrics.range_height_correlation(~isfinite(metrics.range_height_correlation)) = nan;
end
