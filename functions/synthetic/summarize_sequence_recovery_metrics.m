function metrics = summarize_sequence_recovery_metrics(truth_sequence, estimate_sequence, g)
%SUMMARIZE_SEQUENCE_RECOVERY_METRICS Compute 4D synthetic recovery metrics.

truth_value = double(truth_sequence(:));
estimate_value = double(estimate_sequence(:));
residual = estimate_value - truth_value;
truth_rms = sqrt(mean(truth_value .^ 2));

metrics.rmse_4d = sqrt(mean(residual .^ 2));
metrics.nrmse_4d = metrics.rmse_4d ./ max(truth_rms, eps);
metrics.bias_4d = mean(residual);
metrics.correlation_4d = corr(truth_value, estimate_value, Rows="complete");
if ~isfinite(metrics.correlation_4d)
    metrics.correlation_4d = nan;
end
metrics.truth_rms = truth_rms;

% Vertically integrated horizontal-map recovery
vertical_factor_km = g.dz ./ 1000;
truth_integrated_map = squeeze(sum(double(truth_sequence), 3)) .* vertical_factor_km;
estimate_integrated_map = squeeze(sum(double(estimate_sequence), 3)) .* vertical_factor_km;
integrated_map_residual = estimate_integrated_map(:) - truth_integrated_map(:);
integrated_map_truth_rms = sqrt(mean(truth_integrated_map(:) .^ 2));
metrics.integrated_map_rmse = sqrt(mean(integrated_map_residual .^ 2));
metrics.integrated_map_nrmse = metrics.integrated_map_rmse ./ max(integrated_map_truth_rms, eps);
metrics.integrated_map_correlation = corr(truth_integrated_map(:), estimate_integrated_map(:), Rows="complete");
if ~isfinite(metrics.integrated_map_correlation)
    metrics.integrated_map_correlation = nan;
end

% Horizontally averaged vertical-profile recovery
truth_profile = squeeze(mean(double(truth_sequence), [1 2]));
estimate_profile = squeeze(mean(double(estimate_sequence), [1 2]));
profile_residual = estimate_profile(:) - truth_profile(:);
profile_truth_rms = sqrt(mean(truth_profile(:) .^ 2));
metrics.profile_rmse = sqrt(mean(profile_residual .^ 2));
metrics.profile_nrmse = metrics.profile_rmse ./ max(profile_truth_rms, eps);
metrics.profile_correlation = corr(truth_profile(:), estimate_profile(:), Rows="complete");
if ~isfinite(metrics.profile_correlation)
    metrics.profile_correlation = nan;
end

% Physical plausibility and horizontal-detail preservation
negative_values = estimate_value(estimate_value < 0);
positive_values = estimate_value(estimate_value > 0);
metrics.negative_fraction = numel(negative_values) ./ numel(estimate_value);
metrics.negative_mass_fraction = sum(abs(negative_values)) ./ max(sum(positive_values), eps);
truth_dx = diff(double(truth_sequence), 1, 1) ./ g.dx;
truth_dy = diff(double(truth_sequence), 1, 2) ./ g.dy;
estimate_dx = diff(double(estimate_sequence), 1, 1) ./ g.dx;
estimate_dy = diff(double(estimate_sequence), 1, 2) ./ g.dy;
truth_gradient_rms = sqrt((sum(truth_dx(:) .^ 2) + sum(truth_dy(:) .^ 2)) ./ ...
    (numel(truth_dx) + numel(truth_dy)));
estimate_gradient_rms = sqrt((sum(estimate_dx(:) .^ 2) + sum(estimate_dy(:) .^ 2)) ./ ...
    (numel(estimate_dx) + numel(estimate_dy)));
metrics.horizontal_roughness_ratio = estimate_gradient_rms ./ max(truth_gradient_rms, eps);

n_time = size(truth_sequence, 4);
metrics.rmse_time = nan(1, n_time);
metrics.nrmse_time = nan(1, n_time);
metrics.bias_time = nan(1, n_time);
metrics.correlation_time = nan(1, n_time);
metrics.truth_integrated_time = nan(1, n_time);
metrics.estimate_integrated_time = nan(1, n_time);
metrics.integrated_error_time = nan(1, n_time);
for i_t = 1:n_time
    truth_frame = truth_sequence(:, :, :, i_t);
    estimate_frame = estimate_sequence(:, :, :, i_t);
    truth_frame_value = double(truth_frame(:));
    estimate_frame_value = double(estimate_frame(:));
    residual_frame = estimate_frame_value - truth_frame_value;
    truth_frame_rms = sqrt(mean(truth_frame_value .^ 2));

    metrics.rmse_time(i_t) = sqrt(mean(residual_frame .^ 2));
    metrics.nrmse_time(i_t) = metrics.rmse_time(i_t) ./ max(truth_frame_rms, eps);
    metrics.bias_time(i_t) = mean(residual_frame);
    metrics.correlation_time(i_t) = corr(truth_frame_value, estimate_frame_value, Rows="complete");
    metrics.truth_integrated_time(i_t) = sum(truth_frame_value) .* vertical_factor_km;
    metrics.estimate_integrated_time(i_t) = sum(estimate_frame_value) .* vertical_factor_km;
    metrics.integrated_error_time(i_t) = ...
        (metrics.estimate_integrated_time(i_t) - metrics.truth_integrated_time(i_t)) ./ max(abs(metrics.truth_integrated_time(i_t)), eps);
end
metrics.correlation_time(~isfinite(metrics.correlation_time)) = nan;
metrics.median_nrmse_time = median(metrics.nrmse_time, "omitmissing");
metrics.max_nrmse_time = max(metrics.nrmse_time, [], "omitmissing");
metrics.min_correlation_time = min(metrics.correlation_time, [], "omitmissing");
metrics.mean_abs_relative_integrated_error = mean(abs(metrics.integrated_error_time), "omitmissing");
metrics.max_abs_relative_integrated_error = max(abs(metrics.integrated_error_time), [], "omitmissing");
end
