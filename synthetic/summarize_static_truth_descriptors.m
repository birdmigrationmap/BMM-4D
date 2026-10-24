function truth_descriptors = summarize_static_truth_descriptors(truth_fields, g, scans, case_ids, sim_ids)
%SUMMARIZE_STATIC_TRUTH_DESCRIPTORS Compute truth-only descriptors.

if nargin < 4
    case_ids = 1:size(truth_fields, 2);
end
if nargin < 5
    sim_ids = 1:size(truth_fields, 1);
end

truth_descriptors = table('Size', [0 8], ...
    'VariableTypes', ["double" "double" "string" "double" "double" "double" "double" "double"], ...
    'VariableNames', ["i_case" "i_sim" "station" "truth_rms" "horizontal_roughness" ...
    "horizontal_coherence_20km" "mean_altitude_m" "low_altitude_fraction"]);

x_lag = round(20000 ./ median(diff(g.x)));
y_lag = round(20000 ./ median(diff(g.y)));

for i_case = case_ids
    for i_sim = sim_ids
        truth_field = double(truth_fields{i_sim, i_case});
        truth_rms = sqrt(mean(truth_field(:) .^ 2));

        positive_truth = max(truth_field, 0);
        truth_profile = squeeze(sum(positive_truth, [1 2], "omitnan"));
        total_mass = sum(truth_profile, "omitnan");
        mean_altitude_m = sum(truth_profile .* g.z(:), "omitnan") ./ max(total_mass, eps);
        low_altitude_fraction = sum(truth_profile(g.z(:) <= 500), "omitnan") ./ max(total_mass, eps);

        dx_roughness = sqrt(mean(diff(truth_field, 1, 1) .^ 2, "all"));
        dy_roughness = sqrt(mean(diff(truth_field, 1, 2) .^ 2, "all"));
        horizontal_roughness = mean([dx_roughness dy_roughness]) ./ max(truth_rms, eps);

        x_left = truth_field(1:end-x_lag, :, :);
        x_right = truth_field(1+x_lag:end, :, :);
        y_left = truth_field(:, 1:end-y_lag, :);
        y_right = truth_field(:, 1+y_lag:end, :);
        horizontal_coherence_20km = mean([ ...
            corr(x_left(:), x_right(:), Rows="complete") ...
            corr(y_left(:), y_right(:), Rows="complete")], "omitnan");

        truth_descriptors = [truth_descriptors; {i_case, i_sim, scans.station(i_case), truth_rms, ...
            horizontal_roughness, horizontal_coherence_20km, mean_altitude_m, low_altitude_fraction}]; %#ok<AGROW>
    end
end
end
