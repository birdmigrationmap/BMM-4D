function support = compute_static_observation_support(F, bins)
%COMPUTE_STATIC_OBSERVATION_SUPPORT Summarize F support in range-height bins.

data_sensitivity_voxel = full(sum(F .^ 2, 1))';
support_weight_voxel = full(sum(F, 1))';
hit_count_voxel = full(sum(F > 0, 1))';
observable_voxel = double(hit_count_voxel > 0);

support = struct( ...
    "data_sensitivity", reshape(accumarray(bins.bin_id, data_sensitivity_voxel(bins.valid_vector), [bins.n_bin 1], @mean, nan), bins.bin_size), ...
    "support_weight", reshape(accumarray(bins.bin_id, support_weight_voxel(bins.valid_vector), [bins.n_bin 1], @mean, nan), bins.bin_size), ...
    "hit_count", reshape(accumarray(bins.bin_id, hit_count_voxel(bins.valid_vector), [bins.n_bin 1], @mean, nan), bins.bin_size), ...
    "observable_fraction", reshape(accumarray(bins.bin_id, observable_voxel(bins.valid_vector), [bins.n_bin 1], @mean, nan), bins.bin_size));
end
