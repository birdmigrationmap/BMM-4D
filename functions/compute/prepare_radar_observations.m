function observations = prepare_radar_observations(F, obs_matrix, info)
%PREPARE_RADAR_OBSERVATIONS Compact the shared radar observation system.

active_row_id = find(any(obs_matrix, 2) & any(F, 2));
F_active = F(active_row_id, :);
values = obs_matrix(active_row_id, :);

dense_bytes = 8 * numel(values);
sparse_bytes = 16 * nnz(values) + 8 * (size(values, 2) + 1);
if dense_bytes <= sparse_bytes
    values = full(values);
end

observations = struct( ...
    "F", F_active, ...
    "values", values, ...
    "row_ids", uint32(active_row_id), ...
    "obs_grid_size", info.obs_grid_size);

if isfield(info, "row_radar_id")
    observations.row_radar_id = info.row_radar_id(active_row_id);
end
if isfield(info, "row_kernel_id")
    observations.row_kernel_id = info.row_kernel_id(active_row_id);
end
end
