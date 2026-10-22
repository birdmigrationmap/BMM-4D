function observations = prepare_radar_observations(F, obs_matrix, info, varargin)
%PREPARE_RADAR_OBSERVATIONS Compact the shared radar observation system.
%   Kernel rows sum to 1 when the whole beam lies inside the grid. A row cut
%   by the grid top or edges pushes its full measured value into the cells
%   that remain, so real-data cases keep only rows with at least
%   MIN_ROW_MASS of their beam weight inside the grid (0.75 for S09-S11).
%   The default 0 keeps every row: synthetic truth has no birds outside the
%   grid, so truncated rows are unbiased there.

p = inputParser;
addParameter(p, "min_row_mass", 0);
parse(p, varargin{:});

row_mass = full(sum(F, 2));
active_row_id = find(any(obs_matrix, 2) & any(F, 2) & row_mass >= p.Results.min_row_mass);
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
