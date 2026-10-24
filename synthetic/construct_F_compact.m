function [F, obs_matrix, info] = construct_F_compact(radar_scans, elevation_kernels, g)
%CONSTRUCT_F_COMPACT Build a sparse forward operator with only active rows.
% The row space contains only retained range-azimuth gates from the supplied
% scan geometry, avoiding the global range-azimuth-radar-elevation layout.

%% Extract radar sites ----------------------------------------------------
radar_sites = cell(height(radar_scans), 1);
for i_r = 1:height(radar_scans)
    site_kernel_idx = cellfun(@(x) [x.kernel_idx], radar_scans(i_r, :), UniformOutput=false);
    radar_sites{i_r}.kernel_idx = unique(vertcat(site_kernel_idx{:}));
    radar_sites{i_r}.lon = radar_scans{i_r, 1}.lon;
    radar_sites{i_r}.lat = radar_scans{i_r, 1}.lat;
    radar_sites{i_r}.height = radar_scans{i_r, 1}.height;
    radar_sites{i_r}.station = radar_scans{i_r, 1}.station;

    [~, lon_id] = min(abs(g.lon - radar_sites{i_r}.lon));
    [~, lat_id] = min(abs(g.lat - radar_sites{i_r}.lat));

    radar_sites{i_r}.lon_id = lon_id - ceil(elevation_kernels{1}.sz_coarse_grid(1) / 2);
    radar_sites{i_r}.lat_id = lat_id - ceil(elevation_kernels{1}.sz_coarse_grid(2) / 2);
    % Kernel layers are relative to the antenna, so shift them by the radar
    % height in whole layers. A ground-referenced radar (height 0) has no shift.
    kernel_dz = median(diff(elevation_kernels{1}.z(:)));
    radar_sites{i_r}.z_id = round(radar_sites{i_r}.height / kernel_dz);
end

%% Precompute retained kernel entries for the shared synthetic geometry ---
kernel_compact = cell(numel(elevation_kernels), 1);
rows_per_kernel = zeros(numel(elevation_kernels), 1);
for i_f = 1:numel(elevation_kernels)
    sweep_id = find(radar_scans{1, 1}.kernel_idx == i_f, 1);
    if isempty(sweep_id)
        kernel_compact{i_f} = struct("entry_id", [], "local_row", [], "n_row", 0);
        continue
    end

    sweep_data = radar_scans{1, 1}.dens.sweeps(sweep_id).data;
    present_obs_id = find(~isnan(sweep_data));
    if isempty(present_obs_id)
        kernel_compact{i_f} = struct("entry_id", [], "local_row", [], "n_row", 0);
        continue
    end

    [i_range, i_az] = ind2sub(size(sweep_data), present_obs_id);
    obs_raz = sub2ind(elevation_kernels{i_f}.sz_raz, i_range, i_az);
    obs_raz = unique(obs_raz, "stable");
    [is_present, local_row] = ismember(elevation_kernels{i_f}.F_raz_id, obs_raz);

    entry_id = find(is_present);
    local_row = local_row(entry_id);
    active_local_rows = unique(local_row, "stable");
    [~, local_row] = ismember(local_row, active_local_rows);

    rows_per_kernel(i_f) = numel(active_local_rows);
    kernel_compact{i_f} = struct( ...
        "entry_id", entry_id, ...
        "local_row", local_row, ...
        "n_row", rows_per_kernel(i_f));
end

kernel_row_offset = [0; cumsum(rows_per_kernel(1:end - 1))];
rows_per_radar = sum(rows_per_kernel);
n_row = rows_per_radar * numel(radar_sites);

%% Build compact sparse F -------------------------------------------------
F_rows = cell(numel(radar_sites), numel(elevation_kernels));
F_cols = cell(numel(radar_sites), numel(elevation_kernels));
F_vals = cell(numel(radar_sites), numel(elevation_kernels));
coverage_fraction = ones(numel(radar_sites), numel(elevation_kernels));

for i_r = 1:numel(radar_sites)
    radar_row_offset = (i_r - 1) * rows_per_radar;
    for i_f = 1:numel(elevation_kernels)
        entry_id = kernel_compact{i_f}.entry_id;
        if isempty(entry_id)
            F_rows{i_r, i_f} = zeros(0, 1);
            F_cols{i_r, i_f} = zeros(0, 1);
            F_vals{i_r, i_f} = zeros(0, 1);
            continue
        end

        grid_x_id = elevation_kernels{i_f}.F_gr_x(entry_id) + radar_sites{i_r}.lon_id;
        grid_y_id = elevation_kernels{i_f}.F_gr_y(entry_id) + radar_sites{i_r}.lat_id;
        grid_z_id = elevation_kernels{i_f}.F_gr_z(entry_id) + radar_sites{i_r}.z_id;
        id_F_keep = grid_x_id >= 1 & grid_x_id <= g.sz(1) & ...
            grid_y_id >= 1 & grid_y_id <= g.sz(2) & ...
            grid_z_id >= 1 & grid_z_id <= g.sz(3);
        coverage_fraction(i_r, i_f) = mean(id_F_keep);

        F_rows{i_r, i_f} = radar_row_offset + kernel_row_offset(i_f) + kernel_compact{i_f}.local_row(id_F_keep);
        F_cols{i_r, i_f} = sub2ind(g.sz, grid_x_id(id_F_keep), grid_y_id(id_F_keep), grid_z_id(id_F_keep));
        F_vals{i_r, i_f} = double(elevation_kernels{i_f}.Fval(entry_id(id_F_keep)));
    end
end

F = sparse(vertcat(F_rows{:}), vertcat(F_cols{:}), vertcat(F_vals{:}), n_row, prod(g.sz));
obs_matrix = sparse(n_row, width(radar_scans));

row_radar_id = repelem(uint16((1:numel(radar_sites))'), rows_per_radar);
row_within_radar = repmat(uint32((1:rows_per_radar)'), numel(radar_sites), 1);
row_kernel_id = repmat(repelem(uint16((1:numel(elevation_kernels))'), rows_per_kernel), numel(radar_sites), 1);

info = struct( ...
    "obs_grid_size", [rows_per_radar 1 numel(radar_sites) 1], ...
    "obs_block_offset", (0:numel(radar_sites) - 1)' .* rows_per_radar, ...
    "row_radar_id", row_radar_id, ...
    "row_within_radar", row_within_radar, ...
    "row_kernel_id", row_kernel_id, ...
    "rows_per_kernel", rows_per_kernel, ...
    "rows_per_radar", rows_per_radar, ...
    "coverage_fraction", coverage_fraction, ...
    "has_partial_coverage", coverage_fraction < 1, ...
    "compact", true);
end
