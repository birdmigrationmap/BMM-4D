function [elevation_kernels, scans, n_dropped] = attach_kernels(scans, kernel_dir)
%ATTACH_KERNELS Load the forward kernels used by compact scans.
%   Prod copy of the matching/loading part of toolbox/forward/load_F_elev.m
%   without the velocity profiles (computed once per volume in
%   clean_volume). Sweeps whose pattern has no precomputed kernel are
%   dropped instead of failing the radar; N_DROPPED counts them.

kernel_files = struct2table(dir(fullfile(kernel_dir, "F_*.mat")));
kernel_files.name = string(kernel_files.name);
kernel_name_parts = split(kernel_files.name, "_");

% Match each sweep to a kernel file by elevation and azimuth count.
n_dropped = 0;
kernel_id_by_scan = cell(size(scans));
for i_s = 1:numel(scans)
    sweeps = scans{i_s}.dens.sweeps;
    kernel_id = nan(numel(sweeps), 1);
    for i_elev = 1:numel(sweeps)
        match = find("el" + round(sweeps(i_elev).elev * 10) == kernel_name_parts(:, 2) ...
            & "az" + numel(sweeps(i_elev).az) == kernel_name_parts(:, 4), 1);
        if ~isempty(match)
            kernel_id(i_elev) = match;
        end
    end
    keep = ~isnan(kernel_id);
    n_dropped = n_dropped + nnz(~keep);
    scans{i_s}.dens.sweeps = sweeps(keep);
    kernel_id_by_scan{i_s} = kernel_id(keep);
end

% Load each needed kernel once, ordered by elevation.
kernel_id_unique = unique(vertcat(kernel_id_by_scan{:}));
elevation_kernels = cell(numel(kernel_id_unique), 1);
for i_f = 1:numel(kernel_id_unique)
    kernel_data = load(fullfile(kernel_files.folder{kernel_id_unique(i_f)}, kernel_files.name(kernel_id_unique(i_f))));
    kernel = kernel_data.elevation_kernel;
    kernel.sz_raz = [numel(kernel.range) numel(kernel.az)];
    kernel.F_raz_id = sub2ind(kernel.sz_raz, kernel_data.Felevidc(:, 1), kernel_data.Felevidc(:, 2));
    kernel.x = kernel_data.x;
    kernel.y = kernel_data.y;
    kernel.z = kernel_data.z;
    kernel.sz_coarse_grid = [numel(kernel_data.x) numel(kernel_data.y) numel(kernel_data.z)];
    % Kernel cell ids are stored y-first by compute_F_elev.
    [kernel.F_gr_y, kernel.F_gr_x, kernel.F_gr_z] = ind2sub(kernel.sz_coarse_grid([2 1 3]), kernel_data.Felevidc(:, 3));
    kernel.Fval = kernel_data.Felevfc;
    elevation_kernels{i_f} = kernel;
end
[~, order] = sort(cellfun(@(k) k.elev, elevation_kernels));
elevation_kernels = elevation_kernels(order);
kernel_id_sorted = kernel_id_unique(order);
for i_s = 1:numel(scans)
    [~, scans{i_s}.kernel_idx] = ismember(kernel_id_by_scan{i_s}, kernel_id_sorted);
end
end
