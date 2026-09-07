function [elevation_kernels, radar_scans] = load_F_elev(filepath, radar_scans)
%LOAD_F_ELEV Load precomputed forward kernels for the requested scans.
%   [ELEVATION_KERNELS, RADAR_SCANS] = LOAD_F_ELEV(FILEPATH, RADAR_SCANS)
%   loads the elevation-specific kernels stored in FILEPATH and attaches
%   the matching kernel index to each processed radar scan.

%% Load the available precomputed kernel files
kernel_files = struct2table(dir(filepath));
kernel_files.name = string(kernel_files.name);
kernel_files = kernel_files(startsWith(kernel_files.name, "F_"), :);
if isempty(kernel_files)
    error("No precomputed F_elev files found in %s.", filepath)
end
kernel_name_parts = split(kernel_files.name, "_");

%% Match the required scan patterns to saved kernels
kernel_idx_by_scan = cell(size(radar_scans));

dfnew=struct(); u=0;
for i_s=1:numel(radar_scans)
    if ~isempty(radar_scans{i_s})
        kernel_idx_by_scan{i_s} = nan(numel(radar_scans{i_s}.dens.sweeps),1);
        for i_elev=1:numel(radar_scans{i_s}.dens.sweeps)
            % Find the match in the files
            matching_kernel_idx = find( ...
                "el" + num2str(round(radar_scans{i_s}.dens.sweeps(i_elev).elev * 10)) == kernel_name_parts(:, 2) ...
                & "az" + num2str(numel(radar_scans{i_s}.dens.sweeps(i_elev).az)) == kernel_name_parts(:, 4));
            
            if isempty(matching_kernel_idx)
                tmp_name = "el"+num2str(round(radar_scans{i_s}.dens.sweeps(i_elev).elev*10)) + "_r" + num2str(numel(radar_scans{i_s}.dens.sweeps(i_elev).range))+ "_az" + num2str(numel(radar_scans{i_s}.dens.sweeps(i_elev).az));
                warning("scan missing: "+tmp_name)
                u=u+1;
                dfnew(u).elev = radar_scans{i_s}.dens.sweeps(i_elev).elev;
                [~, dfnew(u).range, dfnew(u).az] = sweep2mat(radar_scans{i_s}.dens.sweeps(i_elev));
                dfnew(u).az = round(dfnew(u).az / 0.25)*.25;
                dfnew(u).name = tmp_name;
                i_s
            else
                kernel_idx_by_scan{i_s}(i_elev) = matching_kernel_idx;
            end
        end
    end
end

if u>0
    save('data/reference/scan_pattern_new','dfnew')
    warning("New scan patterns need to be computed. Saved in 'data/reference/scan_pattern_new'.")
    return
end

%% Load the required kernels

kernel_idx_unique = unique(vertcat(kernel_idx_by_scan{:}));

elevation_kernels = cell(size(kernel_idx_unique));
for i_f = 1:numel(elevation_kernels)
        kernel_file = kernel_files(kernel_idx_unique(i_f), :);
        kernel_data = load(fullfile(kernel_file.folder{:}, kernel_file.name{:}));
    
        elevation_kernels{i_f} = kernel_data.elevation_kernel;
        
        elevation_kernels{i_f}.sz_raz = [numel(kernel_data.elevation_kernel.range) numel(kernel_data.elevation_kernel.az)];
        elevation_kernels{i_f}.F_raz_id = sub2ind(elevation_kernels{i_f}.sz_raz, kernel_data.Felevidc(:,1), kernel_data.Felevidc(:,2));
        
        elevation_kernels{i_f}.x = kernel_data.x;
        elevation_kernels{i_f}.y = kernel_data.y;
        elevation_kernels{i_f}.z = kernel_data.z;
        elevation_kernels{i_f}.sz_coarse_grid = [numel(kernel_data.x) numel(kernel_data.y) numel(kernel_data.z)];
        
        [elevation_kernels{i_f}.F_gr_y, elevation_kernels{i_f}.F_gr_x, elevation_kernels{i_f}.F_gr_z] = ind2sub(elevation_kernels{i_f}.sz_coarse_grid, kernel_data.Felevidc(:,3));
    
        elevation_kernels{i_f}.Fval = kernel_data.Felevfc;
    
    %     F = sparse(double(elevation_kernels{i_f}.F_raz_id), double(kernel_data.Felevidc(:,3)), double(kernel_data.Felevfc), prod(elevation_kernels{i_f}.sz_raz), prod(elevation_kernels{i_f}.sz_coarse_grid));
    %     tmp2 = reshape(full(sum(F,1)), elevation_kernels{i_f}.sz_coarse_grid);
end

%% Order

[~,I] = sort(cellfun(@(x) x.elev, elevation_kernels));
elevation_kernels = elevation_kernels(I);
kernel_idx_unique_sorted = kernel_idx_unique(I);

for i_s = 1:numel(kernel_idx_by_scan)
    if ~isempty(kernel_idx_by_scan{i_s})
        [~, radar_scans{i_s}.kernel_idx] = ismember(kernel_idx_by_scan{i_s}, kernel_idx_unique_sorted);
    end
end

assert(numel(unique(cellfun(@(x) numel(x.x), elevation_kernels)))==1)
assert(numel(unique(cellfun(@(x) numel(x.y), elevation_kernels)))==1)
assert(numel(unique(cellfun(@(x) numel(x.z), elevation_kernels)))==1)

%% Compute flightspeed

RMIN_M  = 5000;
RMAX_M  = 37500;
ZSTEP_M = diff(elevation_kernels{1}.z(1:2));
ZMAX_M  = ZSTEP_M*numel(elevation_kernels{1}.z);

parfor i = 1:numel(radar_scans)
    if ~isempty(radar_scans{i})
        [ ~, ~, u, v, rmse, nll, cov] = epvvp(radar_scans{i}, ZSTEP_M, RMIN_M, RMAX_M, ZMAX_M);
        radar_scans{i}.uv = [u v]; % write speed (m/s)
    end
end

end
