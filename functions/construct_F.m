function [F, obs_matrix, info] = construct_F(radar_scans, elevation_kernels, g)
%CONSTRUCT_F Build the sparse forward operator and sparse observation matrix.
%   [F, OBS_MATRIX, INFO] = CONSTRUCT_F(RADAR_SCANS, ELEVATION_KERNELS, GRID)
%   combines the precomputed elevation kernels with the processed radar
%   scans on GRID. The result is a sparse forward operator F and a sparse
%   observation matrix OBS_MATRIX with one column per time step.

%% Extract radar scan 
% A radar scan contains all different scan (elevation) possible for all
% timestep for each radar
bb = load("data/reference/beam_blockage/beam_blockage.mat");
radar_sites = cell(height(radar_scans),1);
for i_r = 1:height(radar_scans)
    site_kernel_idx = cellfun(@(x) [x.kernel_idx], radar_scans(i_r,:), 'UniformOutput', false);
    radar_sites{i_r}.kernel_idx = unique(vertcat(site_kernel_idx{:}));
    radar_sites{i_r}.lon = radar_scans{i_r,1}.lon;
    radar_sites{i_r}.lat = radar_scans{i_r,1}.lat;
    radar_sites{i_r}.height = radar_scans{i_r,1}.height;
    radar_sites{i_r}.station = radar_scans{i_r,1}.station;

    % Find location on g
    [~,lon_id] = min(abs(g.lon-radar_sites{i_r}.lon));
    [~,lat_id] = min(abs(g.lat-radar_sites{i_r}.lat));
    [~,z_id]  =  min(squeeze(abs(g.z-radar_sites{i_r}.height)));

    % Convert the index of position in x,y,z on the radar volume to the
    % position (index) of the multi-radar g. This is done by offest
    % the lat lon sub by the position of the radar. 
    % get the position of the 0,0 of the coarse g
    radar_sites{i_r}.lon_id = lon_id - ceil(elevation_kernels{1}.sz_coarse_grid(1)/2);
    radar_sites{i_r}.lat_id = lat_id - ceil(elevation_kernels{1}.sz_coarse_grid(2)/2);
    % Synthetic scripts set height = 0 to force a ground-referenced radar.
    if radar_sites{i_r}.height == 0
        radar_sites{i_r}.z_id = 0;
        radar_sites{i_r}.z0_m = 0;
    else
        % Use the nearest analysis-g height as the canonical radar level.
        radar_sites{i_r}.z_id = z_id;
        radar_sites{i_r}.z0_m = g.z(z_id);
    end
   

    % Check that the radar is within the g
%     id_lon = radar_sites{i_r}.lon_id - ceil(elevation_kernels.sz_coarse_grid(1)/2);
%     id_lat = radar_sites{i_r}.lat_id - ceil(elevation_kernels.sz_coarse_grid(2)/2);
%     assert(id_lon>=0 & (id_lon+elevation_kernels.sz_coarse_grid(1))<=numel(g.lon))
%     assert(id_lat>=0 & (id_lat+elevation_kernels.sz_coarse_grid(2))<=numel(g.lat))

    % Beam blockage
    % Build a map of lower beam altitude reach based on the elevation from the beam blockage of this radar and corrected
    % for the position of the radar on the current g
    bbr = bb.bb(:,:,radar_scans{i_r,1}.station==bb.bb_radar);


    [AZ, R] = meshgrid(bb.bb_az, bb.bb_range);
    [x,y,z] = sph2cart(deg2rad(AZ), deg2rad(double(bbr)), R);
    bbr_interp = scatteredInterpolant(g.x(lon_id) + x(:), g.y(lat_id) + y(:), radar_sites{i_r}.z0_m + z(:));

    radar_sites{i_r}.bb_alt = bbr_interp({g.x, g.y});
    % figure; imagesc(g.x, g.y, radar_sites{i_r}.bb_alt); colorbar;
end

%% 
% Take the larger size of az and range to have everything on a regular g
kernel_obs_size = cellfun(@(x) x.sz_raz, elevation_kernels, 'UniformOutput', false);
obs_grid_size = [max(vertcat(kernel_obs_size{:})) numel(radar_sites) numel(elevation_kernels)];

% Compute a index of possible value
obs_block_offset = reshape(cumsum(repmat(prod(obs_grid_size(1:2)), prod(obs_grid_size(3:4)),1)), obs_grid_size(3), obs_grid_size(4)) - prod(obs_grid_size(1:2));

% Check if obs_block_offset is correct [id_range, id_az, id_r, id_kernel]=ind2sub(obs_grid_size, obs_block_offset+1);

% We also want to convert the index of F_raz_id into the new g
% (e.g. 900x360 -> 1800 x 720)
kernel_obs_raz_idx = cell(size(elevation_kernels));
kernel_grid_x = cell(size(elevation_kernels));
kernel_grid_y = cell(size(elevation_kernels));
kernel_grid_z = cell(size(elevation_kernels));
for i_f=1:numel(elevation_kernels)  
    [row,col] = ind2sub(elevation_kernels{i_f}.sz_raz, elevation_kernels{i_f}.F_raz_id);
    kernel_obs_raz_idx{i_f} = sub2ind(obs_grid_size(1:2), row, col);
    kernel_grid_x{i_f} = elevation_kernels{i_f}.F_gr_x;
    kernel_grid_y{i_f} = elevation_kernels{i_f}.F_gr_y;
    kernel_grid_z{i_f} = elevation_kernels{i_f}.F_gr_z;
end

%% Compute the observation matrix
obs_entries = cell(numel(radar_sites), numel(elevation_kernels));
obs_presence = cell(numel(radar_sites), numel(elevation_kernels));
for i_r = 1:numel(radar_sites) % through the radars
    for i_f = 1:numel(elevation_kernels) % through the unique scan
        obs_entries_time = cell(width(radar_scans), 1);
        obs_presence_i = false(prod(obs_grid_size(1:2)), 1);
        for i_t = 1:width(radar_scans) % time 
            % Find the index in the radar data of the particular scan.
            i_elev = find(radar_scans{i_r,i_t}.kernel_idx==i_f);

            if isempty(i_elev) % Does not exist in this radar/time
                continue
            end
            
            % Check that the first few range azimuth are similar (i.e., similar
            % start and same resolution) -> assess as the difference to the
            % resolution.
            assert(abs(elevation_kernels{i_f}.elev-radar_scans{i_r,i_t}.dens.sweeps(i_elev).elev)<0.1)
            assert(all(abs(elevation_kernels{i_f}.range(1:10)-radar_scans{i_r,i_t}.dens.sweeps(i_elev).range(1:10))<diff(elevation_kernels{i_f}.range(1:2))))
            assert(all(abs(elevation_kernels{i_f}.az(1:10)-radar_scans{i_r,i_t}.dens.sweeps(i_elev).az(1:10))<diff(elevation_kernels{i_f}.az(1:2))))
            
            sweep_data = radar_scans{i_r,i_t}.dens.sweeps(i_elev).data + eps; % trick to use sparse for zero density.
            present_obs_id = find(~isnan(sweep_data));
            [i_range, i_az] = ind2sub(size(sweep_data), present_obs_id);
            assert(all(i_range <= obs_grid_size(1)))
            assert(all(i_az <= obs_grid_size(2)))
            i_raz = sub2ind(obs_grid_size(1:2), i_range, i_az);
            obs_presence_i(i_raz) = true;
            obs_entries_time{i_t} = [obs_block_offset(i_r,i_f) + i_raz repmat(i_t, numel(i_raz), 1) sweep_data(present_obs_id)];
        end

        obs_entries{i_r,i_f} = vertcat(obs_entries_time{:});
        obs_presence{i_r,i_f} = find(obs_presence_i);
    end
end

obs_triplets = vertcat(obs_entries{:});
obs_matrix = sparse(obs_triplets(:,1), obs_triplets(:,2), obs_triplets(:,3), prod(obs_grid_size), width(radar_scans));


%% Compute F
% Reuse the precomputed kernel triplets for every radar. Only the
% observation-row offset and g-cell offset change by radar site.

F_rows = cell(numel(radar_sites), numel(elevation_kernels));
F_cols = cell(numel(radar_sites), numel(elevation_kernels));
F_vals = cell(numel(radar_sites), numel(elevation_kernels));
coverage_fraction = ones(numel(radar_sites), numel(elevation_kernels));

for i_r = 1:numel(radar_sites) % through the radars
    for i_f = 1:numel(elevation_kernels) % through the unique scan
        id_pres_data = ismember(kernel_obs_raz_idx{i_f}, obs_presence{i_r,i_f});

        grid_x_id = kernel_grid_x{i_f}(id_pres_data) + radar_sites{i_r}.lon_id;
        grid_y_id = kernel_grid_y{i_f}(id_pres_data) + radar_sites{i_r}.lat_id;
        grid_z_id = kernel_grid_z{i_f}(id_pres_data) + radar_sites{i_r}.z_id;
        id_F_keep = grid_x_id >= 1 & grid_x_id <= g.sz(1) & ...
            grid_y_id >= 1 & grid_y_id <= g.sz(2) & ...
            grid_z_id >= 1 & grid_z_id <= g.sz(3);
        if ~isempty(id_F_keep)
            coverage_fraction(i_r, i_f) = mean(id_F_keep);
        end

        F_rows{i_r,i_f} = obs_block_offset(i_r,i_f) + kernel_obs_raz_idx{i_f}(id_pres_data);
        F_rows{i_r,i_f} = F_rows{i_r,i_f}(id_F_keep);
        F_cols{i_r,i_f} = sub2ind(g.sz, grid_x_id(id_F_keep), grid_y_id(id_F_keep), grid_z_id(id_F_keep));
        F_vals{i_r,i_f} = double(elevation_kernels{i_f}.Fval(id_pres_data));
        F_vals{i_r,i_f} = F_vals{i_r,i_f}(id_F_keep);

%         F = sparse(F_rows{i_r,i_f}, F_cols{i_r,i_f}, F_vals{i_r,i_f}, prod(obs_grid_size), prod(g.sz));
%         tmp = reshape(full(sum(F,1)), g.sz);
%         figure;imagesc(log(sum(tmp,3)))
    end
end

%% 
F = sparse(vertcat(F_rows{:}), vertcat(F_cols{:}), vertcat(F_vals{:}), prod(obs_grid_size), prod(g.sz));

%%
% norm = sum(F,2);
% % tmp = find(Fsp>0);
% % [row,~] = ind2sub(size(Fsp),tmp);
% % Fsp(tmp) = Fsp(tmp)./ norm(row);
% [row, col] = find(F>0);
% tmp = F(F>0) ./ norm(row);
% F = sparse(row,col,tmp,size(F,1),size(F,2));

%% 
info.obs_grid_size = obs_grid_size;
info.obs_block_offset = obs_block_offset;
info.coverage_fraction = coverage_fraction;
info.has_partial_coverage = coverage_fraction < 1;
end
