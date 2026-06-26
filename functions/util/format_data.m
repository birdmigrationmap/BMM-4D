function [obs_values, nonzero_obs_idx, azimuth_deg, elevation_deg, range_m, beam_height_m] = format_data(obs_vector, elevation_kernels, info)
%FORMAT_DATA Extract nonzero observations and optional geometry metadata.
%   [OBS_VALUES, NONZERO_OBS_IDX] = FORMAT_DATA(OBS_VECTOR) returns the
%   nonzero entries from one sparse observation vector. Additional outputs
%   require ELEVATION_KERNELS and INFO.

% Retrieve the nonzero observation indices and values.
nonzero_obs_idx = find(obs_vector(:,1));
obs_values = full(obs_vector(nonzero_obs_idx,1));

if nargout > 2
    
    % Find the polar coordinates of the nonzero observations.
    [id_range, id_az, ~, id_kernel] = ind2sub(info.obs_grid_size, nonzero_obs_idx);
        
    range_m = nan(size(id_range));
    azimuth_deg = nan(size(id_range));
    for i_f = unique(id_kernel)'
        range_m(id_kernel == i_f) = elevation_kernels{i_f}.range(id_range(id_kernel == i_f));
        azimuth_deg(id_kernel == i_f) = elevation_kernels{i_f}.az(id_az(id_kernel == i_f));
    end

    kernel_elevation_deg = cellfun(@(x) x.elev, elevation_kernels);
    elevation_deg = kernel_elevation_deg(id_kernel);
end

if nargout > 5
    [~, beam_height_m] = slant2ground(range_m(:), elevation_deg(:));
end
end
