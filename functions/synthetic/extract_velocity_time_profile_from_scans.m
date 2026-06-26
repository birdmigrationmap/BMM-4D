function uv_profiles_time = extract_velocity_time_profile_from_scans(radar_scans)

if isvector(radar_scans)
    radar_scans_time = radar_scans;
else
    radar_scans_time = radar_scans(1, :);
end

n_time = numel(radar_scans_time);
nz = size(radar_scans_time{1}.uv, 1);
uv_profiles_time = nan(nz, 2, n_time);
for i_t = 1:n_time
    uv_profiles_time(:, :, i_t) = radar_scans_time{i_t}.uv;
end
end
