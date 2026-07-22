function radar_scan_sequence = build_single_radar_scan_sequence(template_scan, uv_profile_time, time_offsets, dt)
%BUILD_SINGLE_RADAR_SCAN_SEQUENCE Rebuild compact single-radar synthetic scans.

n_time = numel(time_offsets);
radar_scan_sequence = cell(1, n_time);
uv_profile_time = double(uv_profile_time);
for i_t = 1:n_time
    radar_scan_sequence{i_t} = template_scan;
    radar_scan_sequence{i_t}.datetime = template_scan.datetime + time_offsets(i_t) .* dt;
    radar_scan_sequence{i_t}.uv = uv_profile_time(:, :, i_t);
end
end
