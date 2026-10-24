function status = build_radar_cache(station, radar_files, available_id, n_time, g, grid_key, ...
    read_options, kernel_dir, min_time_fraction, min_row_mass, cache_file)
%BUILD_RADAR_CACHE Clean one radar's night and save its compact observation rows.
%   Runs serially so that UN_B can process whole radars in parallel: clean
%   each volume, attach kernels and velocity profiles, build this radar's
%   forward-operator rows on the shared grid, and save them to CACHE_FILE.
%   STATUS is a one-line message for the log.

radar_tic = tic;
if numel(available_id) < min_time_fraction * n_time
    status = sprintf("%s | skipped, %d/%d times available", station, numel(available_id), n_time);
    return
end

% Clean each volume; a corrupt file only drops its own time.
scans = cell(1, numel(available_id));
for i_f = 1:numel(radar_files)
    try
        scan_cell = bmm4d.radar.read_rd(radar_files(i_f), read_options{:});
        scans{i_f} = scan_cell{1};
    catch err
        warning("UN_B:ReadFailed", "%s: %s", radar_files(i_f), err.message)
    end
end
read_ok = ~cellfun(@isempty, scans);
scans = scans(read_ok);
available_id = available_id(read_ok);
if numel(available_id) < min_time_fraction * n_time
    status = sprintf("%s | skipped, %d/%d times readable", station, numel(available_id), n_time);
    return
end

% Attach kernels and velocity profiles (needs the raw velocity sweeps),
% then keep only what the forward model uses before building the rows.
% load_F_elev fails when a scan pattern has no precomputed kernel.
try
    [elevation_kernels, scans] = bmm4d.forward.load_F_elev(kernel_dir, scans);
    for i_s = 1:numel(scans)
        scans{i_s} = struct( ...
            "station", scans{i_s}.station, "lat", scans{i_s}.lat, "lon", scans{i_s}.lon, ...
            "height", scans{i_s}.height, "datetime", scans{i_s}.datetime, ...
            "dens", scans{i_s}.dens, "kernel_idx", scans{i_s}.kernel_idx, "uv", scans{i_s}.uv);
    end
    [F, obs_matrix, info] = bmm4d.forward.construct_F(scans, elevation_kernels, g);
catch err
    status = sprintf("%s | skipped, forward model failed: %s", station, err.message);
    return
end
radar_obs = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info, min_row_mass=min_row_mass);
clear F obs_matrix info elevation_kernels

% Place observed columns at their requested times; missing times stay zero.
radar_F = radar_obs.F;
[row_id, col_id, obs_value] = find(sparse(double(radar_obs.values)));
radar_values = sparse(row_id, available_id(col_id), obs_value, size(radar_F, 1), n_time);

% Velocity profiles (epvvp bins above the antenna) at all requested
% times, filled linearly in time where a volume is missing.
% radar_uv_on_grid maps them to the grid; NaN heights get no advection.
radar_uv = nan([size(scans{1}.uv) n_time]);
for i_s = 1:numel(scans)
    radar_uv(:, :, available_id(i_s)) = scans{i_s}.uv;
end
radar_uv = fillmissing(radar_uv, "linear", 3, EndValues="nearest");

radar_site = struct("station", station, "lat", scans{1}.lat, "lon", scans{1}.lon, "height", scans{1}.height);
radar_time_id = available_id;
save(cache_file, "radar_F", "radar_values", "radar_uv", "radar_site", "radar_time_id", "grid_key", "radar_files")
status = sprintf("%s | %d/%d times, %.0f s", station, numel(available_id), n_time, toc(radar_tic));
end
