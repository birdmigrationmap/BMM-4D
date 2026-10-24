function rows = build_radar_rows(volume_files, n_time, g, kernel_dir, min_row_mass, aggregate, clutter_mask, clutter_max_elev, near_range_m, no_echo_max_range_m)
%BUILD_RADAR_ROWS Forward-operator rows of one radar for all requested times.
%   VOLUME_FILES(1, n_time) are the radar's cleaned volumes ("" if missing).
%   Builds this radar's rows of the forward operator on the shared grid G,
%   keeps rows with at least MIN_ROW_MASS of their beam weight inside the
%   grid, and places the observations at their requested time columns.
%   Velocity profiles of missing times are filled linearly in time.
%   AGGREGATE [fx fy fz] builds the rows on the coarse solve grid
%   (coarsen_grid); the velocity profiles are then averaged over fz bins so
%   their bin spacing matches the coarse layers.
%   CLUTTER_MASK (static_clutter_mask, or []) removes the flagged gates from
%   sweeps below CLUTTER_MAX_ELEV (nearest range/azimuth lookup).
%   Gates closer than NEAR_RANGE_M are dropped: near the radar the bird echo
%   is suppressed (median 6% / 27% / 89% of the 30-50 km density at 0-5 /
%   5-10 / 10-15 km, us_east), which digs a hole above each radar.
%   No-echo gates (zero density) beyond NO_ECHO_MAX_RANGE_M become missing.
%   ROWS is empty when the radar has no usable volume.

rows = [];
time_id = find(volume_files ~= "");
if isempty(time_id)
    return
end
scans = cell(1, numel(time_id));
for k = 1:numel(time_id)
    loaded = load(volume_files(time_id(k)), "volume");
    scans{k} = loaded.volume;
end

% Drop the suppressed near-range gates and the far no-echo gates
for k = 1:numel(scans)
    for i_elev = 1:numel(scans{k}.dens.sweeps)
        sweep_range = scans{k}.dens.sweeps(i_elev).range(:);
        data = scans{k}.dens.sweeps(i_elev).data;
        data(sweep_range < near_range_m, :) = nan;
        data(data == 0 & sweep_range > no_echo_max_range_m) = nan;
        scans{k}.dens.sweeps(i_elev).data = data;
    end
end

% Remove static clutter (wind farms, towers) from the low sweeps
if ~isempty(clutter_mask)
    for k = 1:numel(scans)
        for i_elev = find([scans{k}.dens.sweeps.elev] < clutter_max_elev)
            sweep = scans{k}.dens.sweeps(i_elev);
            range_id = interp1(clutter_mask.range, 1:numel(clutter_mask.range), sweep.range(:), "nearest", NaN);
            [~, az_id] = min(abs(mod(sweep.az(:) - clutter_mask.az(:)' + 180, 360) - 180), [], 2);
            inside = ~isnan(range_id);
            flagged = false(size(sweep.data));
            flagged(inside, :) = clutter_mask.flag(range_id(inside), az_id);
            scans{k}.dens.sweeps(i_elev).data(flagged) = nan;
        end
    end
end
[elevation_kernels, scans, n_dropped] = bmm4d.forward.attach_kernels(scans, kernel_dir);
[F, obs_matrix, info] = bmm4d.forward.construct_F(scans, elevation_kernels, g, aggregate=aggregate);
radar_obs = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info, min_row_mass=min_row_mass);
% Range (m) and elevation (deg) of each kept row, for data weighting
[range_id, ~, ~, kernel_id] = ind2sub(info.obs_grid_size, double(radar_obs.row_ids));
kernel_range = cellfun(@(k) k.range(:), elevation_kernels, UniformOutput=false);
kernel_elev = cellfun(@(k) k.elev, elevation_kernels);
row_range_m = zeros(numel(range_id), 1, "single");
for i_k = 1:numel(elevation_kernels)
    in_kernel = kernel_id == i_k;
    row_range_m(in_kernel) = kernel_range{i_k}(range_id(in_kernel));
end
row_elev_deg = single(kernel_elev(kernel_id));
clear F obs_matrix info elevation_kernels

[row_id, col_id, obs_value] = find(sparse(double(radar_obs.values)));
uv = nan([size(scans{1}.uv) n_time]);
for k = 1:numel(scans)
    uv(:, :, time_id(k)) = scans{k}.uv;
end
uv = fillmissing(uv, "linear", 3, EndValues="nearest");
if aggregate(3) > 1
    % Mean of each block of aggregate(3) height bins
    n_bin = floor(size(uv, 1) / aggregate(3)) * aggregate(3);
    uv = squeeze(mean(reshape(uv(1:n_bin, :, :), aggregate(3), [], size(uv, 2), n_time), 1, "omitnan"));
    uv = reshape(uv, [], 2, n_time);
end

rows = struct( ...
    "F", radar_obs.F, ...
    "row_range_m", row_range_m, ...
    "row_elev_deg", row_elev_deg(:), ...
    "values", sparse(row_id, time_id(col_id), obs_value, size(radar_obs.F, 1), n_time), ...
    "uv", uv, ...
    "site", struct("station", scans{1}.station, "lat", scans{1}.lat, "lon", scans{1}.lon, "height", scans{1}.height), ...
    "time_id", time_id, ...
    "n_dropped_sweeps", n_dropped);
end
