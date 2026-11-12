function bins = as_rain_bins(raw_files, radar_table, ~, i_frame, out_file)
%AS_RAIN_BINS Overlap bins of rain reflectivity, for radar patterns without birds (dev).
%   BINS = AS_RAIN_BINS(RAW_FILES, RADAR_TABLE, FRAME_TIME, I_FRAME, OUT_FILE)
%   reads the raw Level-II reflectivity of each radar (RAW_FILES [n_radar 1],
%   "" where missing; row = radar id, sites in RADAR_TABLE) and bins it on
%   the cells of as_night_bins (4 km x 4 km x 100 m, heights above sea
%   level). Rain has no aspect, so the azimuthal patterns of the log ratio
%   of two radars in the same cell are radar patterns (blockage,
%   calibration by sector).
%   Gates kept: sweeps up to 1.6 deg, 20-120 km, beam centre below 1200 m
%   above the antenna (below the melting layer on spring and autumn days),
%   15-45 dBZ (stratiform rain; no clear air, no hail cores).
%   No advection to the frame time (the volume nearest the frame).
%   BINS: as as_night_bins with the mean linear reflectivity in "raw"
%   (mm^6 m^-3), u = 0 and v = 1 (a dummy track to the north, so psi is the
%   compass azimuth in as_overlap_design).

bin_size_m = 4000;
dz_m = 100;
earth_radius = 6371200;
lat0 = 37;
lon0 = -85;
bin_parts = {};
for i_radar = 1:numel(raw_files)
    file = raw_files(i_radar);
    if ismissing(file) || file == "" || ~isfile(file) || ~isfinite(radar_table.lat(i_radar))
        continue
    end
    try
        scan = bmm4d.radar.rsl2mat_wrap(file, struct("max_elev", 1.6));
    catch
        continue
    end
    if ~isfield(scan, "dz") || isempty(scan.dz.sweeps)
        continue
    end
    % One reflectivity sweep per elevation (split cuts repeat it: keep the longest range)
    sweeps = scan.dz.sweeps;
    elev_all = round([sweeps.elev], 1);
    parts = {};
    for e = unique(elev_all)
        candidates = find(elev_all == e);
        best_n = -1;
        for k = candidates
            [data_k, range_k, az_k] = sweep2mat(sweeps(k));
            if numel(range_k) > best_n
                [data, range_m, az, best_n] = deal(data_k, range_k, az_k, numel(range_k));
            end
        end
        [range_grid, az_grid] = ndgrid(range_m(:), az(:));
        [ground_m, height_m] = slant2ground(range_grid, e);
        keep = isfinite(data) & data >= 15 & data <= 45 & range_grid >= 20000 & range_grid <= 120000 & ...
            height_m <= 1200;
        parts{end + 1} = [data(keep), range_grid(keep), az_grid(keep), ground_m(keep), height_m(keep), ...
            repmat(e, nnz(keep), 1)]; %#ok<AGROW>
    end
    gates = vertcat(parts{:});
    if isempty(gates)
        continue
    end
    z_lin = 10 .^ (gates(:, 1) / 10);
    az_deg = gates(:, 3);
    % Ground position on the sphere (as as_read_volume)
    delta = gates(:, 4) / earth_radius;
    lat1 = deg2rad(radar_table.lat(i_radar));
    lat2 = asin(sin(lat1) .* cos(delta) + cos(lat1) .* sin(delta) .* cosd(az_deg));
    lon2 = deg2rad(radar_table.lon(i_radar)) + atan2(sind(az_deg) .* sin(delta) .* cos(lat1), cos(delta) - sin(lat1) .* sin(lat2));
    height_asl_m = gates(:, 5) + radar_table.height_m(i_radar);
    x = earth_radius * cosd(lat0) * (lon2 - deg2rad(lon0));
    y = earth_radius * (lat2 - deg2rad(lat0));
    key = ((int64(floor(x / bin_size_m)) + 2^20) * 2^21 + (int64(floor(y / bin_size_m)) + 2^20)) * 2^6 + ...
        int64(floor(height_asl_m / dz_m));
    [key_u, ~, id] = unique(key);
    n = accumarray(id, 1);
    mean_of = @(values) accumarray(id, values) ./ n;
    bin_parts{end + 1} = table(repmat(i_frame, numel(key_u), 1), repmat(i_radar, numel(key_u), 1), key_u, n, ...
        mean_of(z_lin), mod(atan2d(accumarray(id, sind(az_deg)), accumarray(id, cosd(az_deg))), 360), ...
        mean_of(gates(:, 2)), mean_of(gates(:, 6)), mean_of(height_asl_m), zeros(numel(key_u), 1), ...
        ones(numel(key_u), 1), VariableNames=["frame" "radar" "key" "n" "raw" "az" "range_m" "elev" ...
        "height_asl_m" "u" "v"]); %#ok<AGROW>
end
if isempty(bin_parts)
    bins = table();
else
    bins = vertcat(bin_parts{:});
    [~, ~, key_id] = unique(bins.key);
    n_radar_in_cell = accumarray(key_id, 1);
    bins = bins(n_radar_in_cell(key_id) >= 2, :);
    for name = ["raw" "az" "range_m" "elev" "height_asl_m" "u" "v"]
        bins.(name) = single(bins.(name));
    end
    bins.frame = uint8(bins.frame);
    bins.radar = uint8(bins.radar);
    bins.n = uint8(min(bins.n, 255));
end
volumes = struct([]);
if nargin > 4 && strlength(out_file) > 0
    save(out_file, "bins", "volumes")
end
end
