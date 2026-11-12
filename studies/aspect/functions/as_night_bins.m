function [bins, volumes] = as_night_bins(volume_files, frame_time, i_frame, out_file)
%AS_NIGHT_BINS Bin the cleaned volumes of one frame on a common 3D grid (dev).
%   [BINS, VOLUMES] = AS_NIGHT_BINS(VOLUME_FILES, FRAME_TIME, I_FRAME, OUT_FILE)
%   is as_frame_bins for any night cleaned without shape or aspect
%   correction (bmm4d.radar.clean_volume with read option "shape" false):
%   one density version only ("raw"). VOLUME_FILES [n_radar 1] holds the
%   cleaned volume of each radar for this frame ("" where missing); the
%   radar id is its row, so it must follow one radar list for all nights.
%   Gates (beam centre below 3500 m above sea level) are moved to
%   FRAME_TIME with the radar's own velocity profile at their height, then
%   binned on 4 km x 4 km x 100 m cells (equirectangular projection at
%   37 N, 85 W, as as_frame_bins; heights above sea level).
%   BINS: one row per radar and overlap cell (seen by two radars or more)
%   with the mean raw density, vector-mean azimuth, range, elevation,
%   height, the radar's velocity at that height and the mean estimator
%   variance. VOLUMES: one struct per volume with the velocity profile and
%   azimuthal sector sums [36 sectors x 12 sweeps x 5 range bands] of the
%   raw density (as_cell_harmonics, as_time_fit).
%   OUT_FILE ("" to skip) receives bins and volumes.

bin_size_m = 4000;
dz_m = 100;
earth_radius = 6371200;
lat0 = 37;
lon0 = -85;
range_edges = [5 10 30 60 100 150] * 1000;
n_sector = 36;
n_sweep_max = 12;

bin_parts = {};
volume_rows = {};
for i_radar = 1:numel(volume_files)
    file = volume_files(i_radar);
    if ismissing(file) || file == "" || ~isfile(file)
        continue
    end
    [vol, gates] = as_read_volume(file, 3500);
    if isempty(gates) || height(gates) == 0
        continue
    end

    % Velocity of the radar's profile at each gate's height (100 m layers above the antenna)
    uv = vol.uv;
    i_layer = min(size(uv, 1), max(1, floor(gates.height_agl_m / 100) + 1));
    u = uv(i_layer, 1);
    v = uv(i_layer, 2);
    dt_s = seconds(frame_time - vol.datetime);
    u_ok = u;
    v_ok = v;
    u_ok(~isfinite(u_ok)) = 0;
    v_ok(~isfinite(v_ok)) = 0;

    % Position at the frame time, then the cell
    x = earth_radius * cosd(lat0) * deg2rad(gates.lon - lon0) + u_ok * dt_s;
    y = earth_radius * deg2rad(gates.lat - lat0) + v_ok * dt_s;
    ix = int64(floor(x / bin_size_m));
    iy = int64(floor(y / bin_size_m));
    iz = int64(floor(gates.height_asl_m / dz_m));
    key = ((ix + 2^20) * 2^21 + (iy + 2^20)) * 2^6 + iz;
    keep = gates.density >= 0;
    [key_u, ~, id] = unique(key(keep));
    n = accumarray(id, 1);
    mean_of = @(values) accumarray(id, fill_zero(values(keep))) ./ accumarray(id, isfinite(values(keep)));
    bin_rows = table(repmat(i_frame, numel(key_u), 1), repmat(i_radar, numel(key_u), 1), key_u, n, ...
        mean_of(gates.density), ...
        mod(atan2d(accumarray(id, sind(gates.az_deg(keep))), accumarray(id, cosd(gates.az_deg(keep)))), 360), ...
        mean_of(gates.range_m), mean_of(gates.elev), mean_of(gates.height_asl_m), mean_of(u), mean_of(v), ...
        mean_of(gates.sd .^ 2), ...
        VariableNames=["frame" "radar" "key" "n" "raw" "az" "range_m" "elev" "height_asl_m" "u" "v" "sd2"]);
    bin_parts{end + 1} = bin_rows; %#ok<AGROW>

    % Azimuthal sector sums per sweep and range band (single-radar harmonics)
    sector = min(n_sector, floor(gates.az_deg / (360 / n_sector)) + 1);
    band = discretize(gates.range_m, range_edges);
    sweep = min(gates.sweep, n_sweep_max);
    ok = isfinite(band);
    sub = [sector(ok) sweep(ok) band(ok)];
    sz = [n_sector n_sweep_max numel(range_edges) - 1];
    sums = struct("n", accumarray(sub, 1, sz), "raw", accumarray(sub, gates.density(ok), sz), ...
        "height", accumarray(sub, gates.height_agl_m(ok), sz));
    elev = nan(1, n_sweep_max);
    n_el = min(n_sweep_max, numel(vol.dens.sweeps));
    elev(1:n_el) = [vol.dens.sweeps(1:n_el).elev];
    volume_rows{end + 1} = struct("frame", i_frame, "radar", i_radar, "station", vol.station, ...
        "lat", vol.lat, "lon", vol.lon, "height", vol.height, "datetime", vol.datetime, "dt_s", dt_s, ...
        "uv", vol.uv, "elev", elev, "sectors", sums); %#ok<AGROW>
end
if isempty(bin_parts)
    bins = table();
    volumes = struct([]);
    return
end
bins = vertcat(bin_parts{:});
volumes = [volume_rows{:}];

% Keep the cells seen by at least two radars (overlap), compact precision
[~, ~, key_id] = unique(bins.key);
n_radar_in_cell = accumarray(key_id, 1);
bins = bins(n_radar_in_cell(key_id) >= 2, :);
for name = ["raw" "az" "range_m" "elev" "height_asl_m" "u" "v" "sd2"]
    bins.(name) = single(bins.(name));
end
bins.frame = uint8(bins.frame);
bins.radar = uint8(bins.radar);
bins.n = uint8(min(bins.n, 255));
if nargin > 3 && strlength(out_file) > 0
    save(out_file, "bins", "volumes")
end
end

function values = fill_zero(values)
values(~isfinite(values)) = 0;
end
