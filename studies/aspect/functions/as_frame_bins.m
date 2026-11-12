function [bins, volumes] = as_frame_bins(manifest, i_frame, fix_dir, out_file)
%AS_FRAME_BINS Bin all volumes of one frame on a common 3D grid (dev).
%   [BINS, VOLUMES] = AS_FRAME_BINS(MANIFEST, I_FRAME, FIX_DIR, OUT_FILE)
%   reads every volume of frame I_FRAME of the us_east manifest three ways:
%   - raw: the FIX_DIR volume (symmetric aspect correction) with the
%     correction undone (as_read_volume), i.e. no aspect correction;
%   - ell: the production volume (manifest.volume_files, elliptical shape
%     correction); its ratio to raw is constant along each ray, so the
%     ray median of production / raw is the factor of every gate of the ray
%     (gates removed by the old beam-blockage lookup keep a factor);
%   - fix: the FIX_DIR volume as saved (per-volume symmetric fit).
%   Gates (beam centre below 3500 m above sea level) are moved to the frame
%   time with the radar's own velocity profile at their height, then binned
%   on 4 km x 4 km x 100 m cells (equirectangular projection at 37 N, 85 W;
%   heights above sea level). BINS: one row per radar and cell with the
%   mean raw / ell / fix density, mean azimuth (vector mean), range,
%   elevation, height and the radar's velocity at that height.
%   VOLUMES: one row per volume with the velocity profile, the saved
%   symmetric fit, and azimuthal sector sums [36 sectors x 12 sweeps x
%   5 range bands] of raw, ell and fix density (for single-radar harmonics
%   and wedges).

frame_time = manifest.times(i_frame);
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
for i_radar = 1:numel(manifest.radars)
    prod_file = manifest.volume_files(i_radar, i_frame);
    if ismissing(prod_file) || prod_file == ""
        continue
    end
    [~, base] = fileparts(prod_file);
    fix_file = fullfile(fix_dir, base + ".mat");
    if ~isfile(fix_file) || ~isfile(prod_file)
        continue
    end
    [raw, gates] = as_read_volume(fix_file, 3500);
    fixed = load(fix_file, "volume").volume;
    prod = load(prod_file, "volume").volume;

    % Gate factors of the production (ell) and per-volume symmetric (fix) corrections
    f_ell = nan(height(gates), 1);
    f_fix = nan(height(gates), 1);
    for i_elev = 1:numel(raw.dens.sweeps)
        rows = gates.sweep == i_elev;
        if ~any(rows)
            continue
        end
        raw_data = raw.dens.sweeps(i_elev).data;
        az_list = raw.dens.sweeps(i_elev).az(:);
        [~, i_az] = ismember(gates.az_deg(rows), az_list);
        [~, i_range] = ismember(gates.range_m(rows), raw.dens.sweeps(i_elev).range(:));
        index = sub2ind(size(raw_data), i_range, i_az);
        f_fix(rows) = fixed.dens.sweeps(i_elev).data(index) ./ raw_data(index);
        ray_factor = nan(1, numel(az_list));
        if i_elev <= numel(prod.dens.sweeps) && isequal(size(prod.dens.sweeps(i_elev).data), size(raw_data))
            ratio = prod.dens.sweeps(i_elev).data ./ raw_data;
            ratio(~isfinite(ratio) | raw_data <= 0) = nan;
            ray_factor = median(ratio, 1, "omitnan");
            % Rays with no common gate: nearest ray with a factor
            ok = isfinite(ray_factor);
            if any(ok) && ~all(ok)
                az_ok = az_list(ok);
                [~, nearest] = min(abs(mod(az_list(~ok)' - az_ok + 180, 360) - 180), [], 1);
                filled = ray_factor(ok);
                ray_factor(~ok) = filled(nearest);
            end
        end
        f_ell(rows) = ray_factor(i_az);
    end

    % Velocity of the radar's profile at each gate's height (100 m bins above the antenna)
    uv = raw.uv;
    i_layer = min(size(uv, 1), max(1, floor(gates.height_agl_m / 100) + 1));
    u = uv(i_layer, 1);
    v = uv(i_layer, 2);
    dt_s = seconds(frame_time - raw.datetime);
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
    % Mean of the finite values in each cell (sums and counts: fast)
    mean_of = @(values) accumarray(id, fill_zero(values(keep))) ./ accumarray(id, isfinite(values(keep)));
    d_raw = gates.density;
    bin_rows = table(repmat(i_frame, numel(key_u), 1), repmat(i_radar, numel(key_u), 1), key_u, n, ...
        mean_of(d_raw), mean_of(d_raw .* f_ell), mean_of(d_raw .* f_fix), ...
        mod(atan2d(accumarray(id, sind(gates.az_deg(keep))), accumarray(id, cosd(gates.az_deg(keep)))), 360), ...
        mean_of(gates.range_m), mean_of(gates.elev), mean_of(gates.height_asl_m), mean_of(u), mean_of(v), ...
        mean_of(gates.sd .^ 2), ...
        VariableNames=["frame" "radar" "key" "n" "raw" "ell" "fix" "az" "range_m" "elev" "height_asl_m" ...
        "u" "v" "sd2"]);
    bin_parts{end + 1} = bin_rows; %#ok<AGROW>

    % Azimuthal sector sums per sweep and range band
    sector = min(n_sector, floor(gates.az_deg / (360 / n_sector)) + 1);
    band = discretize(gates.range_m, range_edges);
    sweep = min(gates.sweep, n_sweep_max);
    ok = isfinite(band);
    sub = [sector(ok) sweep(ok) band(ok)];
    sz = [n_sector n_sweep_max numel(range_edges) - 1];
    d = gates.density(ok);
    sums = struct("n", accumarray(sub, 1, sz), "raw", accumarray(sub, d, sz), ...
        "ell", accumarray(sub, fill_zero(d .* f_ell(ok)), sz), ...
        "fix", accumarray(sub, fill_zero(d .* f_fix(ok)), sz), ...
        "height", accumarray(sub, gates.height_agl_m(ok), sz));
    elev = nan(1, n_sweep_max);
    n_el = min(n_sweep_max, numel(raw.dens.sweeps));
    elev(1:n_el) = [raw.dens.sweeps(1:n_el).elev];
    aspect = struct("track_deg", nan, "amplitude", 0);
    if isfield(fixed, "aspect")
        aspect = fixed.aspect;
    end
    volume_rows{end + 1} = struct("frame", i_frame, "radar", i_radar, "station", raw.station, ...
        "lat", raw.lat, "lon", raw.lon, "height", raw.height, "datetime", raw.datetime, "dt_s", dt_s, ...
        "uv", raw.uv, "fix_track_deg", aspect.track_deg, "fix_amplitude", aspect.amplitude, ...
        "elev", elev, "sectors", sums); %#ok<AGROW>
end
bins = vertcat(bin_parts{:});
volumes = [volume_rows{:}];

% Keep the cells seen by at least two radars (overlap), compact precision
[~, ~, key_id] = unique(bins.key);
n_radar_in_cell = accumarray(key_id, 1);
bins = bins(n_radar_in_cell(key_id) >= 2, :);
for name = ["raw" "ell" "fix" "az" "range_m" "elev" "height_asl_m" "u" "v" "sd2"]
    bins.(name) = single(bins.(name));
end
bins.frame = uint8(bins.frame);
bins.radar = uint8(bins.radar);
bins.n = uint8(bins.n);
if nargin > 3 && strlength(out_file) > 0
    save(out_file, "bins", "volumes")
end
end

function values = fill_zero(values)
values(~isfinite(values)) = 0;
end
