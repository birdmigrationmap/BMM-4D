function [bins, volumes] = as_frame_polar(manifest, i_frame, fix_dir, out_file)
%AS_FRAME_POLAR Polarimetric sample of one frame: Z_H, Z_V, ZDR, rhoHV (dev).
%   [BINS, VOLUMES] = AS_FRAME_POLAR(MANIFEST, I_FRAME, FIX_DIR, OUT_FILE)
%   reads the raw Level-II file of every volume of frame I_FRAME and takes
%   the raw reflectivity Z_H (dBZ), ZDR (dB) and rhoHV at the gates kept by
%   the cleaning (finite thinned gates of the FIX_DIR volume: rain, clutter
%   and blockage already removed). Z_V = Z_H - ZDR.
%   BINS: as as_frame_bins (same cells, moved to the frame time), with the
%   mean linear Z_H and Z_V of each radar and cell (zh, zv) and the raw
%   cleaned density (raw), for the overlap of the two polarizations.
%   VOLUMES: per volume, 36-sector sums at 10-60 km (all sweeps below
%   3500 m) of linear Z_H, Z_V, density, density x rhoHV, density x ZDR
%   (dB), and counts; the velocity profile.

frame_time = manifest.times(i_frame);
bin_size_m = 4000;
dz_m = 100;
earth_radius = 6371200;
lat0 = 37;
lon0 = -85;
bin_parts = {};
volume_rows = {};
for i_radar = 1:numel(manifest.radars)
    prod_file = manifest.volume_files(i_radar, i_frame);
    if ismissing(prod_file) || prod_file == ""
        continue
    end
    [~, base] = fileparts(prod_file);
    fix_file = fullfile(fix_dir, base + ".mat");
    raw_file = manifest.raw_files(i_radar, i_frame);
    if ~isfile(fix_file) || ismissing(raw_file) || ~isfile(raw_file)
        continue
    end
    [cleaned, gates] = as_read_volume(fix_file, 3500);
    try
        scan = bmm4d.radar.rsl2mat_wrap(raw_file, struct('max_elev', 7));
    catch
        continue
    end

    % Raw Z_H, ZDR, rhoHV at the cleaned gates (same polar grid as sweep2mat of dz)
    [zh, zdr, rho] = deal(nan(height(gates), 1));
    dz_elev = [scan.dz.sweeps.elev];
    dr_elev = [scan.dr.sweeps.elev];
    rh_elev = [scan.rh.sweeps.elev];
    for i_elev = 1:numel(cleaned.dens.sweeps)
        rows = find(gates.sweep == i_elev);
        if isempty(rows)
            continue
        end
        az_list = cleaned.dens.sweeps(i_elev).az(:);
        range_list = cleaned.dens.sweeps(i_elev).range(:);
        [~, i_az] = ismember(gates.az_deg(rows), az_list);
        [~, i_range] = ismember(gates.range_m(rows), range_list);
        elev = cleaned.dens.sweeps(i_elev).elev;
        zh(rows) = field_at(scan.dz.sweeps, dz_elev, elev, az_list, range_list(1), i_range, i_az);
        zdr(rows) = field_at(scan.dr.sweeps, dr_elev, elev, az_list, range_list(1), i_range, i_az);
        rho(rows) = field_at(scan.rh.sweeps, rh_elev, elev, az_list, range_list(1), i_range, i_az);
    end
    zh_lin = 10 .^ (zh / 10);
    zv_lin = 10 .^ ((zh - zdr) / 10);
    zh_lin(~isfinite(zh_lin)) = nan;
    zv_lin(~isfinite(zv_lin)) = nan;

    % Cells at the frame time, as as_frame_bins
    i_layer = min(size(cleaned.uv, 1), max(1, floor(gates.height_agl_m / 100) + 1));
    u = cleaned.uv(i_layer, 1);
    v = cleaned.uv(i_layer, 2);
    dt_s = seconds(frame_time - cleaned.datetime);
    shift_u = u;
    shift_v = v;
    shift_u(~isfinite(shift_u)) = 0;
    shift_v(~isfinite(shift_v)) = 0;
    x = earth_radius * cosd(lat0) * deg2rad(gates.lon - lon0) + shift_u * dt_s;
    y = earth_radius * deg2rad(gates.lat - lat0) + shift_v * dt_s;
    key = ((int64(floor(x / bin_size_m)) + 2^20) * 2^21 + (int64(floor(y / bin_size_m)) + 2^20)) * 2^6 + ...
        int64(floor(gates.height_asl_m / dz_m));
    keep = isfinite(zh_lin) & isfinite(zv_lin) & gates.density >= 0;
    [key_u, ~, id] = unique(key(keep));
    n = accumarray(id, 1);
    mean_of = @(values) accumarray(id, values(keep)) ./ n;
    bin_parts{end + 1} = table(repmat(i_frame, numel(key_u), 1), repmat(i_radar, numel(key_u), 1), key_u, n, ...
        mean_of(gates.density), mean_of(zh_lin), mean_of(zv_lin), ...
        mod(atan2d(accumarray(id, sind(gates.az_deg(keep))), accumarray(id, cosd(gates.az_deg(keep)))), 360), ...
        mean_of(gates.range_m), mean_of(gates.elev), mean_of(gates.height_asl_m), mean_of(u), mean_of(v), ...
        VariableNames=["frame" "radar" "key" "n" "raw" "zh" "zv" "az" "range_m" "elev" "height_asl_m" "u" "v"]); %#ok<AGROW>

    % Sector sums at 10-60 km
    near = keep & gates.range_m >= 10000 & gates.range_m <= 60000 & isfinite(rho);
    sector = min(36, floor(gates.az_deg(near) / 10) + 1);
    d = gates.density(near);
    sums = struct("n", accumarray(sector, 1, [36 1]), "zh", accumarray(sector, zh_lin(near), [36 1]), ...
        "zv", accumarray(sector, zv_lin(near), [36 1]), "density", accumarray(sector, d, [36 1]), ...
        "density_rho", accumarray(sector, d .* rho(near), [36 1]), ...
        "density_zdr", accumarray(sector, d .* zdr(near), [36 1]), ...
        "height", accumarray(sector, gates.height_agl_m(near), [36 1]));
    volume_rows{end + 1} = struct("frame", i_frame, "radar", i_radar, "station", cleaned.station, ...
        "datetime", cleaned.datetime, "uv", cleaned.uv, "vcp", double(scan.vcp), "sectors", sums); %#ok<AGROW>
end
bins = vertcat(bin_parts{:});
volumes = [volume_rows{:}];
[~, ~, key_id] = unique(bins.key);
n_radar_in_cell = accumarray(key_id, 1);
bins = bins(n_radar_in_cell(key_id) >= 2, :);
if nargin > 3 && strlength(out_file) > 0
    save(out_file, "bins", "volumes")
end
end

function values = field_at(sweeps, sweep_elev, elev, az_list, range_first, i_range, i_az)
% Values of the raw sweep at the cleaned gates: the sweep at this elevation
% whose sorted azimuths match the cleaned sweep (split cuts repeat elevations)
values = nan(numel(i_range), 1);
candidates = find(abs(sweep_elev - elev) < 0.05);
for k = candidates
    [data, range_m, az] = sweep2mat(sweeps(k));
    if numel(az) ~= numel(az_list) || max(abs(az(:) - az_list)) > 0.01
        continue
    end
    if abs(range_m(1) - range_first) > 1
        continue
    end
    inside = i_range <= size(data, 1);
    values(inside) = data(sub2ind(size(data), i_range(inside), i_az(inside)));
    return
end
end
