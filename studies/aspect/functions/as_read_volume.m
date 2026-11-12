function [volume, gates] = as_read_volume(file, max_height_asl_m)
%AS_READ_VOLUME Cleaned volume without aspect correction, and its gates (dev).
%   [VOLUME, GATES] = AS_READ_VOLUME(FILE, MAX_HEIGHT_ASL_M) loads a volume
%   saved by bmm4d.radar.clean_volume (or ar_clean_volume / as_clean_volume).
%   When the volume carries an "aspect" struct (symmetric correction
%   T(az) = 1 + a cos 2(az - track - 90 deg), each sweep rescaled to keep
%   its mean), the correction is undone: data .* T is the raw sweep times
%   the sweep's rescale factor, which is recovered as mean(D) / mean(D ./ T)
%   on the thinned gates. The shared aspect function of as_clean_volume
%   (aspect struct with b2_per_km) is undone by multiplying by its template.
%   GATES: one row per finite thinned gate below MAX_HEIGHT_ASL_M (beam
%   centre, 4/3 earth): sweep, elev, range_m, az_deg, ground_m,
%   height_asl_m, height_agl_m (above the antenna), lat, lon, density, sd.

loaded = load(file, "volume");
volume = loaded.volume;
sweeps = volume.dens.sweeps;
if isfield(volume, "aspect") && isfield(volume.aspect, "amplitude") && volume.aspect.amplitude > 0
    side = volume.aspect.track_deg + 90;
    for i_elev = 1:numel(sweeps)
        template = 1 + volume.aspect.amplitude * cosd(2 * (sweeps(i_elev).az(:)' - side));
        rescaled = sweeps(i_elev).data .* template;
        % rescaled ./ template is the stored (corrected) sweep
        scale = mean(rescaled(:), "omitnan") / mean(sweeps(i_elev).data(:), "omitnan");
        if ~isfinite(scale) || scale <= 0
            scale = 1;
        end
        sweeps(i_elev).data = rescaled / scale;
        if isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
            sweeps(i_elev).sd = sweeps(i_elev).sd .* template / scale;
        end
    end
    volume.dens.sweeps = sweeps;
elseif isfield(volume, "aspect") && isstruct(volume.aspect) && isfield(volume.aspect, "b2_per_km")
    % Shared aspect function of as_clean_volume: multiply back by the template
    for i_elev = 1:numel(sweeps)
        [~, height_m] = slant2ground(sweeps(i_elev).range(:), sweeps(i_elev).elev);
        [az, height] = ndgrid(sweeps(i_elev).az(:)', height_m);
        template = as_aspect_template(az', height', volume.uv, volume.aspect);
        sweeps(i_elev).data = sweeps(i_elev).data .* template;
        if isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
            sweeps(i_elev).sd = sweeps(i_elev).sd .* template;
        end
    end
    volume.dens.sweeps = sweeps;
end
if nargout < 2
    return
end

% Gate list: beam centre position of each finite thinned gate
parts = cell(numel(sweeps), 1);
for i_elev = 1:numel(sweeps)
    [i_range, i_az] = find(isfinite(sweeps(i_elev).data));
    range_m = sweeps(i_elev).range(i_range);
    az_deg = sweeps(i_elev).az(i_az);
    [ground_m, height_m] = slant2ground(range_m(:), sweeps(i_elev).elev);
    keep = height_m + volume.height <= max_height_asl_m;
    index = sub2ind(size(sweeps(i_elev).data), i_range(keep), i_az(keep));
    sd = nan(nnz(keep), 1);
    if isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
        sd = sweeps(i_elev).sd(index);
    end
    parts{i_elev} = [repmat([i_elev sweeps(i_elev).elev], nnz(keep), 1) range_m(keep) az_deg(keep) ...
        ground_m(keep) height_m(keep) + volume.height height_m(keep) sweeps(i_elev).data(index) sd];
end
values = cell2mat(parts);
gates = array2table(values, VariableNames=["sweep" "elev" "range_m" "az_deg" "ground_m" "height_asl_m" ...
    "height_agl_m" "density" "sd"]);

% Ground position on the sphere
earth_radius = 6371200;
delta = gates.ground_m / earth_radius;
lat1 = deg2rad(volume.lat);
lat2 = asin(sin(lat1) .* cos(delta) + cos(lat1) .* sin(delta) .* cosd(gates.az_deg));
lon2 = deg2rad(volume.lon) + atan2(sind(gates.az_deg) .* sin(delta) .* cos(lat1), cos(delta) - sin(lat1) .* sin(lat2));
gates.lat = rad2deg(lat2);
gates.lon = rad2deg(lon2);
end
