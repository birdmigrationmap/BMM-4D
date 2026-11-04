function status = clean_volume(raw_file, read_options, vp, filter, out_file)
%CLEAN_VOLUME Clean one radar volume and save a compact copy.
%   Runs read_rd on RAW_FILE at full gate resolution, estimates the EP-VVP
%   velocity profile while the velocity sweeps are still loaded (VP bins:
%   height above the antenna), thins the gates of each density sweep to
%   FILTER.data_res = [range tangential] spacing (m), and saves only what
%   the forward model needs to OUT_FILE. STATUS is a short message; "ok" on
%   success.

try
    scan_cell = bmm4d.radar.read_rd(raw_file, read_options{:});
    scan = scan_cell{1};
    [~, ~, u, v] = epvvp(scan, vp.zstep_m, vp.rmin_m, vp.rmax_m, vp.zmax_m);
    sweeps = scan.dens.sweeps;
    for i_elev = 1:numel(sweeps)
        sweeps(i_elev).data = thin_sweep(sweeps(i_elev).data, sweeps(i_elev).range(:), sweeps(i_elev).az(:), filter.data_res);
    end
    sweeps = rmfield(sweeps, setdiff(fieldnames(sweeps), ["elev" "az" "range" "data"]));
    volume = struct( ...
        "station", string(scan.station), ...
        "lat", scan.lat, ...
        "lon", scan.lon, ...
        "height", scan.height, ...
        "datetime", scan.datetime, ...
        "uv", [u v], ...
        "dens", struct("sweeps", sweeps), ...
        "source_file", string(raw_file), ...
        "mistnet_used", isfield(scan, "mistnet_used") && scan.mistnet_used);
    save(out_file, "volume")
    status = "ok";
catch err
    status = string(err.message);
end
end

function data = thin_sweep(data, range_m, az_deg, data_res)
% Keep one gate every data_res = [range tangential] metres (same rule as
% read_rd's data_res); the others become missing.
res_range = max(1, floor(data_res(1) / median(diff(range_m))));
res_az = max(1, floor(data_res(2) ./ (deg2rad(median(diff(az_deg))) .* range_m)));
thinned = nan(size(data));
for i = 1:res_range:size(data, 1)
    thinned(i, 1:res_az(i):end) = data(i, 1:res_az(i):end);
end
data = thinned;
end
