function status = clean_volume(raw_file, read_options, vp, filter, out_file)
%CLEAN_VOLUME Clean one radar volume and save a compact copy.
%   Runs read_rd on RAW_FILE at full gate resolution, estimates the EP-VVP
%   velocity profile while the velocity sweeps are still loaded (VP bins:
%   height above the antenna), applies the FILTER settings to each density
%   sweep, thins the gates, and saves only what the forward model needs to
%   OUT_FILE. STATUS is a short message; "ok" on success.
%
%   FILTER fields:
%   - max_density: gates above this density (birds km^-3) are removed.
%   - min_patch_km2: connected echo patches smaller than this are removed
%     (point targets, rain fragments, speckle).
%   - data_res: [range tangential] spacing (m) of the kept gates.

try
    scan_cell = read_rd(raw_file, read_options{:});
    scan = scan_cell{1};
    [~, ~, u, v] = epvvp(scan, vp.zstep_m, vp.rmin_m, vp.rmax_m, vp.zmax_m);
    sweeps = scan.dens.sweeps;
    for i_elev = 1:numel(sweeps)
        sweeps(i_elev).data = filter_sweep(sweeps(i_elev).data, sweeps(i_elev).range(:), sweeps(i_elev).az(:), filter);
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

function data = filter_sweep(data, range_m, az_deg, filter)
% Remove implausible densities and small isolated echo patches, then thin
% the gates (same rule as read_rd's data_res).
data(data > filter.max_density) = nan;

if filter.min_patch_km2 > 0
    gate_area_km2 = (range_m * deg2rad(median(diff(az_deg))) * median(diff(range_m))) / 1e6;
    % Patches of echo; zero-density (no-echo) gates separate patches
    patches = bwconncomp(data > 0, 8);
    patch_area = cellfun(@(ids) sum(gate_area_km2(mod(ids - 1, numel(range_m)) + 1)), patches.PixelIdxList);
    small = patches.PixelIdxList(patch_area < filter.min_patch_km2);
    data(vertcat(small{:})) = nan;
end

res_range = max(1, floor(filter.data_res(1) / median(diff(range_m))));
res_az = max(1, floor(filter.data_res(2) ./ (deg2rad(median(diff(az_deg))) .* range_m)));
thinned = nan(size(data));
for i = 1:res_range:size(data, 1)
    thinned(i, 1:res_az(i):end) = data(i, 1:res_az(i):end);
end
data = thinned;
end
