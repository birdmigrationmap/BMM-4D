function status = ar_clean_volume(raw_file, read_options, vp, filter, out_file)
%AR_CLEAN_VOLUME bmm4d.radar.clean_volume with a symmetric aspect correction (dev).
%   STATUS = AR_CLEAN_VOLUME(RAW_FILE, READ_OPTIONS, VP, FILTER, OUT_FILE)
%   cleans one volume as bmm4d.radar.clean_volume (READ_OPTIONS must have
%   "shape", false), then corrects the bird aspect (side-on birds reflect
%   more than head- or tail-on) with a template that is symmetric front /
%   back and side / side and oriented by the birds' track from the
%   volume's own velocity profile (mean over 200-1500 m above the antenna):
%       T(az) = 1 + a cos(2 (az - track - 90 deg)),
%   a = the second azimuthal harmonic of the mean density at 10-60 km (all
%   sweeps pooled) projected on that axis, kept within [0, 0.5]. Each sweep
%   (and its sd) is divided by T and rescaled to keep its mean, as the
%   production shape correction. A gradient across the disc (first
%   harmonic) or a band through the radar at another angle is left alone,
%   so the correction cannot make one-sided wedges.

try
    scan_cell = bmm4d.radar.read_rd(raw_file, read_options{:});
    scan = scan_cell{1};
    [~, z_vp, u, v] = epvvp(scan, vp.zstep_m, vp.rmin_m, vp.rmax_m, vp.zmax_m);
    sweeps = scan.dens.sweeps;

    % Track from the velocity profile, then the aspect amplitude
    layer = z_vp(:) >= 200 & z_vp(:) <= 1500;
    track = atan2d(mean(u(layer), "omitnan"), mean(v(layer), "omitnan"));
    amplitude = 0;
    if isfinite(track)
        side = track + 90;
        az_all = [];
        mean_all = [];
        for i_elev = 1:numel(sweeps)
            rows = sweeps(i_elev).range(:) >= 10000 & sweeps(i_elev).range(:) <= 60000;
            sweep_mean = mean(sweeps(i_elev).data(rows, :), 1, "omitnan");
            az_all = [az_all; sweeps(i_elev).az(:)]; %#ok<AGROW>
            mean_all = [mean_all; sweep_mean(:)]; %#ok<AGROW>
        end
        ok = isfinite(mean_all);
        if nnz(ok) > 30 && sum(mean_all(ok)) > 0
            amplitude = 2 * sum(mean_all(ok) .* cosd(2 * (az_all(ok) - side))) / sum(mean_all(ok));
            amplitude = min(max(amplitude, 0), 0.5);
        end
    end
    for i_elev = 1:numel(sweeps)
        if amplitude == 0
            break
        end
        template = 1 + amplitude * cosd(2 * (sweeps(i_elev).az(:)' - side));
        before = mean(sweeps(i_elev).data(:), "omitnan");
        corrected = sweeps(i_elev).data ./ template;
        scale = before / mean(corrected(:), "omitnan");
        if ~isfinite(scale)
            scale = 1;
        end
        sweeps(i_elev).data = corrected * scale;
        if isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
            sweeps(i_elev).sd = sweeps(i_elev).sd ./ template * scale;
        end
    end

    % Thin and save as bmm4d.radar.clean_volume
    for i_elev = 1:numel(sweeps)
        sweeps(i_elev).data = thin_sweep(sweeps(i_elev).data, sweeps(i_elev).range(:), sweeps(i_elev).az(:), filter.data_res);
        if isfield(filter, "estimator_sd") && filter.estimator_sd && isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
            sweeps(i_elev).sd = thin_sweep(sweeps(i_elev).sd, sweeps(i_elev).range(:), sweeps(i_elev).az(:), filter.data_res);
        end
    end
    keep_fields = ["elev" "az" "range" "data"];
    if isfield(filter, "estimator_sd") && filter.estimator_sd
        keep_fields = [keep_fields "sd"];
    end
    sweeps = rmfield(sweeps, setdiff(fieldnames(sweeps), keep_fields));
    volume = struct("station", string(scan.station), "lat", scan.lat, "lon", scan.lon, "height", scan.height, ...
        "datetime", scan.datetime, "uv", [u v], "dens", struct("sweeps", sweeps), "source_file", string(raw_file), ...
        "mistnet_used", isfield(scan, "mistnet_used") && scan.mistnet_used, ...
        "aspect", struct("track_deg", track, "amplitude", amplitude));
    save(out_file, "volume")
    status = "ok";
catch err
    status = string(err.message);
end
end

function data = thin_sweep(data, range_m, az_deg, data_res)
% As bmm4d.radar.clean_volume: one gate every data_res = [range tangential] m
res_range = max(1, floor(data_res(1) / median(diff(range_m))));
res_az = max(1, floor(data_res(2) ./ (deg2rad(median(diff(az_deg))) .* range_m)));
thinned = nan(size(data));
for i = 1:res_range:size(data, 1)
    thinned(i, 1:res_az(i):end) = data(i, 1:res_az(i):end);
end
data = thinned;
end
