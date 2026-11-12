function status = as_clean_volume(raw_file, read_options, vp, filter, out_file)
%AS_CLEAN_VOLUME bmm4d.radar.clean_volume with a shared bird aspect correction (dev).
%   STATUS = AS_CLEAN_VOLUME(RAW_FILE, READ_OPTIONS, VP, FILTER, OUT_FILE)
%   cleans one volume as bmm4d.radar.clean_volume (READ_OPTIONS must have
%   "shape", false), then divides every gate by the bird aspect factor of
%   as_aspect_template (as_apply_aspect), before thinning:
%       log T = b2(h) (-cos 2 psi) + b1(h) cos psi - log N,
%       psi = az - track(h),  b2(h) = b2_per_km h,
%       b1(h) = b1 min(1, h / b1_ramp_m),
%   h the beam-centre height above the antenna, track(h) the direction of
%   the volume's own velocity profile at that height (EP-VVP, 100 m
%   layers), N the normalisation to an azimuthal mean of one (a uniform
%   field keeps its mean; no rescaling). Side-on birds (psi = +-90 deg)
%   are brighter, tail-on birds (psi = 0: flying away from the radar)
%   dimmer than head-on ones.
%   The parameters are shared (not fitted to the volume, so real gradients
%   and bands are left alone): FILTER.aspect (struct with b2, b2_per_km,
%   b2_per_ms, h0_m, v0_ms, b1, b1_ramp_m; see as_aspect_template), by
%   default the fit to the radar overlap of the us_east night of 21-22 Oct
%   2026 (studies/aspect, AS_B_overlap.m): b2_per_km = 0.092, b1 = -0.21,
%   b1_ramp_m = 1000. FILTER.aspect = [] saves the volume uncorrected.
%   The parameters used are saved in volume.aspect.

try
    aspect = struct("b2", 0, "b2_per_km", 0.092, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, ...
        "b1", -0.21, "b1_ramp_m", 1000);
    if isfield(filter, "aspect")
        aspect = filter.aspect;
    end
    scan_cell = bmm4d.radar.read_rd(raw_file, read_options{:});
    scan = scan_cell{1};
    [~, ~, u, v] = epvvp(scan, vp.zstep_m, vp.rmin_m, vp.rmax_m, vp.zmax_m);
    sweeps = scan.dens.sweeps;
    keep_fields = ["elev" "az" "range" "data"];
    if isfield(filter, "estimator_sd") && filter.estimator_sd
        keep_fields = [keep_fields "sd"];
    end
    sweeps = rmfield(sweeps, setdiff(fieldnames(sweeps), keep_fields));
    volume = struct("station", string(scan.station), "lat", scan.lat, "lon", scan.lon, "height", scan.height, ...
        "datetime", scan.datetime, "uv", [u v], "dens", struct("sweeps", sweeps), "source_file", string(raw_file), ...
        "mistnet_used", isfield(scan, "mistnet_used") && scan.mistnet_used, "aspect", []);

    % Aspect correction on the full-resolution gates
    if ~isempty(aspect)
        volume = as_apply_aspect(volume, aspect);
    end

    % Thin and save as bmm4d.radar.clean_volume
    for i_elev = 1:numel(volume.dens.sweeps)
        sweep = volume.dens.sweeps(i_elev);
        volume.dens.sweeps(i_elev).data = thin_sweep(sweep.data, sweep.range(:), sweep.az(:), filter.data_res);
        if ismember("sd", keep_fields) && ~isempty(sweep.sd)
            volume.dens.sweeps(i_elev).sd = thin_sweep(sweep.sd, sweep.range(:), sweep.az(:), filter.data_res);
        end
    end
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
