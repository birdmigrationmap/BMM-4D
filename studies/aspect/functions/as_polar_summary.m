function P = as_polar_summary(volumes)
%AS_POLAR_SUMMARY Orientation and polarimetric harmonics of each volume (dev).
%   P = AS_POLAR_SUMMARY(VOLUMES) for the volumes of as_frame_polar (sector
%   sums at 10-60 km) returns per volume:
%   - track, speed: velocity profile averaged over 200-1500 m above the antenna;
%   - heading_rho: azimuth of the minimum of the first harmonic of the
%     density-weighted rhoHV (tail-on view: rhoHV is lowest tail-on and
%     highest head-on, Stepanian et al. 2016), rho_amp its amplitude;
%   - drift: heading_rho - track (deg, -180..180);
%   - zdr_c2, zdr_s2: second harmonic (side-on axis from the track) of the
%     sector log(Z_H / Z_V) = ZDR (natural log units, as the density
%     harmonics); gradient-free since Z_H and Z_V see the same birds;
%   - zh_c2, zv_c2: the same for log Z_H and log Z_V separately
%     (single-radar: include the real gradients).

n = numel(volumes);
centre = (5:10:355)';
P = table('Size', [n 13], 'VariableTypes', repmat("double", 1, 13), 'VariableNames', ...
    ["frame" "radar" "track" "speed" "heading_rho" "rho_amp" "drift" "zdr_c2" "zdr_s2" "zdr_c1" "zh_c2" "zv_c2" "zdr_mean"]);
for i = 1:n
    vol = volumes(i);
    s = vol.sectors;
    P.frame(i) = vol.frame;
    P.radar(i) = vol.radar;
    z = (0:size(vol.uv, 1) - 1)' * 100 + 50;
    uv = mean(vol.uv(z >= 200 & z <= 1500, :), 1, "omitnan");
    P.track(i) = mod(atan2d(uv(1), uv(2)), 360);
    P.speed(i) = hypot(uv(1), uv(2));
    ok = s.n > 0 & s.density > 0 & s.zh > 0 & s.zv > 0;
    if nnz(ok) < 30
        P{i, 5:end} = nan;
        continue
    end
    rho = s.density_rho(ok) ./ s.density(ok);
    az = centre(ok);
    first = sum((rho - mean(rho)) .* exp(1i * deg2rad(az))) * 2 / nnz(ok);
    % first = amp * exp(i * az_max): rhoHV maximum (head-on) at az_max
    P.heading_rho(i) = mod(rad2deg(angle(first)) + 180, 360);
    P.rho_amp(i) = abs(first);
    P.drift(i) = mod(P.heading_rho(i) - P.track(i) + 180, 360) - 180;
    psi = az - P.track(i);
    log_zdr = log(s.zh(ok) ./ s.zv(ok));
    P.zdr_mean(i) = mean(log_zdr);
    P.zdr_c2(i) = -2 * mean((log_zdr - mean(log_zdr)) .* cosd(2 * psi));
    P.zdr_s2(i) = 2 * mean((log_zdr - mean(log_zdr)) .* sind(2 * psi));
    P.zdr_c1(i) = 2 * mean((log_zdr - mean(log_zdr)) .* cosd(psi));
    log_zh = log(s.zh(ok) ./ s.n(ok));
    log_zv = log(s.zv(ok) ./ s.n(ok));
    P.zh_c2(i) = -2 * mean((log_zh - mean(log_zh)) .* cosd(2 * psi));
    P.zv_c2(i) = -2 * mean((log_zv - mean(log_zv)) .* cosd(2 * psi));
end
end
