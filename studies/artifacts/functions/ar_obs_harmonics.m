function h = ar_obs_harmonics(S)
%AR_OBS_HARMONICS Azimuthal harmonics of each radar's observations (dev).
%   H = AR_OBS_HARMONICS(S) for a case (ar_load): per radar and frame, the
%   observed rows at 10-60 km are averaged in 10 deg sectors of azimuth
%   (row centroid seen from the radar), and with the track (mean profile
%   velocity 200-1500 m above the radar):
%   - first [n_radar n_time]: amplitude of the first harmonic (a gradient
%     across the disc, plus any front/back aspect effect);
%   - side_on: second harmonic on the side-on axis (track + 90 deg), the
%     side-on aspect effect;
%   - front_back: first harmonic on the track axis, cos(az - track): > 0
%     when the birds are denser where they fly away from the radar.

[X, Y] = ndgrid(S.g.x, S.g.y);
xy = [repmat(X(:), S.g.sz(3), 1) repmat(Y(:), S.g.sz(3), 1)];
centroid = (S.observations.F * xy) ./ full(sum(S.observations.F, 2));
radar_id = double(S.observations.row_radar_id);
azimuth = mod(atan2d(centroid(:, 1) - S.g.f_X(S.radar_selected.lon(radar_id)), ...
    centroid(:, 2) - S.g.f_Y(S.radar_selected.lat(radar_id))), 360);
near = S.observations.row_range_m >= 10000 & S.observations.row_range_m <= 60000;
[row, i_t, value] = find(S.observations.values);
keep = near(row);
n_radar = height(S.radar_selected);
n_time = size(S.observations.values, 2);
sector_mean = accumarray([discretize(azimuth(row(keep)), 0:10:360) radar_id(row(keep)) i_t(keep)], ...
    double(value(keep)), [36 n_radar n_time], @mean, nan);
z = squeeze(S.g.z);
centre = (5:10:355)';
h = struct("first", nan(n_radar, n_time), "side_on", nan(n_radar, n_time), "front_back", nan(n_radar, n_time));
for i_r = 1:n_radar
    layer = z > S.radar_selected.height(i_r) + 200 & z < S.radar_selected.height(i_r) + 1500;
    for k = 1:n_time
        s = sector_mean(:, i_r, k);
        ok = isfinite(s);
        if nnz(ok) < 30
            continue
        end
        profile = bmm4d.prior.radar_uv_on_grid(S.radar_scans{i_r, k}, S.g);
        velocity = mean(profile(layer, :), 1, "omitnan");
        track = atan2d(velocity(1), velocity(2));
        h.first(i_r, k) = 2 * abs(sum(s(ok) .* exp(-1i * deg2rad(centre(ok))))) / sum(s(ok));
        h.side_on(i_r, k) = 2 * sum(s(ok) .* cosd(2 * (centre(ok) - track - 90))) / sum(s(ok));
        h.front_back(i_r, k) = 2 * sum(s(ok) .* cosd(centre(ok) - track)) / sum(s(ok));
    end
end
end
