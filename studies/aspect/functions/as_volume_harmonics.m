function H = as_volume_harmonics(volumes, version, band_ids, layer_m)
%AS_VOLUME_HARMONICS Azimuthal harmonics of each volume relative to the track (dev).
%   H = AS_VOLUME_HARMONICS(VOLUMES, VERSION, BAND_IDS, LAYER_M) uses the
%   sector sums of as_frame_bins (VERSION "raw", "ell" or "fix"; or a
%   function handle f(volume) returning a [36 x 12 x 5] sector-sum array),
%   pools the sweeps and the range bands BAND_IDS (1: 5-10, 2: 10-30,
%   3: 30-60, 4: 60-100, 5: 100-150 km) into the mean density of 36 sectors
%   of 10 deg, and projects it on harmonics of psi = az - track, with the
%   track from the volume's velocity profile averaged over LAYER_M (m above
%   the antenna, default [200 1500]):
%       c1 = 2 sum(s cos psi) / sum(s)   (> 0: denser downstream, tail-on side)
%       s1 = 2 sum(s sin psi) / sum(s)
%       c2 = -2 sum(s cos 2 psi) / sum(s) (> 0: denser side-on)
%       s2 = 2 sum(s sin 2 psi) / sum(s)
%   first = magnitude of the first harmonic (gradient across the disc).
%   Also the wedge score (deepest sector against the lower of the sectors
%   30 deg to each side; as ar_wedge_score but on the polar data at the
%   pooled bands) and the track and speed. Volumes with fewer than 30
%   sectors with data get NaN.

if nargin < 4
    layer_m = [200 1500];
end
% LAYER_M = "local": each sweep and band uses the track at its own height
local_track = isstring(layer_m) && layer_m == "local";
if local_track
    layer_m = [200 1500];
end
n = numel(volumes);
centre = (5:10:355)';
H = table('Size', [n 12], 'VariableTypes', repmat("double", 1, 12), ...
    'VariableNames', ["frame" "radar" "track" "speed" "c1" "s1" "c2" "s2" "first" "wedge" "wedge_az" "mean_density"]);
for i = 1:n
    vol = volumes(i);
    H.frame(i) = vol.frame;
    H.radar(i) = vol.radar;
    z = (0:size(vol.uv, 1) - 1)' * 100 + 50;
    layer = z >= layer_m(1) & z <= layer_m(2);
    uv = mean(vol.uv(layer, :), 1, "omitnan");
    H.track(i) = mod(atan2d(uv(1), uv(2)), 360);
    H.speed(i) = hypot(uv(1), uv(2));
    if isa(version, "function_handle")
        sums = version(vol);
    else
        sums = vol.sectors.(version);
    end
    counts = vol.sectors.n;
    s = sum(sums(:, :, band_ids), [2 3]) ./ sum(counts(:, :, band_ids), [2 3]);
    ok = isfinite(s);
    H.mean_density(i) = mean(s(ok));
    if nnz(ok) < 30 || ~isfinite(H.track(i)) || sum(s(ok)) <= 0
        H{i, 5:12} = nan;
        continue
    end
    if local_track
        % Each sweep and band with the track at its own mean height, then
        % the harmonics averaged with weights sum(s) of each sweep and band
        [c, total] = deal(zeros(1, 4), 0);
        for i_sweep = 1:size(sums, 2)
            for i_band = band_ids
                n_c = counts(:, i_sweep, i_band);
                s_c = sums(:, i_sweep, i_band) ./ n_c;
                ok_c = n_c > 0;
                if nnz(ok_c) < 30 || sum(s_c(ok_c)) <= 0
                    continue
                end
                h = sum(vol.sectors.height(:, i_sweep, i_band)) / sum(n_c);
                i_layer = min(size(vol.uv, 1), max(1, floor(h / 100) + 1));
                track_c = atan2d(vol.uv(i_layer, 1), vol.uv(i_layer, 2));
                if ~isfinite(track_c)
                    continue
                end
                psi = centre(ok_c) - track_c;
                w = s_c(ok_c);
                c = c + 2 * [sum(w .* cosd(psi)) sum(w .* sind(psi)) -sum(w .* cosd(2 * psi)) sum(w .* sind(2 * psi))];
                total = total + sum(w);
            end
        end
        c = c / total;
        [H.c1(i), H.s1(i), H.c2(i), H.s2(i)] = deal(c(1), c(2), c(3), c(4));
    else
        psi = centre(ok) - H.track(i);
        total = sum(s(ok));
        H.c1(i) = 2 * sum(s(ok) .* cosd(psi)) / total;
        H.s1(i) = 2 * sum(s(ok) .* sind(psi)) / total;
        H.c2(i) = -2 * sum(s(ok) .* cosd(2 * psi)) / total;
        H.s2(i) = 2 * sum(s(ok) .* sind(2 * psi)) / total;
    end
    H.first(i) = hypot(H.c1(i), H.s1(i));
    ratio = s ./ min(circshift(s, 3), circshift(s, -3));
    [H.wedge(i), k] = min(ratio);
    H.wedge_az(i) = centre(k);
end
end
