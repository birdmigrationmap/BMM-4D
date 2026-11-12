function C = as_cell_harmonics(volumes, version)
%AS_CELL_HARMONICS Single-radar harmonics of each sweep and range band (dev).
%   C = AS_CELL_HARMONICS(VOLUMES, VERSION) computes for every volume, sweep
%   and range band of the sector sums of as_frame_bins (VERSION "raw" or a
%   function handle f(volume) returning corrected sums) with at least 30 of
%   the 36 sectors seen, the harmonics of the sector means relative to the
%   track at the cell's mean height (velocity profile of the volume):
%   c2 = -2 sum(s cos 2 psi) / sum(s) (> 0: denser side-on), s2, c1
%   (> 0: denser downstream, tail-on), s1; with the cell's height above the
%   antenna, range band, sweep, the track and the speed at that height. Real gradients
%   add noise to each cell; their population mean estimates the aspect.

centre = (5:10:355)';
parts = cell(numel(volumes), 1);
for i = 1:numel(volumes)
    vol = volumes(i);
    if isa(version, "function_handle")
        sums = version(vol);
    else
        sums = vol.sectors.(version);
    end
    counts = vol.sectors.n;
    [n_sweep, n_band] = deal(size(counts, 2), size(counts, 3));
    rows = nan(n_sweep * n_band, 12);
    k = 0;
    for i_sweep = 1:n_sweep
        for i_band = 1:n_band
            n_c = counts(:, i_sweep, i_band);
            ok = n_c > 0;
            if nnz(ok) < 30
                continue
            end
            s = sums(ok, i_sweep, i_band) ./ n_c(ok);
            if sum(s) <= 0
                continue
            end
            h = sum(vol.sectors.height(:, i_sweep, i_band)) / sum(n_c);
            layer = min(size(vol.uv, 1), max(1, floor(h / 100) + 1));
            track = atan2d(vol.uv(layer, 1), vol.uv(layer, 2));
            if ~isfinite(track)
                continue
            end
            psi = centre(ok) - track;
            total = sum(s);
            k = k + 1;
            rows(k, :) = [vol.frame vol.radar i_sweep i_band h hypot(vol.uv(layer, 1), vol.uv(layer, 2)) ...
                -2 * sum(s .* cosd(2 * psi)) / total, 2 * sum(s .* sind(2 * psi)) / total, ...
                2 * sum(s .* cosd(psi)) / total, 2 * sum(s .* sind(psi)) / total, mean(s), mod(track, 360)];
        end
    end
    parts{i} = rows(1:k, :);
end
C = array2table(vertcat(parts{:}), VariableNames=["frame" "radar" "sweep" "band" "height_m" "speed" ...
    "c2" "s2" "c1" "s1" "mean_density" "track"]);
end
