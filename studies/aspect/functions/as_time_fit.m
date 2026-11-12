function fit = as_time_fit(C, rows, n_boot)
%AS_TIME_FIT Aspect harmonics from the rotation of the track over the night (dev).
%   FIT = AS_TIME_FIT(C, ROWS, N_BOOT) uses the single-radar cell harmonics
%   of as_cell_harmonics (ROWS: logical subset) in the geographic frame:
%       A1 = (c1 + i s1) exp(i track),  A2 = (-c2 + i s2) exp(2 i track).
%   An aspect effect turns with the track: A1 = b1 exp(i track) and
%   A2 = -b2 exp(2 i track), while a radar's persistent pattern (blockage,
%   terrain, a lasting gradient) stays in the geographic frame. A fixed
%   effect per radar, sweep and range band removes the persistent part, so
%   b1 and b2 come from the change of the track during the night only
%   (complex least squares; real part: aspect along the track axis,
%   imaginary part: axis turned). Standard errors: bootstrap over radars
%   (N_BOOT draws, default 200).

if nargin < 3
    n_boot = 200;
end
C = C(rows, :);
[~, ~, group] = unique([C.radar C.sweep C.band], "rows");
x1 = exp(1i * deg2rad(C.track));
x2 = -exp(2i * deg2rad(C.track));
a1 = (C.c1 + 1i * C.s1) .* x1;
a2 = (-C.c2 + 1i * C.s2) .* exp(2i * deg2rad(C.track));
n_g = accumarray(group, 1);
keep = n_g(group) >= 3;
estimate = @(sel) solve(a1(sel), a2(sel), x1(sel), x2(sel), group(sel));
[b1, b2, spread] = estimate(keep);
radar_list = unique(C.radar(keep));
boot = nan(n_boot, 2);
stream = RandStream("mt19937ar", "Seed", 1);
for k = 1:n_boot
    draw = radar_list(randi(stream, numel(radar_list), numel(radar_list), 1));
    sel = [];
    for r = draw'
        sel = [sel; find(keep & C.radar == r)]; %#ok<AGROW>
    end
    % Duplicated radars get their own groups
    [~, ~, g_boot] = unique([group(sel) cumsum([1; diff(sel) <= 0])], "rows");
    [bb1, bb2] = solve(a1(sel), a2(sel), x1(sel), x2(sel), g_boot);
    boot(k, :) = [real(bb1) real(bb2)];
end
fit = struct("b1", real(b1), "b1_turned", imag(b1), "b2", real(b2), "b2_turned", imag(b2), ...
    "b1_se", std(boot(:, 1)), "b2_se", std(boot(:, 2)), "n", nnz(keep), "n_radar", numel(radar_list), ...
    "track_spread_deg", spread);
end

function [b1, b2, spread] = solve(a1, a2, x1, x2, group)
% Complex least squares after removing the mean of each group
[~, ~, group] = unique(group);
n = accumarray(group, 1);
centre = @(values) values - index(complex(accumarray(group, real(values)), accumarray(group, imag(values))) ./ n, group);
x1c = centre(x1);
x2c = centre(x2);
b1 = sum(conj(x1c) .* centre(a1)) / sum(abs(x1c) .^ 2);
b2 = sum(conj(x2c) .* centre(a2)) / sum(abs(x2c) .^ 2);
% Within-group spread of the track (deg), what identifies the fit
spread = rad2deg(sqrt(mean(abs(x1c) .^ 2)));
end

function values = index(values, id)
values = values(id);
end
