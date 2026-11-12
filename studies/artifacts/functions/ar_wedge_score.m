function score = ar_wedge_score(S, column)
%AR_WEDGE_SCORE Azimuthal wedges around the radars in the column maps (dev).
%   SCORE = AR_WEDGE_SCORE(S, COLUMN) for column maps [nx ny n_time]
%   (ar_columns): per radar and frame, the mean column in 10 deg sectors at
%   10-70 km, and the deepest sector against the lower of the two sectors
%   30 deg to each side (1 = no wedge; 0.7 = a sector 30% below both
%   sides). Returns per_radar [n_radar n_time] and the summary: median
%   over radars and frames, share below 0.8, and the deepest-sector
%   azimuth per_radar_azimuth (deg).

g = S.g;
[X, Y] = ndgrid(g.x, g.y);
radar_x = g.f_X(S.radar_selected.lon);
radar_y = g.f_Y(S.radar_selected.lat);
n_radar = numel(radar_x);
n_time = size(column, 3);
edges = 0:10:360;
per_radar = nan(n_radar, n_time);
per_radar_azimuth = nan(n_radar, n_time);
for i_r = 1:n_radar
    dx = X - radar_x(i_r);
    dy = Y - radar_y(i_r);
    r = hypot(dx, dy);
    near = r >= 10000 & r <= 70000;
    sector = discretize(mod(atan2d(dx(near), dy(near)), 360), edges);
    for i_t = 1:n_time
        values = column(:, :, i_t);
        values = values(near);
        valid = isfinite(values);
        s = accumarray(sector(valid), values(valid), [numel(edges) - 1 1], @mean, nan);
        if nnz(isfinite(s)) < 30
            continue
        end
        ratio = s ./ min(circshift(s, 3), circshift(s, -3));
        [per_radar(i_r, i_t), k] = min(ratio);
        per_radar_azimuth(i_r, i_t) = edges(k) + 5;
    end
end
score = struct("per_radar", per_radar, "per_radar_azimuth", per_radar_azimuth, ...
    "median", median(per_radar, "all", "omitnan"), "share_below_08", mean(per_radar(isfinite(per_radar)) < 0.8));
end
