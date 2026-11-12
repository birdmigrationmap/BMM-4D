function score = ar_terrain_score(S, column)
%AR_TERRAIN_SCORE How much the column follows the terrain (dev).
%   SCORE = AR_TERRAIN_SCORE(S, COLUMN) for column maps [nx ny n_time]
%   (ar_columns), over the frames:
%   - ridge_valley: mean column where the ground is above 1000 m against
%     where it is below 450 m (1 = no deficit over the ridges);
%   - per_km: column ratio per km of ground, from the regression of the
%     small-scale log column on the small-scale ground (both minus their
%     50 km Gaussian mean): 1 = the column ignores the terrain.

dem = max(S.g.dem, 0);
active_column = any(S.active_mask_3d, 3);
mean_column = mean(column, 3, "omitnan");
ridge = active_column & dem > 1000;
valley = active_column & dem < 450;
ridge_valley = mean(mean_column(ridge)) / mean(mean_column(valley));

% Small scales: minus the 50 km Gaussian mean (normalized over active columns)
sigma_cells = 50000 / S.g.dx;
smooth = @(m) imgaussfilt(m .* active_column, sigma_cells) ./ max(imgaussfilt(double(active_column), sigma_cells), eps);
log_column = log(max(mean_column, 1));
log_column(~active_column) = 0;
small_column = log_column - smooth(log_column);
small_dem = dem / 1000 - smooth(dem / 1000);
use = active_column & imgaussfilt(double(active_column), sigma_cells) > 0.9;
slope = small_dem(use) \ small_column(use);
score = struct("ridge_valley", ridge_valley, "per_km", exp(slope));
end
