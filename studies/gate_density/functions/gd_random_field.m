function density = gd_random_field(range_m, az_step_deg, mean_density, log_sd, length_m)
%GD_RANDOM_FIELD Smooth lognormal bird density on a sweep's polar gates.
%   White noise on a 500 m Cartesian grid is smoothed with a Gaussian of
%   standard deviation LENGTH_M, scaled to unit variance, and turned into
%   a lognormal field with mean MEAN_DENSITY and log standard deviation
%   LOG_SD. DENSITY is [n_range n_az] at the gate centres.

range_m = range_m(:);
az = 0:az_step_deg:360 - az_step_deg;
extent = max(range_m) + 3 * length_m;
step = 500;
axis_m = -extent:step:extent;
field = imgaussfilt(randn(numel(axis_m), numel(axis_m)), length_m / step);
field = field / std(field(:));
[x, y] = ndgrid(axis_m, axis_m);
gate_x = range_m .* sind(az);
gate_y = range_m .* cosd(az);
z = interpn(x, y, field, gate_x, gate_y, "linear");
density = mean_density * exp(log_sd * z - log_sd ^ 2 / 2);
end
