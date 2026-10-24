function mean_density = gd_window_mean(density, range_m, az_step_deg, radius_m)
%GD_WINDOW_MEAN Linear mean of the gates within a radius of each gate.
%   Same window as GD_ESTIMATE_SWEEP (GD_WINDOW_SUM). DENSITY is
%   [n_range n_az] (birds km^-3, no echo = 0, NaN = left out).
valid = ~isnan(density);
density(~valid) = 0;
mean_density = gd_window_sum(density, range_m, az_step_deg, radius_m) ./ ...
    gd_window_sum(double(valid), range_m, az_step_deg, radius_m);
end
