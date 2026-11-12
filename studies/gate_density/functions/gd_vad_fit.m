function [coef, fit_rms] = gd_vad_fit(vel, az, use)
%GD_VAD_FIT Velocity azimuth display fit of one ring of gates.
%   VEL is [n_gate n_az] (m/s), AZ [1 n_az] (deg), USE the gates to fit
%   (unfiltered, with a velocity). Fits v_r = c1 + c2 sin(az) + c3 cos(az),
%   then refits without residuals beyond 3 rms. (c2, c3) is the horizontal
%   motion (east, north) times the cosine of the elevation. No dealiasing:
%   rings with folded velocities show up as a large FIT_RMS.

az_grid = repmat(az(:)', size(vel, 1), 1);
design = [ones(nnz(use), 1) sind(az_grid(use)) cosd(az_grid(use))];
v = vel(use);
coef = design \ v;
inlier = abs(v - design * coef) < 3 * rms(v - design * coef);
coef = design(inlier, :) \ v(inlier);
fit_rms = rms(v(inlier) - design(inlier, :) * coef);
end
