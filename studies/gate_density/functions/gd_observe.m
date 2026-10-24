function dbz = gd_observe(eta, detection_dbz, background_dbz)
%GD_OBSERVE What the radar outputs for simulated bird reflectivity.
%   Adds a weak background echo (insects, noise: lognormal around
%   BACKGROUND_DBZ with a 5 dB spread; -Inf for none) to ETA (cm^2 km^-3,
%   [n_range n_az]), converts to dBZ and sets gates below DETECTION_DBZ
%   (one value per range) to NaN (no echo).

background = z_to_refl(idb(background_dbz + 5 * randn(size(eta))));
dbz = db(refl_to_z(eta + background));
dbz(dbz < detection_dbz(:)) = nan;
end
