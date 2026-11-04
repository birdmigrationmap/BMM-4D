function [log10_lo, log10_hi] = bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2)
%BOUNDS What each gate tells about its summed bird contribution.
%   Converts a sweep DBZ [n_range n_az] into intervals [LOG10_LO, LOG10_HI]
%   of s = eta * V_eff (log10 cm^2, BMM4D.RADAR.GATE.BEAM):
%   - a value: the value +-0.25 dB (quantization);
%   - no echo (NO_ECHO true): below the detection limit DETECTION_DBZ (one
%     per range);
%   - other missing gates (NaN, e.g. rain or blocked beam): no
%     information, [-Inf Inf];
%   - anything below S_FLOOR_CM2 (a small fraction of one bird) cannot be
%     told from no bird or a background echo (insects, noise): lower end
%     -Inf.

range_m = range_m(:);
volume = beam.volume_eff_km3(range_m);
log10_s = log10(z_to_refl(idb(dbz)) .* volume);
log10_det = log10(z_to_refl(idb(detection_dbz(:))) .* volume) .* ones(1, size(dbz, 2));
floor_value = log10(s_floor_cm2);
log10_lo = log10_s - 0.025;
log10_hi = log10_s + 0.025;
log10_lo(no_echo) = -Inf;
log10_hi(no_echo) = log10_det(no_echo);
log10_lo(log10_lo < floor_value) = -Inf;
log10_hi(log10_hi < floor_value) = floor_value;
unknown = isnan(dbz) & ~no_echo;
log10_lo(unknown) = -Inf;
log10_hi(unknown) = Inf;
end
