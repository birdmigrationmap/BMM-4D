function [log10_lo, log10_hi] = gd_gate_bounds(dbz, cfp_code, range_m, beam, detection_dbz, s_floor_cm2)
%GD_GATE_BOUNDS What each gate tells about its summed bird contribution.
%   Converts a sweep (DBZ [n_range n_az], NaN = no echo; CFP_CODE the
%   clutter filter power removed codes, or zeros) into intervals
%   [LOG10_LO, LOG10_HI] of s = eta * V_eff (log10 cm^2, GD_BEAM):
%   - echo, filter not applied (code 0): the value, +-0.25 dB (quantization);
%   - echo, power removed (code >= 8, CFP = code - 8 dB): between the output
%     and the output plus CFP (the removed power is clutter and/or birds);
%   - echo, point clutter filter (code 1): at least the output (the
%     radar replaced a point target, the value removed is unknown);
%   - no echo: below the detection limit DETECTION_DBZ (one per range),
%     plus CFP if power was removed;
%   - anything below S_FLOOR_CM2 (a small fraction of one bird) cannot be
%     told from no bird or a background echo (insects, noise): lower end
%     -Inf.

range_m = range_m(:);
volume = beam.volume_eff_km3(range_m);
log10_s = log10(z_to_refl(idb(dbz)) .* volume);
log10_det = log10(z_to_refl(idb(detection_dbz(:))) .* volume) .* ones(1, size(dbz, 2));
removed_db = max(double(cfp_code) - 8, 0);
log10_lo = log10_s - 0.025;
log10_hi = log10_s + 0.025 + removed_db / 10;
no_echo = isnan(dbz);
log10_lo(no_echo) = -Inf;
log10_hi(no_echo) = log10_det(no_echo) + removed_db(no_echo) / 10;
log10_hi(cfp_code == 1 & ~no_echo) = Inf;
floor_value = log10(s_floor_cm2);
log10_lo(log10_lo < floor_value) = -Inf;
log10_hi(log10_hi < floor_value) = floor_value;
end
