function scan = rain_dilation(scan, parms)
%RAIN_DILATION Grow the rain mask by rainDilation (m) to remove rain edges.
%   Weak rain edges look like biology. The mask of the rain filters run
%   before is grown by rainDilation along range, and in azimuth by the
%   same distance at 50 km (where most bird echo lies), wrapping around
%   north; the grown gates are removed.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    n_range = ceil(parms.rainDilation / median(diff(gate.range)));
    n_az = ceil(parms.rainDilation / (deg2rad(median(diff(gate.az))) * 50000));
    n_rows = min(size(gate.rain_mask, 1), numel(gate.range));
    % Pad in azimuth so the dilation wraps around north
    padded = gate.rain_mask(1:n_rows, [end-n_az+1:end 1:end 1:n_az]);
    padded = imdilate(padded, true(2*n_range+1, 2*n_az+1));
    mask = gate.rain_mask;
    mask(1:n_rows, :) = padded(:, n_az+1:end-n_az);
    scan.gates(i_elev) = bmm4d.radar.filter.remove_rain(gate, mask);
end
end
