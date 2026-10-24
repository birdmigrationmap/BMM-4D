function scan = range_limits(scan, parms)
%RANGE_LIMITS Remove gates closer than minRange, beyond maxRange or above maxAlt (m).
%   The altitude is the beam-centre height above the antenna
%   (slant2ground).

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    gate.data(gate.range < parms.minRange, :) = nan;
    gate.data(gate.range > parms.maxRange, :) = nan;
    [~, sweep_height] = slant2ground(gate.range, gate.elev);
    gate.data(sweep_height > parms.maxAlt, :) = nan;
    scan.gates(i_elev) = gate;
end
end
