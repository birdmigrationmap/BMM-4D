function gate = remove_rain(gate, mask)
%REMOVE_RAIN Add MASK to a sweep's rain mask and remove those gates.
%   Shared by the rain filters: rain gates become missing (NaN) and are no
%   longer no-echo gates.

gate.rain_mask = gate.rain_mask | mask;
gate.data(gate.rain_mask) = nan;
gate.no_echo(gate.rain_mask) = false;
end
