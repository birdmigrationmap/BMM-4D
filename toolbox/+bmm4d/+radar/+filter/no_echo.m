function scan = no_echo(scan, parms)
%NO_ECHO Keep no-echo gates within noEchoRange (m) as zero reflectivity.
%   Gates without echo in the raw sweep are set to -Inf dBZ (zero density):
%   an empty unblocked beam is evidence of no birds, not missing data.
%   Gates removed since (rain, blockage) are no longer no-echo gates.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    gate.data(gate.no_echo & gate.range <= parms.noEchoRange) = -inf;
    scan.gates(i_elev) = gate;
end
end
