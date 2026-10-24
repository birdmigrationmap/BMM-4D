function scan = beam_blockage(scan, ~)
%BEAM_BLOCKAGE Remove gates whose beam is blocked by terrain.
%   A gate is blocked when the terrain horizon (precomputed beam-blockage
%   map, data/reference/beam_blockage) is above the lower half-power edge
%   of the beam (elevation minus half the vertical beam width). Blocked
%   gates are missing, also when they had no echo.

% The map is loaded once per MATLAB process (client or worker).
persistent bb
if isempty(bb)
    bb = load(bmm4d.util.data_path("reference", "beam_blockage", "beam_blockage.mat"));
end
bbr_f = griddedInterpolant({bb.bb_range, bb.bb_az}, bb.bb(:,:,scan.station==bb.bb_radar),'nearest','nearest');
for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    blocked = bbr_f({gate.range, gate.az}) > (gate.elev - gate.vert_half_bw);
    gate.data(blocked) = nan;
    gate.no_echo(blocked) = false;
    scan.gates(i_elev) = gate;
end
end
