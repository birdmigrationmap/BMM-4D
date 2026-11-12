function scan = beam_blockage(scan, ~)
%BEAM_BLOCKAGE Remove gates whose beam is blocked by terrain (dev override).
%   Study copy of the toolbox filter (studies/artifacts), put first on the
%   path only by the artefact cleaning jobs (ar_clean_volumes). The
%   beam-blockage maps (studies/beam_blockage/beam_blockage_script.m) were
%   built with pol2cart(azimuth, range), i.e. with the angle counted
%   counter-clockwise from east, while radar azimuths are compass bearings
%   (clockwise from north). The map column for map angle a holds the
%   terrain at compass bearing 90 - a, so a gate at compass azimuth az is
%   looked up at map angle 90 - az. The toolbox filter looks up az itself:
%   it removed gates in the mirrored direction and kept the blocked ones,
%   whose missing echo then became zero-density observations.

persistent bb
if isempty(bb)
    bb = load(bmm4d.util.data_path("reference", "beam_blockage", "beam_blockage.mat"));
end
bbr_f = griddedInterpolant({bb.bb_range, bb.bb_az}, bb.bb(:,:,scan.station==bb.bb_radar),'nearest','nearest');
for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    blocked = bbr_f({gate.range, mod(90 - gate.az, 360)}) > (gate.elev - gate.vert_half_bw);
    gate.data(blocked) = nan;
    gate.no_echo(blocked) = false;
    scan.gates(i_elev) = gate;
end
end
