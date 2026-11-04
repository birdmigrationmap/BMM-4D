function scan = estimate_density(scan, parms)
%ESTIMATE_DENSITY Bird density around every gate from all gates of its window.
%   Replaces the gate-by-gate cleaning (max_density, outlier_mad,
%   outlier_continuity, no_echo, to_density, smooth_median). Near the
%   radar a gate is so small that most gates hold no bird and the few
%   that hold one look very dense; filters that drop the high values or
%   take a median read the empty gates and lose most of the birds. Here
%   each gate is an interval of its summed bird echo (BMM4D.RADAR.GATE.BOUNDS)
%   and the density is the one that makes all gates within
%   +-smoothing/2 (m) most probable together, given how many birds a gate
%   of that size holds (BMM4D.RADAR.GATE.ESTIMATE_SWEEP;
%   studies/gate_density/report.html).
%   Gates without echo within noEchoRange are below the detection limit,
%   taken per range as the median over 9 range bins of the 1st percentile
%   of the sweep's values. Missing gates (rain, blocked beam, no echo
%   beyond noEchoRange) carry no information and stay missing; the
%   estimate around them uses the other gates of the window.
%   Works on reflectivity (dBZ) and returns density (birds km^-3), and in
%   each gate's SD field the standard deviation of that density from the
%   window likelihood (estimate_sweep).


% The bird-sum table is built once per MATLAB process (client or worker)
persistent table_s
if isempty(table_s)
    table_s = bmm4d.radar.gate.bird_sum_table(bmm4d.radar.gate.beam(), 2000, 20000);
end
s_floor_cm2 = 0.1;                                   % 1% of one bird on the beam axis
radius_m = parms.smoothing / 2;
if ~isfield(scan.gates, "sd")
    [scan.gates.sd] = deal([]);
end
for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    range_m = gate.range(:);
    beam = bmm4d.radar.gate.beam(gate_length_m=median(diff(range_m)));
    % Only the rows that range_limits keeps, and the window around them
    [~, height_m] = slant2ground(range_m, gate.elev);
    rows = range_m >= parms.minRange - radius_m & range_m <= parms.maxRange + radius_m & height_m <= parms.maxAlt + radius_m;
    gate.data(~rows, :) = nan;
    gate.no_echo(~rows, :) = false;
    no_echo = gate.no_echo & range_m <= parms.noEchoRange;
    detection_dbz = fillmissing(movmedian(prctile(gate.data, 1, 2), 9, "omitnan"), "nearest");
    [lo, hi] = bmm4d.radar.gate.bounds(gate.data, no_echo, range_m, beam, detection_dbz, s_floor_cm2);
    estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, 360 / numel(gate.az), beam, table_s, radius_m=radius_m);
    density = estimate.density;
    density(isnan(gate.data) & ~no_echo) = nan;
    gate.data = density;
    gate.sd = estimate.sd;
    gate.sd(isnan(density)) = nan;
    gate.unit = "birds km^-3";
    scan.gates(i_elev) = gate;
end
end
