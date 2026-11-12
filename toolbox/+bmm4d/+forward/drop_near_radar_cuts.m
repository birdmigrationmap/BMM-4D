function [scans, n_dropped] = drop_near_radar_cuts(scans, waveforms, max_range_m)
%DROP_NEAR_RADAR_CUTS Remove the gates of chosen cuts close to the radar.
%   [SCANS, N_DROPPED] = DROP_NEAR_RADAR_CUTS(SCANS, WAVEFORMS, MAX_RANGE_M)
%   sets to missing (NaN) the gates closer than MAX_RANGE_M (m, slant
%   range) of every sweep whose waveform (clean_volume) is one of
%   WAVEFORMS, in each cleaned volume of SCANS. The steeper batch cuts read
%   less than the split cuts within ~20 km of the radar (studies/near_radar);
%   removing them lets the split cuts, the neighbouring radars and the
%   advection fill the column there. The split cuts (0.5-1.3 deg, which
%   carry the lowest layer) are never removed. N_DROPPED counts the
%   observed gates (data not NaN) removed. An empty WAVEFORMS or a zero
%   MAX_RANGE_M removes nothing.

arguments
    scans cell
    waveforms string
    max_range_m (1,1) double {mustBeNonnegative}
end
assert(~any(waveforms == "split"), "bmm4d:drop_near_radar_cuts:Split", ...
    "The split cuts carry the lowest layer and are never dropped.")
n_dropped = 0;
if isempty(waveforms) || max_range_m == 0
    return
end
for k = 1:numel(scans)
    sweeps = scans{k}.dens.sweeps;
    if ~isfield(sweeps, "waveform")
        error("bmm4d:drop_near_radar_cuts:NoWaveform", ...
            "%s has no sweep waveform: clean it again (cleaning_version 3).", scans{k}.station)
    end
    for i_elev = find(ismember([sweeps.waveform], waveforms))
        near = sweeps(i_elev).range(:) < max_range_m;
        n_dropped = n_dropped + nnz(~isnan(sweeps(i_elev).data(near, :)));
        sweeps(i_elev).data(near, :) = nan;
        if isfield(sweeps, "sd") && ~isempty(sweeps(i_elev).sd)
            sweeps(i_elev).sd(near, :) = nan;
        end
    end
    scans{k}.dens.sweeps = sweeps;
end
end
