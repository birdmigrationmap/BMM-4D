function scan = rain_rho(scan, parms)
%RAIN_RHO Remove rain: gates whose rhoHV exceeds rainRhoThreshold.
%   Rain is round and uniform, so its correlation coefficient is close to
%   1. Beyond the rhoHV sweep's last bin every gate counts as rain. The
%   mask is added to gates(i).rain_mask (for rain_dilation) and the gates
%   are removed, including from the no-echo gates.

for i_elev = 1:numel(scan.gates)
    [min_elev, i_elev_rh] = min(abs(scan.gates(i_elev).elev - vertcat(scan.rh.sweeps.elev)));
    assert(min_elev < .01)
    rh_data = sweep2mat(scan.rh.sweeps(i_elev_rh));
    mask = true(size(scan.gates(i_elev).data));
    min_nbins = min(size(rh_data, 1), size(mask, 1));
    mask(1:min_nbins, :) = rh_data(1:min_nbins, :) > parms.rainRhoThreshold;
    scan.gates(i_elev) = bmm4d.radar.filter.remove_rain(scan.gates(i_elev), mask);
end
end
