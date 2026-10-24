function scan = rain_dr(scan, parms)
%RAIN_DR Remove rain: gates whose depolarization ratio is below rainDrThreshold.
%   Depolarization ratio from ZDR (dB) and rhoHV for simultaneous H/V
%   transmission (Kilambi et al. 2018):
%   DR = 10 log10((1 + Zdr - 2 rho sqrt(Zdr)) / (1 + Zdr + 2 rho sqrt(Zdr))),
%   Zdr linear. Rain (rhoHV ~ 1, ZDR ~ 0) gives very low DR. Sweeps whose
%   rhoHV or ZDR grid does not match the reflectivity are left unchanged.

zdr_sweeps = unique_elev_sweeps(scan, 'dr');
for i_elev = 1:numel(scan.gates)
    dz_data = scan.gates(i_elev).data;
    mask = false(size(dz_data));
    [~, i_elev_rh] = min(abs(scan.gates(i_elev).elev - [scan.rh.sweeps.elev]));
    [~, i_elev_zdr] = min(abs(scan.gates(i_elev).elev - [zdr_sweeps.elev]));
    rh_data = sweep2mat(scan.rh.sweeps(i_elev_rh));
    zdr_data = sweep2mat(zdr_sweeps(i_elev_zdr));
    if size(rh_data, 2) == size(dz_data, 2) && size(zdr_data, 2) == size(dz_data, 2)
        n_bins = min([size(dz_data, 1) size(rh_data, 1) size(zdr_data, 1)]);
        zdr_linear = 10 .^ (zdr_data(1:n_bins, :) / 10);
        cross_term = 2 * rh_data(1:n_bins, :) .* sqrt(zdr_linear);
        depolarization = 10 * log10(max(1 + zdr_linear - cross_term, eps) ./ (1 + zdr_linear + cross_term));
        mask(1:n_bins, :) = depolarization < parms.rainDrThreshold;
    end
    scan.gates(i_elev) = bmm4d.radar.filter.remove_rain(scan.gates(i_elev), mask);
end
end
