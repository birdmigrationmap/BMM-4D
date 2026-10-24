function scan = max_density(scan, parms)
%MAX_DENSITY Remove gates whose bird density exceeds thr_dens (birds km^-3).
%   Works on reflectivity (dBZ); 3000 birds km^-3 is about 40 dBZ.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    sweep_density = z_to_refl(idb(gate.data)) / 11;   % bird density [bird/km^3]
    if ismember("outliers_threshold", parms.plot_opts) && i_elev == 1
        fig = bmm4d.plot.create_export_figure(1200, 550);
        tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
        nexttile;
        z = bmm4d.plot.plot_sweep(gate.data, gate.az, gate.range, max_range=parms.maxRange, show_colorbar=false);
        c_lim = [0 quantile(z(:),.999)];
        clim(c_lim); title("Before threshold of high value")
        nexttile;
        tmp = gate.data; tmp(sweep_density > parms.thr_dens) = nan;
        bmm4d.plot.plot_sweep(tmp, gate.az, gate.range, max_range=parms.maxRange, ...
            c_lim=c_lim, ttl="After threshold of high value " + num2str(parms.thr_dens,1) + " bird/km^3", show_colorbar=false);
        bmm4d.radar.filter.finish_plot(fig, parms, "_outliers_threshold.png")
    end
    gate.data(sweep_density > parms.thr_dens) = nan;
    scan.gates(i_elev) = gate;
end
end
