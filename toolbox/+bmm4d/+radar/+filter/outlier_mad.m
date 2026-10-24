function scan = outlier_mad(scan, parms)
%OUTLIER_MAD Remove gates far from the sweep's median reflectivity.
%   Gates outside median +- removeOutliers x MAD (normal-consistent
%   median absolute deviation) of the whole sweep, in dBZ, are removed.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    mdata = median(gate.data(:),'omitnan');
    MAD = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(gate.data(:) - mdata), 'omitnan');
    thr = mdata + [-1 1]* parms.removeOutliers*MAD;
    id = gate.data < thr(1) | gate.data > thr(2);
    if ismember("outliers_mad", parms.plot_opts) && i_elev == 1
        fig = bmm4d.plot.create_export_figure(1400, 500);
        tiledlayout(1,3,'TileSpacing','tight','Padding','tight');
        nexttile; histogram(gate.data(:));
        title("Histogram of all data")
        xline(mdata,'-k',"Median", 'LineWidth',2); xline(thr,'--k',"Outlier threshold", 'LineWidth',2)
        xlabel("reflectivity");
        nexttile;
        z = bmm4d.plot.plot_sweep(gate.data, gate.az, gate.range, max_range=parms.maxRange, show_colorbar=false);
        c_lim = [0 quantile(z(:),.999)];
        clim(c_lim); title("Before threshold of MAD")
        nexttile;
        tmp = gate.data; tmp(id) = nan;
        bmm4d.plot.plot_sweep(tmp, gate.az, gate.range, max_range=parms.maxRange, ...
            c_lim=c_lim, ttl="After threshold of MAD between " + num2str(round(thr(1))) + " and " + num2str(round(thr(2))) + " bird/km^3", show_colorbar=false);
        bmm4d.radar.filter.finish_plot(fig, parms, "_outliers_mad.png")
    end
    gate.data(id) = nan;
    scan.gates(i_elev) = gate;
end
end
