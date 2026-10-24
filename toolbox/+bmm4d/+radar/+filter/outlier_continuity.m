function scan = outlier_continuity(scan, parms)
%OUTLIER_CONTINUITY Remove gates that break spatial continuity.
%   Gates outside the 5 x 5 (range x azimuth) moving median +- 3 x the
%   sweep's mean local MAD, in dBZ, are removed. The last two range bins
%   are removed too (the neighbourhood wraps around there).

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    sweep_dbz = gate.data;
    sweep_dbz(end-1:end,:) = nan;
    dd = repmat(sweep_dbz,1,1,5*5); u=1;
    for i=-2:2
        for j=-2:2
            dd(:,:,u) = circshift(sweep_dbz,[i j]);
            u=u+1;
        end
    end
    % Medians only where the 5x5 neighbourhood has data (elsewhere they
    % are NaN anyway); same result, much faster on mostly empty sweeps
    dd = reshape(dd, [], 5*5);
    has_data = any(~isnan(dd), 2);
    data_median = nan(size(sweep_dbz));
    data_median(has_data) = median(dd(has_data,:), 2, 'omitnan');
    data_median(end-1:end,:)=nan;

    MAD = nan(size(sweep_dbz));
    MAD(has_data) = median(abs(dd(has_data,:) - data_median(has_data)), 2, 'omitnan');
    MAD = -1/(sqrt(2)*erfcinv(3/2)) * mean(MAD(:), "omitnan");
    thr = data_median + permute([-1 1],[3 1 2]).* 3.*MAD;
    id = sweep_dbz < thr(:,:,1) | sweep_dbz > thr(:,:,2);
    if ismember("outliers_continuity", parms.plot_opts) && i_elev == 1
        fig = bmm4d.plot.create_export_figure(600, 600);
        tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
        nexttile; hold on; axis tight; box on;
        bmm4d.plot.plot_sweep(sweep_dbz, gate.az, gate.range, coordinate_system="polar", ...
            max_range=100000, alpha_data=~isnan(sweep_dbz), show_colorbar=false);
        [x,y]=ind2sub(size(id),find(id)); plot(gate.az(y),gate.range(x),'.r');
        nexttile; hold on; box on;
        histogram(sweep_dbz)
        histogram(sweep_dbz(id));
        xlabel("Reflectivity")
        legend("All data before filtering", "Data filtered")
        nexttile;
        z = bmm4d.plot.plot_sweep(sweep_dbz, gate.az, gate.range, max_range=100000, show_colorbar=false);
        c_lim = [0 quantile(z(:),.999)];
        clim(c_lim); title("Before filter of spatial continuity")
        nexttile;
        tmp = sweep_dbz; tmp(id) = nan;
        bmm4d.plot.plot_sweep(tmp, gate.az, gate.range, max_range=100000, ...
            c_lim=c_lim, ttl="After filter of spatial continuity", show_colorbar=false);
        bmm4d.radar.filter.finish_plot(fig, parms, "_outliers_continuity.png")
    end
    sweep_dbz(id) = nan;
    gate.data = sweep_dbz;
    scan.gates(i_elev) = gate;
end
end
