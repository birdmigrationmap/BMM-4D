function scan = smooth_median(scan, parms)
%SMOOTH_MEDIAN Moving median of density over smoothing (m) in range and azimuth.
%   The window spans smoothing metres along range, then the same arc
%   length along azimuth (more radials close to the radar), ignoring
%   missing gates.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    dens = gate.data;
    window_meter = parms.smoothing;

    % Smoothing along range
    window = window_meter / median(diff(gate.range));
    dens1 = smoothdata(dens, "movmedian", window, "omitnan");

    % Smoothing along azimuth, rows grouped by window length
    window = max(1,round(window_meter ./ (deg2rad(median(diff(gate.az)))*gate.range)));
    [G, ID] = findgroups(window);
    dens2 = dens1;
    for u=1:max(G)
        if ID(u)>1
            dens2(G==u,:) = smoothdata(dens1(G==u,:), 2, "movmedian", ID(u), "omitnan");
        end
    end

    if ismember("smoothing", parms.plot_opts) && i_elev == 1
        fig = bmm4d.plot.create_export_figure(900, 600);
        nexttile;
        imagesc(dens,'AlphaData',~isnan(dens));
        c_lim = [0 quantile(dens(:),.99)]; clim(c_lim);
        y_lim = [0 find(mean(isnan(dens),2)<.5,1,"last")]; ylim(y_lim);
        nexttile; imagesc(dens2,'AlphaData',~isnan(dens2)); clim(c_lim); ylim(y_lim);
        bmm4d.radar.filter.finish_plot(fig, parms, "_smoothing.png")
    end
    gate.data = dens2;
    scan.gates(i_elev) = gate;
end
end
