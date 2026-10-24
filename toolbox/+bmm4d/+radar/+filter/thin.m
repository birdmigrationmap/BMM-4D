function scan = thin(scan, parms)
%THIN Keep one gate every data_res = [range tangential] metres.
%   Every k-th range bin is kept, and in each kept bin every m-th radial,
%   with m growing closer to the radar so that kept gates are about
%   data_res(2) apart; the others become missing.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    dens = gate.data;
    res_range = max(1,floor(parms.data_res(1)./median(diff(gate.range))));
    res_az = max(1,floor(parms.data_res(2)./(deg2rad(median(diff(gate.az)))*gate.range)));
    dens_dns = nan(size(dens));
    for i=1:res_range:size(dens_dns,1)
        dens_dns(i, 1:res_az(i):end) = dens(i, 1:res_az(i):end);
    end

    if ismember("downscale", parms.plot_opts) && i_elev == 1
        data0 = sweep2mat(scan.dz.sweeps(i_elev));
        dens0 = z_to_refl(idb(data0))/11;
        c_lim = [0 quantile(dens0(:),.99)];
        fig = bmm4d.plot.create_export_figure(900, 900);
        tiledlayout(3,1,'TileSpacing','tight','Padding','tight')
        nexttile; hold on;
        bmm4d.plot.plot_sweep(dens0, gate.az, gate.range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);
        nexttile; hold on;
        bmm4d.plot.plot_sweep(dens, gate.az, gate.range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);
        yline(parms.minRange,'-r')
        yline(parms.maxRange,'-r')
        nexttile; hold on;
        tmp = fillmissing2([dens_dns dens_dns(:,1)], 'linear');
        tmp = tmp(:,1:end-1);
        bmm4d.plot.plot_sweep(tmp, gate.az, gate.range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);
        bmm4d.radar.filter.finish_plot(fig, parms, "_downscale.png")
    end
    gate.data = dens_dns;
    scan.gates(i_elev) = gate;
end
end
