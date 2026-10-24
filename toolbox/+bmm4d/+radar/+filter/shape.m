function scan = shape(scan, parms)
%SHAPE Correct the azimuthal anisotropy of bird reflectivity (density).
%   Birds seen side-on reflect more than head- or tail-on. The orientation
%   axis is fitted to the azimuthal rhoHV pattern (sinusoid), then an
%   elliptical template with different front and back axes is fitted to
%   the smoothed azimuthal mean density of all sweeps. Each sweep is
%   divided by the template and rescaled to keep its mean density. The
%   template fit is stored in scan.shape_fit ([] when skipped: fewer than
%   10 azimuths with echo).

scan.shape_fit = [];
dens_elev = arrayfun(@(gate) gate.data, scan.gates, UniformOutput=false);
dmse = cell(numel(scan.gates),1);
rhmse = cell(numel(scan.gates),1);
for i_elev=1:numel(scan.gates)
    % Filter potential outlier by removing values below quantile .05 and above .95
    tmp = dens_elev{i_elev};
    tmp_q = quantile(tmp, [.05 .95]);
    tmp(tmp<tmp_q(1,:) | tmp>tmp_q(2,:))=nan;

    % Compute the mean
    tmp = mean(tmp, "omitnan");
    % To perform smoothing along azimuth, replicate data to have
    % cyclical behaviour
    w = 100;
    tmp_wrapped = [tmp(end-w+1:end) tmp tmp(1:w)];

    % Apply a moving median to remove outliers and a moving mean to
    % smooth the curve along azimuth
    dens_elev_sm = movmean(movmedian(tmp_wrapped,5),21);
    dens_elev_sm = dens_elev_sm(w+1:end-w);

    sweep_az = scan.gates(i_elev).az;

    % figure; hold on;
    % plot(sweep_az,quantile(dens_elev{i_elev}, [.2 .5 .8])');
    % plot(sweep_az,mean(dens_elev{i_elev}, "omitnan"));
    % plot(sweep_az,dens_elev_sm);

    dmse{i_elev} = [i_elev*ones(size(sweep_az)) sweep_az dens_elev_sm'];

    % Extract rho_HV from the sweep at the same elevation; dz and rh
    % sweep lists can differ (e.g. an extra low rh sweep).
    [~, i_elev_rh] = min(abs(scan.gates(i_elev).elev - [scan.rh.sweeps.elev]));
    [rh, ~, sweep_az] = sweep2mat(scan.rh.sweeps(i_elev_rh));

    % Remove when we don't have enough bird
    tmp = min(size(rh,1), size(dens_elev{i_elev},1));
    rh(dens_elev{i_elev}(1:tmp,:)<10)=nan;

    % Compute the average rh along range weighted by density
    w = sum(dens_elev{i_elev}(1:tmp,:) .* ~isnan(rh(1:tmp,:)), "omitnan");
    tmp = sum(dens_elev{i_elev}(1:tmp,:) .* rh(1:tmp,:), "omitnan") ./ w;
    tmp(isnan(tmp))=0;

    rhmse{i_elev} = [i_elev*ones(size(sweep_az)) sweep_az tmp' w'];

end
dmse = cell2mat(dmse);
rhmse = cell2mat(rhmse);

% The shape model needs echo in enough azimuths; otherwise keep the
% densities uncorrected (e.g. a volume almost entirely removed as rain).
if nnz(rhmse(:,4) > 0) < 10 || nnz(dmse(:,3) > 0) < 10
    warning("read_rd:ShapeSkipped", "%s: too little echo for the shape correction.", parms.filename)
else

    % figure; nexttile; box on; scatter(dmse(:,2), dmse(:,3), 10, dmse(:,1), 'filled')
    % nexttile; box on; scatter(rhmse(:,2), rhmse(:,3), 10, rhmse(:,1), 'filled')

    % Perform nonlinear curve fitting
    mdl2 = fit(deg2rad(rhmse(:,2)), rhmse(:,3), fittype('A*sin(x + phi) + C', 'independent', 'x'), StartPoint=[range(rhmse(:,3))/2, mean(rhmse(:,3)), 0], Lower=[0,0,-2*pi], Upper=[1,1,2*pi], Weights = rhmse(:,4));

    % radius/distance from center of an ellipse
    % R = @(a,b,theta,x) a*b./sqrt(a^2*sind(x-theta).^2+b^2*cosd(x-theta).^2);
    % Rfit = @(alpha, frontback, xm, x) xm .* R(alpha, 1/alpha, rad2deg(-mdl2.phi+pi), x);

    % Upgrade to have different small axis
    b_func = @(xt, b_left, b_right) b_left + (b_right - b_left) * (mod(xt,360) >= 180);

    R = @(a, b_left, b_right, theta, x) a .* b_func(x - theta, b_left, b_right) ./ sqrt(a^2 .* sind(x - theta).^2 + b_func(x - theta, b_left, b_right).^2 .* cosd(x - theta).^2);
    Rfit = @(alpha, beta, xm, x) xm .* R(alpha, beta*1/alpha, (1-beta)*1/alpha, rad2deg(-mdl2.phi+pi), x);
    % figure; x=0:360; plot(x, R(1, 0.8, 0.6, 0, x))


    % Fit the anisotropy model on azimuth sectors with positive smoothed density.
    id = dmse(:,3) > 0;

    % Perform fit
    mdl = fit(dmse(id,2), dmse(id,3), fittype(Rfit), StartPoint=[0 .5 mean(dmse(id,3),"omitnan")], Lower=[.1 .1 0], Upper=[5 .9 10000], Weights = dmse(id,3));

    % figure; hold on; plot(mdl, dmse(id,2), dmse(id,3));

    dens_elev_corr = dens_elev;

    % Apply the correction to each elevation angle separately
    for i_elev=1:numel(scan.gates)
        sweep_az = scan.gates(i_elev).az;

        % Compute the mean density
        mdens = mean(dens_elev_corr{i_elev}(:),'omitnan');

        % Apply correction
        dens_elev_corr{i_elev} = dens_elev_corr{i_elev} ./ R(mdl.alpha, mdl.beta*1/mdl.alpha, (1-mdl.beta)*1/mdl.alpha, rad2deg(-mdl2.phi+pi), sweep_az)';

        % Normalize to keep the overall same density
        dens_elev_corr{i_elev} = dens_elev_corr{i_elev} / mean(dens_elev_corr{i_elev}(:),'omitnan') .* mdens;
    end


    if ismember("shape", parms.plot_opts)
        i_elev=1;
        [rh, sweep_range, sweep_az] = sweep2mat(scan.rh.sweeps(i_elev));

        fig = bmm4d.plot.create_export_figure(2400, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
        theme(fig,"light")
        tiledlayout(3,6,'TileSpacing','tight','Padding','tight');

        % Plot the data and the fitted curve
        ax_summary_rh=nexttile([1 3]); hold on; grid(gca, "on"); box on;
        id = rhmse(:,3)>0;
        scatter(deg2rad(rhmse(id,2)),rhmse(id,3),rhmse(id,4)/500, rhmse(id,1),'filled',MarkerFaceAlpha=0.3)
        x_fit = linspace(0, 2*pi, 361)';
        plot(x_fit, mdl2(x_fit), "k-", DisplayName="Fitted sinusoid");
        xline(mod(-mdl2.phi+3/2*pi,2*pi), LineWidth=1.5, Label="Orientation axis", DisplayName="Orientation axis")
        xlabel('Azimuth (rad)');  ylabel('\rho_{HV}');
        lg_rh=legend('Azimuthal summary', 'Sinusoid', 'Orientation axis', 'Location', 'southoutside');
        axis tight; xlim([0 2*pi]); c=colorbar;c.Label.String="Elevation (°)";cb_angle_rh=c;

        ax_summary_density=nexttile([1 3]);hold on; grid(gca, "on"); box on;
        scatter(dmse(:,2),dmse(:,3),dmse(:,3)/10, dmse(:,1),'filled',MarkerFaceAlpha=0.3)
        x_fit = linspace(0, 360, 361)';
        plot(x_fit, mdl(x_fit), "k-", DisplayName="Fitted anisotropy template");
        xline(rad2deg(mod(-mdl2.phi+3/2*pi,2*pi)), LineWidth=1.5)
        xlabel('Azimuth (deg)');  ylabel('Bird density (birds km^{-3})');
        lg_density=legend('Azimuthal summary', 'Anisotropy template', 'Location', 'southoutside'); axis tight;
        xlim([0 360]); c=colorbar;c.Label.String="Elevation (°)";cb_angle_density=c;


        ax = nexttile([2 2]);ax_map_rh=ax;
        bmm4d.plot.plot_sweep(rh, sweep_az, sweep_range, max_range=100000, show_colorbar=true, ...
            ttl="\rho_{HV} and orientation axis");
        plot([0 -50000.*cos(mdl2.phi)], [0 -50000.*sin(mdl2.phi)], "k", LineWidth=4)
        xticks([-100 0 100]*1000); xticklabels([-100 0 100]); yticks([-100 0 100]*1000); yticklabels([-100 0 100]); xlabel("x (km)"); ylabel("y (km)")
        % clim([0 quantile(z(:),.999)])
        colormap(ax, bmm4d.plot.colormap_wct("rh"))

        sweep_range = scan.gates(i_elev).range;
        sweep_az = scan.gates(i_elev).az;

        %x=0:365;
        %polarplot(deg2rad(dmse(:,2)), dmse(:,3), '.',Color=[.8 .8 .8]); hold on;
        % polarplot(deg2rad(dmse(id,2)), dmse(id,3), '.k'); hold on;
        %polarplot(deg2rad(dmse(:,2)), dmse(:,3) ./ R(sqrt(mdl.ani), 1/sqrt(mdl.ani), mdl.theta, dmse(:,2)), '.r');
        %polarplot(deg2rad(x), Rfit(mdl.ani, mdl.theta, mdl.xm, x), 'r');
        %legend("Before shape correction", "After shape correction", "Correction template", 'Location', 'southwest')
        %ax = gca;
        %ax.ThetaZeroLocation = 'top';
        %ax.ThetaDir = 'clockwise';

        ax_map_before=nexttile([2 2]);
        z = bmm4d.plot.plot_sweep(dens_elev{i_elev}, sweep_az, sweep_range, max_range=100000, show_colorbar=false);
        c_lim = [0 quantile(z(:),.99)];
        xticks([-100 0 100]*1000); xticklabels([-100 0 100]); yticks([-100 0 100]*1000); yticklabels([-100 0 100]); xlabel("x (km)"); ylabel("y (km)")
        clim(c_lim); title("Before shape correction")

        ax_map_after=nexttile([2 2]);
        bmm4d.plot.plot_sweep(dens_elev_corr{i_elev}, sweep_az, sweep_range, max_range=100000, ...
            c_lim=c_lim, ttl="After shape correction", show_colorbar=true);
        xticks([-100 0 100]*1000); xticklabels([-100 0 100]); yticks([-100 0 100]*1000); yticklabels([-100 0 100]); xlabel("x (km)"); ylabel("y (km)")
        cb_map_rh=colorbar(ax_map_rh);cb_map_density=colorbar(ax_map_after);
        shape_axes=[ax_summary_rh ax_summary_density ax_map_rh ax_map_before ax_map_after];
        shape_positions=[0.06 0.65 0.34 0.29;0.55 0.65 0.34 0.29;0.065 0.13 0.20 0.34;0.39 0.13 0.20 0.34;0.715 0.13 0.20 0.34];
        for i_axis=1:numel(shape_axes)
            shape_axes(i_axis).Parent=fig;shape_axes(i_axis).Position=shape_positions(i_axis,:);
        end
        cb_angle_rh.Position=[0.415 0.65 0.012 0.29];cb_angle_density.Position=[0.91 0.65 0.012 0.29];
        cb_map_rh.Position=[0.28 0.13 0.012 0.34];cb_map_density.Position=[0.93 0.13 0.012 0.34];
        cb_map_density.Label.String="";
        lg_rh.Units="normalized";lg_rh.Position=[0.08 0.52 0.29 0.09];lg_rh.Box="off";
        lg_density.Units="normalized";lg_density.Position=[0.57 0.52 0.29 0.09];lg_density.Box="off";
        bmm4d.plot.style_paper_figure(fig,9.5)
        lg_rh.String=["Data" "Sinusoid" "Orientation"];lg_rh.NumColumns=3;lg_rh.FontSize=8;
        lg_rh.Position=[0.055 0.52 0.37 0.035];
        lg_density.String=["Data" "Anisotropy template"];lg_density.NumColumns=2;lg_density.FontSize=8;
        lg_density.Position=[0.555 0.52 0.34 0.035];
        bmm4d.radar.filter.finish_plot(fig, parms, "_shape.png")
    end
    scan.shape_fit = mdl;
    dens_elev = dens_elev_corr;
end

for i_elev = 1:numel(scan.gates)
    scan.gates(i_elev).data = dens_elev{i_elev};
end
end
