
function [radar_scans, shape_fit] = read_rd(scans_filepath, varargin)
%READ_RD Read and process weather radar scan files.
%   [RADAR_SCANS, SHAPE_FIT] = READ_RD(SCANS_FILEPATH, ...) reads and
%   processes the radar scan files specified in SCANS_FILEPATH (cell array
%   of strings or string array). Returns a cell array RADAR_SCANS with
%   processed radar data and SHAPE_FIT with shape correction models (if
%   enabled).
%
%   Optional name-value pair arguments:
%     'removeRain'    - Method to remove rain ('rho', 'mistnet', 'none'). Default: 'rho'.
%     'rainRhoThreshold' - rhoHV threshold for rain removal when removeRain='rho'. Default: 0.95.
%     'removeOutliers' - MAD value for outlier removal in reflectivity units. Default: 3.
%     'thr_dens'      - Max bird density (bird/km^3) accepted. Default: 1e9.
%     'shape'         - Apply shape correction (true/false). Default: true.
%     'smoothing'     - Smoothing window (meters) for bird density. Default: 2000.
%     'limitElev'     - Max elevation scan to consider (degrees). Default: 5.
%     'minRange'      - Minimum range from radar (meters). Default: 5000.
%     'maxRange'      - Maximum range from radar (meters). Default: Inf.
%     'maxAlt'        - Maximum altitude to consider (meters). Default: 3000.
%     'data_res'      - Data resolution [range, azimuth*range] in meters. Default: [2000 2000].
%     'waitbar'       - Show waitbar (true/false). Default: true.
%     'plot'          - Plot options (cell array of strings). Default: "".
%     'exportDir'     - Output folder for exported figures. Default: "export/read_rd".
%     'parallel'      - Clean volumes in parallel (no plot/waitbar). Default: true.
%     'maxWorkers'    - Maximum volumes cleaned at the same time. Default: 6.
%
%   Example:
%     [radar_scans, shape_fit] = read_rd({'scan1.mat','scan2.mat'}, 'removeRain', 'rho', 'plot', {'rain','export'}, 'exportDir', fullfile("03_simulator_fit", "figures", "read_rd"));

p = inputParser;
addParameter(p,'removeRain','rho'); % Method to remove rain ('rho', 'mistnet', 'none')
addParameter(p,'rainRhoThreshold', 0.95); % rhoHV threshold for rain removal
addParameter(p,'removeOutliers', 3); % MAD value for outlier removal in reflectivity units
addParameter(p,'thr_dens', 1e9); % Max density (bird/km^3) accepted
addParameter(p,'shape', true); % Apply shape correction
addParameter(p,'smoothing', 2000); % Smoothing window (meters)
addParameter(p,'limitElev', 5); % Max elevation scan to consider (degrees)
addParameter(p,'minRange', 5000); % Minimum range from radar (meters)
addParameter(p,'maxRange', Inf); % Maximum range from radar (meters)
addParameter(p,'maxAlt', 3000); % Maximum altitude to consider (meters)
addParameter(p,'data_res', [2000 2000]); % Data resolution [range, azimuth*range] in meters
addParameter(p,'waitbar', false); % Show waitbar
addParameter(p,'plot', ""); % Plot options
addParameter(p,'exportDir', fullfile("export", "read_rd")); % Output folder for exported figures
addParameter(p,'parallel', true); % Clean volumes in parallel when no plot or waitbar is requested
addParameter(p,'maxWorkers', 6); % Maximum number of volumes cleaned at the same time

parse(p,varargin{:});
parms = p.Results;
plot_opts = string(parms.plot);
export_dir = parms.exportDir;
if ismember("export", plot_opts) && ~exist(export_dir, "dir")
    mkdir(export_dir);
end

radar_scans = cell(size(scans_filepath));
shape_fit = cell(size(scans_filepath));

if parms.waitbar
    q = parallel.pool.DataQueue;
    afterEach(q,@UpdateWaitbar);
    i_s_wb = 0;
    N = numel(scans_filepath);
    wb = waitbar(0,"Clean and load weather radar data: "+ num2str(0) + "/"+num2str(N));
end

% Check that all files exist
for i_s = 1:numel(scans_filepath)
    if(~exist(scans_filepath(i_s),'file'))
        error("scans_filepath does not exist: %s", scans_filepath(i_s))
    end
end

% Volumes are independent: clean them in parallel (at most maxWorkers at
% a time to bound memory) unless plots or the waitbar need the serial loop.
use_parallel = parms.parallel && ~parms.waitbar && all(plot_opts == "") && numel(scans_filepath) > 1;
if use_parallel
    parfor (i_s = 1:numel(scans_filepath), parms.maxWorkers)
        [radar_scans{i_s}, shape_fit{i_s}] = clean_scan(scans_filepath(i_s), parms, plot_opts, export_dir);
    end
else
    for i_s = 1:numel(scans_filepath)
        [radar_scans{i_s}, shape_fit{i_s}] = clean_scan(scans_filepath(i_s), parms, plot_opts, export_dir);
        if parms.waitbar
            send( q, i_s );
        end
    end
end
    function UpdateWaitbar(~)
        waitbar(i_s_wb/N,wb,"Clean and load weather radar data: "+ num2str(i_s_wb) + "/"+num2str(N));
        i_s_wb = i_s_wb + 1;
    end
if parms.waitbar && isvalid(wb)
    close(wb);
end
end


function [radar_scan, scan_shape_fit] = clean_scan(scan_file, parms, plot_opts, export_dir)
%CLEAN_SCAN Read and clean one radar volume (the body of the read_rd loop).

% Reference maps are loaded once per MATLAB process (client or worker).
persistent bb radar_sites
if isempty(bb)
    bb = load("data/reference/beam_blockage/beam_blockage");
end
if isempty(radar_sites)
    radar_sites = readtable("data/reference/radar.csv", TextType="string");
end
scan_shape_fit = [];

% Extract filename for figures
[~, filename] = fileparts(scan_file);

radar_scan = rsl2mat_wrap(scan_file, struct('max_elev', parms.limitElev));

% The file header height is the ground elevation of the site, but beam
% heights (slant2ground) are relative to the antenna. Use the antenna
% height (ground + tower) from radar.csv, as for the beam-blockage maps.
site_id = find(radar_sites.name == string(radar_scan.station), 1);
if isempty(site_id)
    warning("read_rd:UnknownStation", "%s not in radar.csv; using the file height.", radar_scan.station)
else
    radar_scan.height = radar_sites.height(site_id);
end

bbr_f = griddedInterpolant({bb.bb_range, bb.bb_az}, bb.bb(:,:,radar_scan.station==bb.bb_radar),'nearest','nearest');

% Get unique sweeps
radar_scan.dz.sweeps = unique_elev_sweeps(radar_scan, 'dz');
radar_scan.rh.sweeps = unique_elev_sweeps(radar_scan, 'rh');

% Limit elevation angle
radar_scan.dz.sweeps = radar_scan.dz.sweeps([radar_scan.dz.sweeps.elev]< parms.limitElev);
radar_scan.rh.sweeps = radar_scan.rh.sweeps([radar_scan.rh.sweeps.elev]< parms.limitElev);

% Build rain masks on the sorted polar grids returned by sweep2mat.
rain_masks = cell(numel(radar_scan.dz.sweeps), 1);
if parms.removeRain=="rho"
    for i_elev=1:numel(radar_scan.dz.sweeps)
        [min_elev,i_elev_rh] = min(abs(radar_scan.dz.sweeps(i_elev).elev-vertcat(radar_scan.rh.sweeps.elev)));
        assert(min_elev<.01)
        [rh_data, ~, ~] = sweep2mat(radar_scan.rh.sweeps(i_elev_rh));
        dz_data = sweep2mat(radar_scan.dz.sweeps(i_elev));
        rain_masks{i_elev} = true(size(dz_data));
        min_nbins = min(size(rh_data, 1), size(dz_data, 1));
        rain_masks{i_elev}(1:min_nbins,:) = rh_data(1:min_nbins,:) > parms.rainRhoThreshold;
    end
elseif parms.removeRain=="mistnet"
    [preds, classes] = mistnet_polar_radar(radar_scan);
    for i_elev=1:numel(radar_scan.dz.sweeps)
        rain_masks{i_elev} = preds{i_elev}==classes.rain;
    end
elseif parms.removeRain=="none"
    for i_elev=1:numel(radar_scan.dz.sweeps)
        rain_masks{i_elev} = false(size(sweep2mat(radar_scan.dz.sweeps(i_elev))));
    end
else
    warning("Unknown rain removal method: %s", parms.removeRain)
end

if ismember("rain", plot_opts)
    [rain_dbz, rain_range, rain_az] = sweep2mat(radar_scan.dz.sweeps(1));
    rain_dbz_filtered = rain_dbz;
    rain_dbz_filtered(rain_masks{1}) = nan;
    rain_range_limit = min(parms.maxRange, max(rain_range));

    fig = create_export_figure(1200, 450);
    tiledlayout(1, 3, "TileSpacing", "tight", "Padding", "tight")
    nexttile
    plot_sweep(rain_dbz, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Reflectivity before rain removal", show_colorbar=false);

    nexttile
    plot_sweep(rain_masks{1}, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Rain mask: " + string(parms.removeRain), show_colorbar=false);

    nexttile
    plot_sweep(rain_dbz_filtered, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Reflectivity after rain removal", show_colorbar=false);

    if ismember("export", plot_opts)
        export_fixed_figure(fig, fullfile(export_dir, filename + "_rain.png"), Resolution=300)
        close(fig)
    else
        set(fig, "Visible", "on")
    end
end

% Remove unused fields
radar_scan = rmfield(radar_scan,{'sw','dr','ph','spulse','lpulse','constants','params'});
% radar_scan = rmfield(radar_scan,{'rh'});

% Create bird density structure
radar_scan.dens = struct();
radar_scan.dens.type = "Bird density";
radar_scan.dens.sweeps = struct([]);

% Basic cleaning by elevation and conversion to density
dens_elev = cell(numel(radar_scan.dz.sweeps), 1);
for i_elev=1:numel(radar_scan.dz.sweeps)

    % Extract reflectivity in dBZ for one sweep on a regular polar g.
    [sweep_dbz, sweep_range, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));
    sweep_dbz(rain_masks{i_elev}) = nan;

    % Discard values that overlap with the precomputed beam blockage
    % map. Overlap is considered if the ground is within the half-beam
    % width
    sweep_dbz(bbr_f({sweep_range, sweep_az}) > (radar_scan.dz.sweeps(i_elev).elev-radar_scan.dz.sweeps(i_elev).vert_half_bw)) = nan;

    if parms.removeOutliers>0

        % 1. Threshold high value
        % density of bird too high
        sweep_density = z_to_refl(idb(sweep_dbz)) / 11; % inverse decibel transform --> bird density [bird/km^3]
        % parms.thr_dens = z_to_refl(idb(40))/11; % 3000 bird per km^3
        if ismember("outliers_threshold", plot_opts) && i_elev==1
            fig = create_export_figure(1200, 550);
            tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
            nexttile;
            z = plot_sweep(sweep_dbz, sweep_az, sweep_range, max_range=p.Results.maxRange, show_colorbar=false);
            c_lim = [0 quantile(z(:),.999)];
            clim(c_lim); title("Before threshold of high value")
            nexttile;
            tmp = sweep_dbz; tmp(sweep_density > parms.thr_dens) = nan;
            plot_sweep(tmp, sweep_az, sweep_range, max_range=p.Results.maxRange, ...
                c_lim=c_lim, ttl="After threshold of high value " + num2str(parms.thr_dens,1) + " bird/km^3", show_colorbar=false);
            if ismember("export", plot_opts)
                export_fixed_figure(fig, fullfile(export_dir, filename + "_outliers_threshold.png"), Resolution=300)
                close(fig)
            else
                set(fig, "Visible", "on")
            end
        end
        sweep_dbz(sweep_density > parms.thr_dens) = nan;

        % 2. Threshold on the marginal distribution.
        % grubbs marginal distribution
        % tmp = isoutlier(sweep_dbz(:),'grubbs','ThresholdFactor',1);
        % n = sum(~isnan(sweep_dbz(:)));
        % alpha=.05;
        % tval = tinv(alpha./n/2,repmat(n-2,1,numel(alpha))).^2;
        % G = (n-1)./sqrt(n) .* sqrt(tval./(n-2+tval));
        % thr = mean(sweep_dbz(:),"omitnan")+[-1 1]*G.*std(sweep_dbz(:),"omitnan")

        mdata = median(sweep_dbz(:),'omitnan');
        MAD = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(sweep_dbz(:) - mdata), 'omitnan');
        thr = mdata + [-1 1]* parms.removeOutliers*MAD;
        id = sweep_dbz < thr(1) | sweep_dbz > thr(2);
        if ismember("outliers_mad", plot_opts) && i_elev==1
            fig = create_export_figure(1400, 500);
            tiledlayout(1,3,'TileSpacing','tight','Padding','tight');
            nexttile; histogram(sweep_dbz(:));
            title("Histogram of all data")
            xline(mdata,'-k',"Median", 'LineWidth',2); xline(thr,'--k',"Outlier threshold", 'LineWidth',2)
            xlabel("reflectivity");
            nexttile;
            z = plot_sweep(sweep_dbz, sweep_az, sweep_range, max_range=p.Results.maxRange, show_colorbar=false);
            c_lim = [0 quantile(z(:),.999)];
            clim(c_lim); title("Before threshold of MAD")
            nexttile;
            tmp = sweep_dbz; tmp(id) = nan;
            plot_sweep(tmp, sweep_az, sweep_range, max_range=p.Results.maxRange, ...
                c_lim=c_lim, ttl="After threshold of MAD between " + num2str(round(thr(1))) + " and " + num2str(round(thr(2))) + " bird/km^3", show_colorbar=false);
            if ismember("export", plot_opts)
                export_fixed_figure(fig, fullfile(export_dir, filename + "_outliers_mad.png"), Resolution=300)
                close(fig)
            else
                set(fig, "Visible", "on")
            end
        end
        sweep_dbz(id) = nan;

        % 3. Spatial continuity (moving median)
        sweep_dbz(end-1:end,:) = nan;
        dd = repmat(sweep_dbz,1,1,5*5); u=1;
        for i=-2:2
            for j=-2:2
                dd(:,:,u) = circshift(sweep_dbz,[i j]);
                u=u+1;
            end
        end
        data_median = median(dd,3,'omitnan');
        data_median(end-1:end,:)=nan;

        MAD = -1/(sqrt(2)*erfcinv(3/2))*median(abs(dd-data_median),3,'omitnan');
        MAD = mean(MAD(:), "omitnan");
        thr = data_median + permute([-1 1],[3 1 2]).* 3.*MAD;
        id = sweep_dbz < thr(:,:,1) | sweep_dbz > thr(:,:,2);
        if ismember("outliers_continuity", plot_opts)  && i_elev==1

            fig = create_export_figure(600, 600);
            tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

            nexttile; hold on; axis tight; box on;
            plot_sweep(sweep_dbz, sweep_az, sweep_range, coordinate_system="polar", ...
                max_range=100000, alpha_data=~isnan(sweep_dbz), show_colorbar=false);
            [x,y]=ind2sub(size(id),find(id)); plot(sweep_az(y),sweep_range(x),'.r');

            nexttile; hold on; box on;
            histogram(sweep_dbz)
            histogram(sweep_dbz(id));
            xlabel("Reflectivity")
            legend("All data before filtering", "Data filtered")


            nexttile;
            z = plot_sweep(sweep_dbz, sweep_az, sweep_range, max_range=100000, show_colorbar=false);
            c_lim = [0 quantile(z(:),.999)];
            clim(c_lim); title("Before filter of spatial continuity")
            nexttile;
            tmp = sweep_dbz; tmp(id) = nan;
            plot_sweep(tmp, sweep_az, sweep_range, max_range=100000, ...
                c_lim=c_lim, ttl="After filter of spatial continuity", show_colorbar=false);

            if ismember("export", plot_opts)
                export_fixed_figure(fig, fullfile(export_dir, filename + "_outliers_continuity.png"), Resolution=300)
                close(fig)
            else
                set(fig, "Visible", "on")
            end
        end

        sweep_dbz(id) = nan;

    end

    % Set nodata pulse volumes to 0 reflectivity (-inf on decibel scale)
    % sweep_dbz(isnan(sweep_dbz)) = -inf;

    % Convert the cleaned dBZ sweep into linear reflectivity and density.
    sweep_reflectivity = idb(sweep_dbz);    % inverse decibel transform --> reflectivity factor
    sweep_reflectivity = z_to_refl(sweep_reflectivity); % reflectivity factor (mm^6/m^3) --> reflectivity (cm^2/km^3)
    dens_elev{i_elev} = sweep_reflectivity / 11; % reflectivity --> bird density [bird/km^3]
end

% Perform shape correction for all elevations together
if parms.shape

    % Compute a smoothed mean density profile along azimuth for each elevation angle
    dmse = cell(numel(radar_scan.dz.sweeps),1);
    rhmse = cell(numel(radar_scan.dz.sweeps),1);
    for i_elev=1:numel(radar_scan.dz.sweeps)
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

        [~, ~, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));

        % figure; hold on;
        % plot(sweep_az,quantile(dens_elev{i_elev}, [.2 .5 .8])');
        % plot(sweep_az,mean(dens_elev{i_elev}, "omitnan"));
        % plot(sweep_az,dens_elev_sm);

        dmse{i_elev} = [i_elev*ones(size(sweep_az)) sweep_az dens_elev_sm'];

        % Extract rho_HV from the sweep at the same elevation; dz and rh
        % sweep lists can differ (e.g. an extra low rh sweep).
        [~, i_elev_rh] = min(abs(radar_scan.dz.sweeps(i_elev).elev - [radar_scan.rh.sweeps.elev]));
        [rh, ~, sweep_az] = sweep2mat(radar_scan.rh.sweeps(i_elev_rh));

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
    for i_elev=1:numel(radar_scan.dz.sweeps)
        [~, ~, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));

        % Compute the mean density
        mdens = mean(dens_elev_corr{i_elev}(:),'omitnan');

        % Apply correction
        dens_elev_corr{i_elev} = dens_elev_corr{i_elev} ./ R(mdl.alpha, mdl.beta*1/mdl.alpha, (1-mdl.beta)*1/mdl.alpha, rad2deg(-mdl2.phi+pi), sweep_az)';

        % Normalize to keep the overall same density
        dens_elev_corr{i_elev} = dens_elev_corr{i_elev} / mean(dens_elev_corr{i_elev}(:),'omitnan') .* mdens;
    end


    if ismember("shape", plot_opts)
        i_elev=1;
        [rh, ~, sweep_az] = sweep2mat(radar_scan.rh.sweeps(i_elev));

        fig = create_export_figure(2400, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
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
        plot_sweep(rh, sweep_az, sweep_range, max_range=100000, show_colorbar=true, ...
            ttl="\rho_{HV} and orientation axis");
        plot([0 -50000.*cos(mdl2.phi)], [0 -50000.*sin(mdl2.phi)], "k", LineWidth=4)
        xticks([-100 0 100]*1000); xticklabels([-100 0 100]); yticks([-100 0 100]*1000); yticklabels([-100 0 100]); xlabel("x (km)"); ylabel("y (km)")
        % clim([0 quantile(z(:),.999)])
        colormap(ax, colormap_wct("rh"))

        [~, sweep_range, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));

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
        z = plot_sweep(dens_elev{i_elev}, sweep_az, sweep_range, max_range=100000, show_colorbar=false);
        c_lim = [0 quantile(z(:),.99)];
        xticks([-100 0 100]*1000); xticklabels([-100 0 100]); yticks([-100 0 100]*1000); yticklabels([-100 0 100]); xlabel("x (km)"); ylabel("y (km)")
        clim(c_lim); title("Before shape correction")

        ax_map_after=nexttile([2 2]);
        plot_sweep(dens_elev_corr{i_elev}, sweep_az, sweep_range, max_range=100000, ...
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
        style_paper_figure(fig,9.5)
        lg_rh.String=["Data" "Sinusoid" "Orientation"];lg_rh.NumColumns=3;lg_rh.FontSize=8;
        lg_rh.Position=[0.055 0.52 0.37 0.035];
        lg_density.String=["Data" "Anisotropy template"];lg_density.NumColumns=2;lg_density.FontSize=8;
        lg_density.Position=[0.555 0.52 0.34 0.035];
        if ismember("export", plot_opts)
            export_fixed_figure(fig, fullfile(export_dir, filename + "_shape.png"), Resolution=300)
            close(fig)
        else
            set(fig, "Visible", "on")
        end
    end
     scan_shape_fit = mdl;
    dens_elev = dens_elev_corr;
end



for i_elev=1:numel(radar_scan.dz.sweeps)

    % retrieve elevation values needed
    [~, sweep_range, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));
    dens = dens_elev{i_elev};

    if parms.smoothing>0

        smooth_method = "movmedian";

        % Use smoothing as the std of the gaussian kernel in m
        window_meter = parms.smoothing;

        % Smoothing along range
        window = window_meter / median(diff(sweep_range));

        % Option Gaussian filter manual accounting for nan
        % assuming window = std
        % ss = 2*ceil(2*window)+1;
        % W = kernelGaussian1D(ss, sigma_unit);
        % W = reshape(W(:),1,1,[]);
        % dd = nan(size(dens,1), size(dens,2), ss);
        % for u=1:ss
        %     dd(:,:,u) = circshift(dens,u-((ss-1)/2),1);
        % end
        % dens1 = sum(dd.*W,3,'omitnan') ./ sum(~isnan(dd).*W,3);
        % % Set as nan is more than 50% of the gaussian kernel is also nan
        % dens1(sum(W.*~isnan(dd),3)<.5) = nan;
        % % dens1(isnan(dens)) = nan;
        %
        % Option 2: Moving median much faster!
        dens1 = smoothdata(dens, smooth_method, window, "omitnan");


        % Smoothing along azimuth
        window = max(1,round(window_meter ./ (deg2rad(median(diff(sweep_az)))*sweep_range)));
        [G, ID] = findgroups(window);

        % Option Gaussian filter manual accounting for nan
        % ss = 2*ceil(2*window)+1;
        % W = kernelGaussian1D(max(ss), sigma_unit);
        % W = permute(W,[1 3 2]);
        % dd = nan(size(dens1,1), size(dens1,2), max(ss));
        % for u=1:max(ss)
        %     dd(:,:,u) = circshift(dens1,u-((max(ss)-1)/2),2);
        % end
        % dens2 = sum(dd.*W,3,'omitnan') ./ sum(~isnan(dd).*W,3);
        % dens2(sum(W.*~isnan(dd),3)<.5) = nan;

        % Option 2: Moving median much faster!
        dens2 = dens1;
        for u=1:max(G)
            if ID(u)>1
                dens2(G==u,:) = smoothdata(dens1(G==u,:), 2, smooth_method, ID(u), "omitnan");
            end
        end

        if ismember("smoothing", plot_opts) && i_elev==1
            fig = create_export_figure(900, 600);
            nexttile;
            imagesc(dens,'AlphaData',~isnan(dens));
            c_lim = [0 quantile(dens(:),.99)]; clim(c_lim);
            y_lim = [0 find(mean(isnan(dens),2)<.5,1,"last")]; ylim(y_lim);
            nexttile; imagesc(dens2,'AlphaData',~isnan(dens2)); clim(c_lim); ylim(y_lim);

            if ismember("export", plot_opts)
                export_fixed_figure(fig, fullfile(export_dir, filename + "_smoothing.png"), Resolution=300)
                close(fig)
            else
                set(fig, "Visible", "on")
            end
        end

        dens = dens2;
    end


    % Remove datapoint too close or too high altitude
    dens(sweep_range<parms.minRange,:)=nan;
    dens(sweep_range>parms.maxRange,:)=nan;
    [ ~, sweep_height ] = slant2ground( sweep_range, radar_scan.dz.sweeps(i_elev).elev );
    dens(sweep_height>parms.maxAlt,:)=nan;

    % Downscale
    res_range = max(1,floor(parms.data_res(1)./median(diff(sweep_range))));
    res_az = max(1,floor(parms.data_res(2)./(deg2rad(median(diff(sweep_az)))*sweep_range)));

    dens_dns = nan(size(dens));
    for i=1:res_range:size(dens_dns,1)
        dens_dns(i, 1:res_az(i):end) = dens(i, 1:res_az(i):end);
    end

    if ismember("downscale", plot_opts) && i_elev==1

        data0 = sweep2mat(radar_scan.dz.sweeps(i_elev));
        dens0 = z_to_refl(idb(data0))/11;
        c_lim = [0 quantile(dens0(:),.99)];

        fig = create_export_figure(900, 900);
        tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

        nexttile; hold on;
        plot_sweep(dens0, sweep_az, sweep_range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);

        nexttile; hold on;
        plot_sweep(dens, sweep_az, sweep_range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);
        yline(parms.minRange,'-r')
        yline(parms.maxRange,'-r')

        nexttile; hold on;
        %imagescnan(sweep_az, sweep_range,dens_dns);
        tmp = fillmissing2([dens_dns dens_dns(:,1)], 'linear');
        tmp = tmp(:,1:end-1);
        plot_sweep(tmp, sweep_az, sweep_range, coordinate_system="polar", ...
            max_range=parms.maxRange, c_lim=c_lim, show_colorbar=false);
        % tmp = find(~isnan(dens_dns));
        %[tmp1, tmp2] = meshgrid(sweep_az, sweep_range);
        % plot(tmp1(tmp), tmp2(tmp),'.k')

        if ismember("export", plot_opts)
            export_fixed_figure(fig, fullfile(export_dir, filename + "_downscale.png"), Resolution=300)
            close(fig)
        else
            set(fig, "Visible", "on")
        end
    end

    radar_scan.dens.sweeps(i_elev).data = dens_dns;
    radar_scan.dens.sweeps(i_elev).az = sweep_az;
    radar_scan.dens.sweeps(i_elev).range = sweep_range;
    radar_scan.dens.sweeps(i_elev).elev = radar_scan.dz.sweeps(i_elev).elev;

    % some additional data to keep useful for e.g. sweep2mat
    radar_scan.dens.sweeps(i_elev).range_bin1 = radar_scan.dz.sweeps(i_elev).range_bin1;
    radar_scan.dens.sweeps(i_elev).gate_size = radar_scan.dz.sweeps(i_elev).gate_size;
    radar_scan.dens.sweeps(i_elev).nbins = radar_scan.dz.sweeps(i_elev).nbins;
    % Need to rewrite azim_v because sweep2mat reorders data
    % for increasing azimuth
    radar_scan.dens.sweeps(i_elev).azim_v = sweep_az;
    radar_scan.dens.sweeps(i_elev).prf = radar_scan.dz.sweeps(i_elev).prf;
end
radar_scan = rmfield(radar_scan,{'dz'});
end
