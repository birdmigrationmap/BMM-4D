
function [radar_scans, shape_fit] = read_rd(scans_filepath, varargin)
%READ_RD Read and process weather radar scan files.
%   [RADAR_SCANS, SHAPE_FIT] = READ_RD(SCANS_FILEPATH, ...) reads and
%   processes the radar scan files specified in SCANS_FILEPATH (cell array
%   of strings or string array). Returns a cell array RADAR_SCANS with
%   processed radar data and SHAPE_FIT with shape correction models (if
%   enabled).
%
%   Optional name-value pair arguments:
%     'removeRain'    - Method to remove rain ('rho', 'dr', 'mistnet', 'mistnet_rho', 'none').
%                       'dr' removes gates whose depolarization ratio, estimated
%                       from ZDR and rhoHV (Kilambi et al. 2018), is below
%                       rainDrThreshold (rain is round and uniform: low DR).
%                       'mistnet_rho' removes MistNet rain and gates whose
%                       median-smoothed rhoHV exceeds rainRhoThreshold. Default: 'rho'.
%     'rainRhoThreshold' - rhoHV threshold for rain removal ('rho', 'mistnet_rho'). Default: 0.95.
%     'rainDrThreshold' - Depolarization ratio (dB) below which a gate is rain ('dr'). Default: -12.
%     'rainDilation'  - Grow the rain mask by this distance (meters) to remove
%                       rain edges. Default: 0.
%     'mistnetScreen' - Run MistNet only on volumes with rain-like echo (see
%                       rain_like_area); others get no rain mask. Default: false.
%     'removeOutliers' - MAD value for outlier removal in reflectivity units. Default: 3.
%     'thr_dens'      - Max bird density (bird/km^3) accepted. Default: 1e9.
%     'shape'         - Apply shape correction (true/false). Default: true.
%     'smoothing'     - Smoothing window (meters) for bird density. Default: 2000.
%     'limitElev'     - Max elevation scan to consider (degrees). Default: 5.
%     'minRange'      - Minimum range from radar (meters). Default: 5000.
%     'maxRange'      - Maximum range from radar (meters). Default: Inf.
%     'maxAlt'        - Maximum altitude to consider (meters). Default: 3000.
%     'noEchoRange'   - Gates without echo within this range (meters) are kept
%                       as zero density instead of missing. Default: 0 (off).
%     'data_res'      - Data resolution [range, azimuth*range] in meters. Default: [2000 2000].
%     'waitbar'       - Show waitbar (true/false). Default: true.
%     'plot'          - Plot options (cell array of strings). Default: "".
%     'exportDir'     - Output folder for exported figures. Default: "export/read_rd".
%     'parallel'      - Clean volumes in parallel (no plot/waitbar). Default: true.
%     'maxWorkers'    - Maximum volumes cleaned at the same time. Default: 6.
%     'filters'       - Cleaning steps, applied in this order: names of
%                       functions in bmm4d.radar.filter, each taking and
%                       returning the volume (see default_filters). Default:
%                       built from the options above, which reproduces the
%                       cleaning before the split into filters:
%                       rain_<method>, rain_dilation, beam_blockage,
%                       max_density, outlier_mad, outlier_continuity,
%                       no_echo, to_density, shape, smooth_median,
%                       range_limits, thin. Filters before to_density work
%                       on reflectivity (dBZ), the others on density.
%
%   Each volume also gets FILTER_LOG, a table of the net number of gates
%   each filter made missing; it is negative for filters that fill gates
%   (no_echo sets them to zero, smooth_median fills gaps).
%
%   Example:
%     [radar_scans, shape_fit] = read_rd({'scan1.mat','scan2.mat'}, 'removeRain', 'rho', 'plot', {'rain','export'}, 'exportDir', fullfile("03_simulator_fit", "figures", "read_rd"));

p = inputParser;
addParameter(p,'removeRain','rho'); % Method to remove rain ('rho', 'mistnet', 'none')
addParameter(p,'rainRhoThreshold', 0.95); % rhoHV threshold for rain removal
addParameter(p,'rainDrThreshold', -12); % Depolarization ratio (dB) threshold for rain removal
addParameter(p,'rainDilation', 0); % Grow the rain mask by this distance (m)
addParameter(p,'mistnetScreen', false); % Run MistNet only on volumes with rain-like echo
addParameter(p,'removeOutliers', 3); % MAD value for outlier removal in reflectivity units
addParameter(p,'thr_dens', 1e9); % Max density (bird/km^3) accepted
addParameter(p,'shape', true); % Apply shape correction
addParameter(p,'smoothing', 2000); % Smoothing window (meters)
addParameter(p,'limitElev', 5); % Max elevation scan to consider (degrees)
addParameter(p,'minRange', 5000); % Minimum range from radar (meters)
addParameter(p,'maxRange', Inf); % Maximum range from radar (meters)
addParameter(p,'maxAlt', 3000); % Maximum altitude to consider (meters)
addParameter(p,'noEchoRange', 0); % Range (m) within which no-echo gates are zero density
addParameter(p,'data_res', [2000 2000]); % Data resolution [range, azimuth*range] in meters
addParameter(p,'waitbar', false); % Show waitbar
addParameter(p,'plot', ""); % Plot options
addParameter(p,'exportDir', fullfile("export", "read_rd")); % Output folder for exported figures
addParameter(p,'parallel', true); % Clean volumes in parallel when no plot or waitbar is requested
addParameter(p,'maxWorkers', 6); % Maximum number of volumes cleaned at the same time
addParameter(p,'filters', strings(0)); % Cleaning steps in order (default_filters when empty)

parse(p,varargin{:});
parms = p.Results;
plot_opts = string(parms.plot);
export_dir = parms.exportDir;
if isempty(parms.filters)
    parms.filters = default_filters(parms);
end
parms.filters = string(parms.filters);
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



function filters = default_filters(parms)
%DEFAULT_FILTERS Cleaning steps implied by the read_rd options.
filters = strings(0);
switch parms.removeRain
    case "rho"
        filters(end + 1) = "rain_rho";
    case "dr"
        filters(end + 1) = "rain_dr";
    case {"mistnet", "mistnet_rho"}
        filters(end + 1) = "rain_mistnet";
    case "none"
    otherwise
        warning("Unknown rain removal method: %s", parms.removeRain)
end
if parms.rainDilation > 0
    filters(end + 1) = "rain_dilation";
end
filters(end + 1) = "beam_blockage";
if parms.removeOutliers > 0
    filters = [filters "max_density" "outlier_mad" "outlier_continuity"];
end
filters = [filters "no_echo" "to_density"];
if parms.shape
    filters(end + 1) = "shape";
end
if parms.smoothing > 0
    filters(end + 1) = "smooth_median";
end
filters = [filters "range_limits" "thin"];
end


function [radar_scan, scan_shape_fit] = clean_scan(scan_file, parms, plot_opts, export_dir)
%CLEAN_SCAN Read one radar volume and apply the cleaning filters in order.

% Radar sites are loaded once per MATLAB process (client or worker).
persistent radar_sites
if isempty(radar_sites)
    radar_sites = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
end

% Settings and plot options seen by every filter
[~, filename] = fileparts(scan_file);
parms.filename = string(filename);
parms.plot_opts = plot_opts;
parms.export_dir = export_dir;

radar_scan = bmm4d.radar.rsl2mat_wrap(scan_file, struct('max_elev', parms.limitElev));

% The file header height is the ground elevation of the site, but beam
% heights (slant2ground) are relative to the antenna. Use the antenna
% height (ground + tower) from radar.csv, as for the beam-blockage maps.
site_id = find(radar_sites.name == string(radar_scan.station), 1);
if isempty(site_id)
    warning("read_rd:UnknownStation", "%s not in radar.csv; using the file height.", radar_scan.station)
else
    radar_scan.height = radar_sites.height(site_id);
end

% Unique sweeps below the elevation limit
radar_scan.dz.sweeps = unique_elev_sweeps(radar_scan, 'dz');
radar_scan.rh.sweeps = unique_elev_sweeps(radar_scan, 'rh');
radar_scan.dz.sweeps = radar_scan.dz.sweeps([radar_scan.dz.sweeps.elev]< parms.limitElev);
radar_scan.rh.sweeps = radar_scan.rh.sweeps([radar_scan.rh.sweeps.elev]< parms.limitElev);

% Gates of each sweep on a regular polar grid (sweep2mat sorts azimuths).
% Filters update data (dBZ, then density after to_density), remember the
% gates without echo in the raw sweep (no_echo) and the rain mask.
gates = struct("elev", {}, "vert_half_bw", {}, "range", {}, "az", {}, "data", {}, "unit", {}, ...
    "no_echo", {}, "rain_mask", {});
for i_elev = 1:numel(radar_scan.dz.sweeps)
    [sweep_dbz, sweep_range, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));
    gates(i_elev) = struct("elev", radar_scan.dz.sweeps(i_elev).elev, ...
        "vert_half_bw", radar_scan.dz.sweeps(i_elev).vert_half_bw, "range", sweep_range, "az", sweep_az, ...
        "data", sweep_dbz, "unit", "dBZ", "no_echo", isnan(sweep_dbz), "rain_mask", false(size(sweep_dbz)));
end
radar_scan.gates = gates;

% Apply the filters, counting the gates each one removes
n_removed = zeros(numel(parms.filters), 1);
for i_filter = 1:numel(parms.filters)
    missing_before = arrayfun(@(gate) nnz(isnan(gate.data)), radar_scan.gates);
    filter_fn = str2func("bmm4d.radar.filter." + parms.filters(i_filter));
    radar_scan = filter_fn(radar_scan, parms);
    n_removed(i_filter) = sum(arrayfun(@(gate) nnz(isnan(gate.data)), radar_scan.gates) - missing_before);
end
radar_scan.filter_log = table(parms.filters(:), n_removed, VariableNames=["filter" "n_removed"]);

if ismember("rain", plot_opts)
    [rain_dbz, rain_range, rain_az] = sweep2mat(radar_scan.dz.sweeps(1));
    rain_dbz_filtered = rain_dbz;
    rain_dbz_filtered(radar_scan.gates(1).rain_mask) = nan;
    rain_range_limit = min(parms.maxRange, max(rain_range));

    fig = bmm4d.plot.create_export_figure(1200, 450);
    tiledlayout(1, 3, "TileSpacing", "tight", "Padding", "tight")
    nexttile
    bmm4d.plot.plot_sweep(rain_dbz, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Reflectivity before rain removal", show_colorbar=false);
    nexttile
    bmm4d.plot.plot_sweep(radar_scan.gates(1).rain_mask, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Rain mask: " + string(parms.removeRain), show_colorbar=false);
    nexttile
    bmm4d.plot.plot_sweep(rain_dbz_filtered, rain_az, rain_range, coordinate_system="polar", ...
        max_range=rain_range_limit, ttl="Reflectivity after rain removal", show_colorbar=false);
    bmm4d.radar.filter.finish_plot(fig, parms, "_rain.png")
end

scan_shape_fit = [];
if isfield(radar_scan, "shape_fit")
    scan_shape_fit = radar_scan.shape_fit;
    radar_scan = rmfield(radar_scan, "shape_fit");
end

% Bird density sweeps, with the sweep metadata kept for e.g. sweep2mat
radar_scan.dens = struct();
radar_scan.dens.type = "Bird density";
radar_scan.dens.sweeps = struct([]);
for i_elev=1:numel(radar_scan.gates)
    radar_scan.dens.sweeps(i_elev).data = radar_scan.gates(i_elev).data;
    radar_scan.dens.sweeps(i_elev).az = radar_scan.gates(i_elev).az;
    radar_scan.dens.sweeps(i_elev).range = radar_scan.gates(i_elev).range;
    radar_scan.dens.sweeps(i_elev).elev = radar_scan.dz.sweeps(i_elev).elev;
    radar_scan.dens.sweeps(i_elev).range_bin1 = radar_scan.dz.sweeps(i_elev).range_bin1;
    radar_scan.dens.sweeps(i_elev).gate_size = radar_scan.dz.sweeps(i_elev).gate_size;
    radar_scan.dens.sweeps(i_elev).nbins = radar_scan.dz.sweeps(i_elev).nbins;
    % azim_v follows the azimuth order of sweep2mat
    radar_scan.dens.sweeps(i_elev).azim_v = radar_scan.gates(i_elev).az;
    radar_scan.dens.sweeps(i_elev).prf = radar_scan.dz.sweeps(i_elev).prf;
end
radar_scan = rmfield(radar_scan, intersect(fieldnames(radar_scan), ...
    ["dz" "gates" "sw" "dr" "ph" "spulse" "lpulse" "constants" "params"]));
end
