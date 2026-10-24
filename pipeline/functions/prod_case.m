function c = prod_case(case_name)
%PROD_CASE Settings of one production test case.
%   C = PROD_CASE(CASE_NAME) returns the radars, times, grid, cleaning and
%   inversion settings shared by all pipeline/scripts stages. Edit this file to
%   try other settings; every stage records the settings it used.
%
%   Cases (append _rho, _dr, _mistnet or _none to choose the rain removal, and
%   _coarse for a 4 km x 200 m solve grid):
%   - "keywest": KBYX, 23 April 2018, 01:00-07:00 UTC every 15 min (fast test,
%     same radar/times/weights as the paper case 09).
%   - "midwest": 17 radars, 16 May 2020, 01:30-12:00 UTC every 30 min (the
%     paper case 11; scaling test for a network).
%   - "us_east": radars east of -104 deg, one night (default 21-22 Oct 2026).

% Settings shared by all cases ---------------------------------------------
c = struct();
c.name = string(case_name);
c.r_lim = 150000;                 % kernel set (F_elev/150km)
c.heightmax = 3000;               % grid top above sea level (m)
c.min_row_mass = 0.75;            % beam weight required inside the grid
c.flight_layer_agl_m = 3000;      % solved cells: ground to this height above it
c.solve_aggregate = [1 1 1];      % solve grid: blocks of [x y z] kernel-grid cells
c.active_range_m = 250000;        % solved columns: within this distance of a radar
                                  % (> half the radar spacing, so gaps are filled)
c.closest_accepted_min = 10;      % max offset between requested and scan time
c.sun_limit_deg = -6;             % night: sun below civil twilight

% Cleaning: MistNet removes precipitation; other read_rd options as paper
% except maxRange: MistNet's Cartesian grid ends at 150 km, so rain beyond
% it would stay unclassified (seen at 150-160 km east of KBYX). Unlike the
% paper, gates without echo are zero-density observations (not missing),
% so beams above the bird layer constrain it (cleaning_study SC_D).
c.read_options = { ...
    "removeRain", "mistnet", ...
    "limitElev", 7, ...
    "minRange", 3000, ...
    "maxRange", 150000, ...       % MistNet only classifies within 150 km
    "maxAlt", 3500, ...
    "noEchoRange", 150000, ...    % no-echo gates are observations of zero density
    "rainDilation", 2000, ...     % grow MistNet rain by 2 km (blocky rain edges)
    "mistnetScreen", true, ...    % skip MistNet on volumes without rain-like echo
    "removeOutliers", 3, ...
    "shape", true, ...
    "smoothing", 2000, ...
    "data_res", [0 0], ...        % full gate resolution; clean_volume thins after filtering
    "waitbar", false, ...
    "parallel", false};           % volumes are parallelized one level up
c.vp = struct("zstep_m", 100, "rmin_m", 5000, "rmax_m", 37500, "zmax_m", 3000);

% Filters applied after read_rd (clean_volume), on full-resolution sweeps:
% - max_density: implausible bird density (birds km^-3); 3000 ~ 40 dBZ.
% - min_patch_km2: echo patches smaller than this are removed (point
%   targets, rain fragments, speckle); 0 disables.
% - data_res: then thin gates to this spacing [range tangential] (m).
c.filter = struct("max_density", 3000, "min_patch_km2", 25, "data_res", [1000 1000]);

% Static clutter (P2b): gates > threshold_db above their window_m local
% median in more than min_fraction of the night's volumes (wind farms,
% towers) are removed from sweeps below max_elev in P3.
% Near-range gates are dropped in P3: bird echo is suppressed close to the
% radar (clutter filtering; lowest sweep at 27% of the 30-50 km density at
% 5-10 km). 25 km also removes the dark centre seen with 15 km (radar
% distance bias, P4b).
c.near_range_m = 25000;
% No-echo gates (zero density) are kept only within this range in P3: far
% out the beam is mostly above the birds and its weak lower edge falls
% below the radar sensitivity, so "no echo" is not "no birds" there.
% 75 km halves the bright ring at 25-50 km and flattens 75-150 km (P4b).
c.no_echo_max_range_m = 75000;
c.static_clutter = struct("threshold_db", 8, "window_m", 20000, "min_fraction", 0.5, ...
    "max_range_m", 150000, "max_elev", 2);

% Inversion (paper Midwest weights by default).
c.inversion = struct( ...
    "lambda_profile", 10, ...
    "lambda_time_profile", 1, ...
    "density_prior_idw_neighbors", 3, ...
    "density_prior_idw_power", 2.25, ...
    "prior_interpolation", "gaussian", ...  % smooth radar weights (paper: "idw", 3 nearest radars)
    "prior_length_scale_m", 150000, ...     % Gaussian length scale of the prior
    "prior_include_water", true, ...        % fit profiles on water cells too
    "velocity_idw_neighbors", 10, ...
    "velocity_idw_power", 1.6, ...
    "prior", 0.3, ...             % 0.1 in the paper Midwest case; 0.3 halves negative cells and iterations (P4b)
    "smoothness", 0.3, ...
    "advection", 1, ...
    "coeff_horver", 0, ...
    "shoreline_crossing_coeff", 1, ...
    "maxit", 500, ...
    "tol", 1e-4, ...              % approximate solve (paper: 1e-6)
    "advection_relax_exponent", 2, ...  % weaken advection where the prior changes (takeoff, landing)
    "nonnegative", true, ...      % projected FISTA after CG: densities >= 0
    "projected_maxit", 100);
c.water_override = [];            % [lon; lat] polygon forced to water
c.domain = [];                    % geographic solve domain (build_domain_mask); [] = radar ranges

% A suffix selects the rain removal, e.g. "keywest_rho" or "midwest_none";
% without suffix MistNet is used.
name_parts = split(c.name, "_");
base_name = name_parts(1);
if base_name == "us"
    base_name = "us_east";
end
rain_method = "mistnet";
rain_part = name_parts(ismember(name_parts, ["rho" "dr" "mistnet" "none"]));
if ~isempty(rain_part)
    rain_method = rain_part(1);
end
c.read_options = set_option(c.read_options, "removeRain", rain_method);
% Suffix "coarse" solves on 4 km x 200 m blocks, e.g. "midwest_coarse".
use_coarse = any(name_parts == "coarse");

% Case-specific settings -------------------------------------------------
switch base_name
    case "keywest"
        c.radars = "KBYX";
        % Paper case 09 cleaning options except the rain method.
        c.read_options = set_option(c.read_options, "removeOutliers", 2);
        c.times = datetime("2018-04-23 01:00:00"):minutes(15):datetime("2018-04-23 07:00:00");
        c.closest_accepted_min = 15;
        c.grid_bounds = [-83.2, -80.2, 22.7, 25.6];
        c.figure_bounds = [-83.2, -80.0, 22.8, 25.7];
        c.water_override = [-82.05 -81.98 -81.93 -81.97 -82.04; 24.02 24.07 24.02 23.97 23.98];
        % Published Key West weights (paper case 09).
        c.inversion.velocity_idw_neighbors = 1;
        c.inversion.prior = 0.003;
        c.inversion.smoothness = 1;
        c.inversion.advection = 3;
        c.inversion.coeff_horver = -0.5;
        c.inversion.shoreline_crossing_coeff = 0.3;
        c.inversion.maxit = 300;
    case "midwest"
        c.grid_bounds = [-97, -85, 38, 45];
        c.figure_bounds = [-97, -86, 38.5, 45];
        c.radars = select_radars_in_bbox(data_path("reference", "radar.csv"), c.grid_bounds).name';
        c.times = datetime("2020-05-16") + (1.5:0.5:12) ./ 24;
    case "us_east"
        night_date = datetime("2026-10-21");
        radar_bounds = [-104 -66 24 50];
        c.radars = select_radars_in_bbox(data_path("reference", "radar.csv"), radar_bounds).name';
        c.times = (night_date + hours(23.5)):minutes(30):(night_date + hours(36));
        radar = select_radars_in_bbox(data_path("reference", "radar.csv"), radar_bounds);
        pad_lat = 160 / 111.32;
        pad_lon = 160 / (111.32 * cosd(max(radar.lat)));
        c.grid_bounds = [min(radar.lon) - pad_lon, max(radar.lon) + pad_lon, min(radar.lat) - pad_lat, max(radar.lat) + pad_lat];
        c.figure_bounds = c.grid_bounds;
        c.filter.data_res = [2000 2000];
        % 4 km x 100 m solve grid (4x fewer unknowns than 2 km x 100 m)
        c.solve_aggregate = [2 2 1];
        % Fixed geographic domain: US east of the High Plains and 50 km offshore
        % (coasts, Gulf of Mexico and Great Lakes alike)
        c.domain = struct("west_lon", -104, "offshore_km", 50);
    case "sandbox"
        % Fast prototyping case (studies/sandbox): Midwest part of the us_east
        % night, 6 frames at the start of the night, same cleaning (shared
        % volume cache), no data cuts in P3 (decided at run time).
        radar_bounds = [-100 -85 36 46];
        c.radars = select_radars_in_bbox(data_path("reference", "radar.csv"), radar_bounds).name';
        c.times = datetime("2026-10-22 00:30"):minutes(30):datetime("2026-10-22 03:00");
        c.grid_bounds = radar_bounds;
        c.figure_bounds = radar_bounds;
        c.filter.data_res = [2000 2000];
        c.solve_aggregate = [2 2 1];
        c.domain = struct("west_lon", -180, "offshore_km", 50);
        c.near_range_m = 0;
        c.no_echo_max_range_m = Inf;
    otherwise
        error("prod_case:UnknownCase", "Unknown case %s.", case_name)
end

if use_coarse
    c.solve_aggregate = [2 2 2];
end
% Suffix "nearNN" drops gates within NN km (default 15), e.g. "us_east_near25".
near_part = regexp(c.name, "near(\d+)", "tokens", "once");
if ~isempty(near_part)
    c.near_range_m = str2double(near_part{1}) * 1000;
end
% Suffix "zeroNN" keeps no-echo gates only within NN km, e.g. "us_east_zero75".
zero_part = regexp(c.name, "zero(\d+)", "tokens", "once");
if ~isempty(zero_part)
    c.no_echo_max_range_m = str2double(zero_part{1}) * 1000;
end

% Night window: keep frames where the sun is below sun_limit_deg at one
% radar at least, and use each radar only while it is night there
% (radar_night, [radar time]). Daytime echo is insects or landed birds,
% and takeoff/landing are not migration.
radar_table = readtable(data_path("reference", "radar.csv"), TextType="string");
[~, radar_id] = ismember(c.radars, radar_table.name);
c.radar_night = false(numel(c.radars), numel(c.times));
for i_t = 1:numel(c.times)
    c.radar_night(:, i_t) = sun_elevation(radar_table.lat(radar_id), radar_table.lon(radar_id), c.times(i_t)) < c.sun_limit_deg;
end
night_frame = any(c.radar_night, 1);
c.times = c.times(night_frame);
c.radar_night = c.radar_night(:, night_frame);

% Folders: products outside the synced project (local_settings.json).
c.data_dir = fullfile(local_data_dir(), "prod", c.name);
c.figure_dir = fullfile(local_folder("pipeline_figure_dir", fullfile("pipeline", "figures")), c.name);
c.kernel_dir = data_path("derived", "F_elev", num2str(c.r_lim / 1000) + "km");
end

function options = set_option(options, name, value)
% Set one name-value pair in a read_rd option cell.
id = find(cellfun(@(x) string(x) == name, options(1:2:end)));
options{2 * id} = value;
end
