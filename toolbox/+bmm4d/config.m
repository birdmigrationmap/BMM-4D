function c = config(name, varargin)
%CONFIG Settings of one reconstruction: radars, times, grid, cleaning, inversion.
%   C = BMM4D.CONFIG(NAME, "radars", R, "times", T, "grid_bounds", B, ...)
%   returns the settings used by every stage (get_volumes, clean_volumes,
%   static_clutter, build_observations, invert, export_maps). NAME names
%   the case and its product folder. Any field below can be given as a
%   name-value pair: a struct value (vp, filter, static_clutter,
%   inversion, domain) only replaces the fields it contains, and a
%   read_options cell only replaces the read_rd options it names.
%
%   Required: radars (string array), times (datetime), grid_bounds
%   ([lon_min lon_max lat_min lat_max]). figure_bounds defaults to
%   grid_bounds.
%
%   Derived at the end:
%   - radar_night [radar time]: the sun is below sun_limit_deg at the
%     radar; frames where no radar is in the night are dropped from times.
%   - data_dir: <local_data_dir>/prod/<name> (products: manifest, case,
%     inversion, stage log); figure_dir: <figure_dir setting>/<name>;
%     kernel_dir: the beam kernels of r_lim in the data folder.

% Grid and domain -----------------------------------------------------------
c = struct();
c.name = string(name);
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

% Static clutter (static_clutter stage): gates > threshold_db above their
% window_m local median in more than min_fraction of the night's volumes
% (wind farms, towers) are removed from sweeps below max_elev.
% Near-range gates are dropped in build_observations: bird echo is
% suppressed close to the radar (clutter filtering; lowest sweep at 27% of
% the 30-50 km density at 5-10 km). 25 km also removes the dark centre
% seen with 15 km (radar distance bias).
c.near_range_m = 25000;
% No-echo gates (zero density) are kept only within this range: far out
% the beam is mostly above the birds and its weak lower edge falls below
% the radar sensitivity, so "no echo" is not "no birds" there. 75 km
% halves the bright ring at 25-50 km and flattens 75-150 km.
c.no_echo_max_range_m = 75000;
c.static_clutter = struct("threshold_db", 8, "window_m", 20000, "min_fraction", 0.5, ...
    "max_range_m", 150000, "max_elev", 2);

% Inversion (paper Midwest weights by default)
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
    "prior", 0.3, ...             % 0.1 in the paper Midwest case; 0.3 halves negative cells and iterations
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
c.radars = strings(0);
c.times = datetime.empty;
c.grid_bounds = [];
c.figure_bounds = [];

% Overrides -----------------------------------------------------------------
for k = 1:2:numel(varargin)
    field = string(varargin{k});
    value = varargin{k + 1};
    if ~isfield(c, field)
        error("bmm4d:config:UnknownSetting", "Unknown setting %s.", field)
    end
    if field == "read_options"
        for i = 1:2:numel(value)
            c.read_options = set_option(c.read_options, value{i}, value{i + 1});
        end
    elseif isstruct(value) && isstruct(c.(field))
        for sub = string(fieldnames(value))'
            c.(field).(sub) = value.(sub);
        end
    else
        c.(field) = value;
    end
end
if isempty(c.radars) || isempty(c.times) || isempty(c.grid_bounds)
    error("bmm4d:config:MissingSetting", "radars, times and grid_bounds are required.")
end
if isempty(c.figure_bounds)
    c.figure_bounds = c.grid_bounds;
end

% Night window ----------------------------------------------------------------
% Keep frames where the sun is below sun_limit_deg at one radar at least,
% and use each radar only while it is night there. Daytime echo is insects
% or landed birds, and takeoff/landing are not migration.
radar_table = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
[~, radar_id] = ismember(c.radars, radar_table.name);
c.radar_night = false(numel(c.radars), numel(c.times));
for i_t = 1:numel(c.times)
    c.radar_night(:, i_t) = bmm4d.util.sun_elevation(radar_table.lat(radar_id), radar_table.lon(radar_id), c.times(i_t)) < c.sun_limit_deg;
end
night_frame = any(c.radar_night, 1);
c.times = c.times(night_frame);
c.radar_night = c.radar_night(:, night_frame);

% Folders: products outside the synced project (local_settings.json)
c.data_dir = fullfile(bmm4d.util.local_data_dir(), "prod", c.name);
c.figure_dir = fullfile(bmm4d.util.local_folder("figure_dir", "figures"), c.name);
c.kernel_dir = bmm4d.util.data_path("derived", "F_elev", num2str(c.r_lim / 1000) + "km");
end

function options = set_option(options, name, value)
% Set one name-value pair in a read_rd option cell (appended if absent).
id = find(cellfun(@(x) string(x) == string(name), options(1:2:end)));
if isempty(id)
    options(end + 1:end + 2) = {string(name), value};
else
    options{2 * id} = value;
end
end
