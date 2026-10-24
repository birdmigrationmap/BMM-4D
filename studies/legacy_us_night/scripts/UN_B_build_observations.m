%% US night reconstruction: clean scans and build observations
% Process the night radar by radar so that a radar's full volumes are only
% in memory while it is processed: clean the volumes, estimate velocity
% profiles, build the forward-operator rows on the shared grid, and keep
% only the compact observation rows. Each radar is cached so the script can
% resume, and radars run in parallel on a bounded pool.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "us_night", "functions"))
script_tic = tic;

night_date = datetime("2026-10-21");
night_tag = string(night_date, "yyyyMMdd");
night_dir = fullfile(local_data_dir(), "us_night", "data", night_tag);
radar_cache_dir = fullfile(night_dir, "radar_cache");
if ~exist(radar_cache_dir, "dir")
    mkdir(radar_cache_dir)
end
load(fullfile(night_dir, "manifest.mat"), "radar_selected", "scan_file", "night_controls")
scan_time_requested = night_controls.scan_time_requested;

%% Processing controls ----------------------------------------------------
r_lim = 150000;
% Heights are above sea level. Ground within radar coverage reaches 1773 m
% at the 99.9th percentile (2156 m max, Black Hills); 3500 m keeps a
% 2 km flight layer above almost all of it. Cells below ground and above
% the flight layer are excluded later by the active mask (UN_C).
heightmax = 3500;
min_time_fraction = 0.5;
min_row_mass = 0.75;  % as for the real-data paper cases (S09-S11)
delete_raw = false;
n_workers = 6;
read_options = { ...
    "removeRain", "rho", ...
    "limitElev", 7, ...
    "minRange", 3000, ...
    "maxRange", 160000, ...
    "maxAlt", 3500, ...
    "removeOutliers", 3, ...
    "shape", true, ...
    "smoothing", 2000, ...
    "data_res", [2000 2000], ...
    "waitbar", false};
kernel_dir = data_path("derived", "F_elev", num2str(r_lim / 1000) + "km");

n_radar = height(radar_selected);
n_time = numel(scan_time_requested);

%% Build the shared grid --------------------------------------------------
% The grid covers every selected radar plus its maximum range. Only the
% kernel x/y/z spacing is needed here, so one kernel file is enough.
kernel_files = dir(fullfile(kernel_dir, "F_*.mat"));
kernel_ref = load(fullfile(kernel_files(1).folder, kernel_files(1).name), "x", "y", "z");
grid_kernels = {struct("x", kernel_ref.x, "y", kernel_ref.y, "z", kernel_ref.z)};

pad_lat_deg = 160 / 111.32;
pad_lon_deg = 160 / (111.32 * cosd(max(radar_selected.lat)));
grid_bounds = [min(radar_selected.lon) - pad_lon_deg, max(radar_selected.lon) + pad_lon_deg, ...
    min(radar_selected.lat) - pad_lat_deg, max(radar_selected.lat) + pad_lat_deg];

% Per-radar caches are only valid for this exact grid.
grid_key = [grid_bounds heightmax];

radar_meta = cell(n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        radar_meta{i_r, i_t} = struct("lon", radar_selected.lon(i_r), "lat", radar_selected.lat(i_r), ...
            "datetime", scan_time_requested(i_t));
    end
end
g = create_grid(grid_kernels, radar_meta, grid_bounds(1), grid_bounds(2), grid_bounds(3), grid_bounds(4), heightmax);
clear radar_meta
fprintf("[UN/B] grid %d x %d x %d (%.3g cells), %d times\n", g.sz, prod(g.sz), n_time)

%% Process radars in parallel ---------------------------------------------
% Each worker cleans one whole radar (build_radar_cache), which keeps all
% workers busy; the serial per-volume work dominates. A worker holds one
% radar's raw volumes, kernels and rows (about 2-3 GB), so the pool is
% limited to n_workers to stay within laptop memory.
% A cache is reused only for the same grid and the same volumes, so a
% re-run after more of the night is downloaded rebuilds that radar.
todo = false(n_radar, 1);
for i_r = 1:n_radar
    cache_file = fullfile(radar_cache_dir, radar_selected.name(i_r) + ".mat");
    available_id = find(scan_file(i_r, :) ~= "");
    todo(i_r) = true;
    if isfile(cache_file)
        cached = load(cache_file, "grid_key", "radar_files");
        todo(i_r) = ~(isfield(cached, "grid_key") && isequal(cached.grid_key, grid_key) ...
            && isequal(cached.radar_files, scan_file(i_r, available_id)));
    end
end
todo_id = find(todo);
fprintf("[UN/B] %d radars to process, %d cached\n", numel(todo_id), n_radar - numel(todo_id))

pool = gcp("nocreate");
if isempty(pool) || pool.NumWorkers ~= n_workers
    delete(pool);
    parpool("Processes", n_workers);
end
status = strings(numel(todo_id), 1);
stations = radar_selected.name;
parfor k = 1:numel(todo_id)
    i_r = todo_id(k);
    available_id = find(scan_file(i_r, :) ~= "");
    status(k) = build_radar_cache(stations(i_r), scan_file(i_r, available_id), available_id, n_time, g, grid_key, ...
        read_options, kernel_dir, min_time_fraction, min_row_mass, fullfile(radar_cache_dir, stations(i_r) + ".mat"));
end
fprintf("[UN/B] %s\n", status)
if delete_raw
    raw_files = scan_file(scan_file ~= "");
    delete(raw_files{:})
end

%% Assemble the compact observation system --------------------------------
cache_list = dir(fullfile(radar_cache_dir, "*.mat"));
station_cached = erase(string({cache_list.name}), ".mat");
radar_kept = ismember(radar_selected.name, station_cached);
radar_selected = radar_selected(radar_kept, :);
n_radar = height(radar_selected);

F_parts = cell(n_radar, 1);
value_parts = cell(n_radar, 1);
row_radar_parts = cell(n_radar, 1);
radar_scans = cell(n_radar, n_time);
scan_available = false(n_radar, n_time);
for i_r = 1:n_radar
    cached = load(fullfile(radar_cache_dir, radar_selected.name(i_r) + ".mat"));
    F_parts{i_r} = cached.radar_F;
    value_parts{i_r} = cached.radar_values;
    row_radar_parts{i_r} = repmat(uint16(i_r), size(cached.radar_F, 1), 1);
    scan_available(i_r, cached.radar_time_id) = true;
    for i_t = 1:n_time
        radar_scans{i_r, i_t} = struct( ...
            "station", cached.radar_site.station, ...
            "lat", cached.radar_site.lat, ...
            "lon", cached.radar_site.lon, ...
            "height", cached.radar_site.height, ...
            "datetime", scan_time_requested(i_t), ...
            "uv", cached.radar_uv(:, :, i_t));
    end
end

% Rows are stacked radar by radar; row_radar_id replaces the global row
% layout of construct_F when fitting the density prior.
observations = struct( ...
    "F", vertcat(F_parts{:}), ...
    "values", vertcat(value_parts{:}), ...
    "row_radar_id", vertcat(row_radar_parts{:}));
clear F_parts value_parts row_radar_parts cached

case_controls = struct( ...
    "night_controls", night_controls, ...
    "r_lim", r_lim, ...
    "heightmax", heightmax, ...
    "grid_bounds", grid_bounds, ...
    "min_time_fraction", min_time_fraction, ...
    "min_row_mass", min_row_mass, ...
    "read_options", {read_options}, ...
    "scan_available", scan_available);

case_file = fullfile(night_dir, "case.mat");
save(case_file, "observations", "g", "radar_scans", "radar_selected", "case_controls", "-v7.3")
fprintf("[UN/B] %d radars kept, %.3g observation rows, nnz(F) = %.3g\n", ...
    n_radar, size(observations.F, 1), nnz(observations.F))
fprintf("[UN/B] saved %s in %.1f min\n", case_file, toc(script_tic) / 60)
