%% Production test 3: build the observation system
% Build the shared grid, then each radar's forward-operator rows from its
% cleaned volumes (radars in parallel), and stack them into one compact
% observation system with per-radar velocity profiles.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"))
% Case: "keywest", "midwest" or "us_east"; suffix _rho or _none for other rain
% removal. A case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
n_workers = 6;
load(fullfile(c.data_dir, "manifest.mat"), "volume_files")
n_radar = numel(c.radars);
n_time = numel(c.times);

%% Build the shared grid ------------------------------------------------------
% Only the kernel x/y/z spacing is needed, so one kernel file is enough.
stage_tic = tic;
kernel_files = dir(fullfile(c.kernel_dir, "F_*.mat"));
kernel_ref = load(fullfile(kernel_files(1).folder, kernel_files(1).name), "x", "y", "z");
radar = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
[~, radar_id] = ismember(c.radars, radar.name);
radar_meta = cell(n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        radar_meta{i_r, i_t} = struct("lon", radar.lon(radar_id(i_r)), "lat", radar.lat(radar_id(i_r)), "datetime", c.times(i_t));
    end
end
g = bmm4d.grid.create_grid({struct("x", kernel_ref.x, "y", kernel_ref.y, "z", kernel_ref.z)}, radar_meta, ...
    c.grid_bounds(1), c.grid_bounds(2), c.grid_bounds(3), c.grid_bounds(4), c.heightmax);
if ~isempty(c.water_override)
    [lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
    g.mask_water(inpolygon(lon_grid, lat_grid, c.water_override(1, :), c.water_override(2, :))) = true;
end
clear radar_meta
% Solve grid: blocks of c.solve_aggregate fine cells (forward-operator
% columns are summed per block in construct_F)
g_fine = g;
g = bmm4d.grid.coarsen_grid(g_fine, c.solve_aggregate);
stage_log(c, "P3_grid", stage_tic, sprintf("grid %dx%dx%d (fine %dx%dx%d)", g.sz, g_fine.sz))

%% Build each radar's rows in parallel ----------------------------------------
stage_tic = tic;
pool = gcp("nocreate");
if isempty(pool) || pool.NumWorkers ~= n_workers
    delete(pool);
    parpool("Processes", n_workers);
end
rows = cell(n_radar, 1);
kernel_dir = c.kernel_dir;
min_row_mass = c.min_row_mass;
aggregate = c.solve_aggregate;
% Static clutter masks from P2b, if computed
clutter_masks = cell(n_radar, 1);
clutter_file = fullfile(c.data_dir, "static_clutter.mat");
if isfile(clutter_file)
    loaded = load(clutter_file, "clutter_masks", "clutter_radars");
    [found, id] = ismember(c.radars, loaded.clutter_radars);
    clutter_masks(found) = loaded.clutter_masks(id(found));
end
clutter_max_elev = c.static_clutter.max_elev;
near_range_m = c.near_range_m;
no_echo_max_range_m = c.no_echo_max_range_m;
parfor i_r = 1:n_radar
    try
        rows{i_r} = build_radar_rows(volume_files(i_r, :), n_time, g_fine, kernel_dir, min_row_mass, aggregate, ...
            clutter_masks{i_r}, clutter_max_elev, near_range_m, no_echo_max_range_m);
    catch err
        warning("P3:RadarFailed", "%s: %s", c.radars(i_r), err.message)
    end
end
delete(gcp("nocreate"))   % release worker memory before stacking and the inversion
clear g_fine
kept = ~cellfun(@isempty, rows);
rows = rows(kept);
radar_selected = radar(radar_id(kept), :);
stage_log(c, "P3_radar_rows", stage_tic, sprintf("%d/%d radars, %d sweeps without kernel", ...
    nnz(kept), n_radar, sum(cellfun(@(r) r.n_dropped_sweeps, rows))))

%% Stack the observation system ------------------------------------------------
% Rows are stacked radar by radar; row_radar_id gives each row's radar.
stage_tic = tic;
n_radar = numel(rows);
radar_scans = cell(n_radar, n_time);
scan_available = false(n_radar, n_time);
for i_r = 1:n_radar
    scan_available(i_r, rows{i_r}.time_id) = true;
    for i_t = 1:n_time
        radar_scans{i_r, i_t} = struct( ...
            "station", rows{i_r}.site.station, "lat", rows{i_r}.site.lat, "lon", rows{i_r}.site.lon, ...
            "height", rows{i_r}.site.height, "datetime", c.times(i_t), "uv", rows{i_r}.uv(:, :, i_t));
    end
end
F_parts = cellfun(@(r) r.F, rows, UniformOutput=false);
value_parts = cellfun(@(r) r.values, rows, UniformOutput=false);
rows_range = cellfun(@(r) r.row_range_m, rows, UniformOutput=false);
rows_elev = cellfun(@(r) r.row_elev_deg, rows, UniformOutput=false);
observations = struct( ...
    "F", vertcat(F_parts{:}), ...
    "values", vertcat(value_parts{:}), ...
    "row_radar_id", reshape(repelem(uint16(1:n_radar), cellfun(@(r) size(r.F, 1), rows)), [], 1), ...
    "row_range_m", vertcat(rows_range{:}), ...
    "row_elev_deg", vertcat(rows_elev{:}));
clear F_parts value_parts rows_range rows_elev
clear rows
case_controls = c;
save(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected", ...
    "scan_available", "case_controls", "-v7.3")
stage_log(c, "P3_save_case", stage_tic, sprintf("%.3g rows, nnz(F) %.3g", size(observations.F, 1), nnz(observations.F)))
