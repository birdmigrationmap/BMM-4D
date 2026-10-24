%% Scan cleaning study: screen random migration volumes
% Draw random radar volumes from spring and fall migration nights, clean
% each one with every candidate rain-removal method, and record how much
% signal each method removes. Volumes where the methods disagree are the
% ones worth reviewing by eye (SC_B, SC_C).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "cleaning_study", "functions"))
data_dir = fullfile(project_root, "prod", "cleaning_study", "data");
if ~exist(data_dir, "dir")
    mkdir(data_dir)
end
script_tic = tic;

%% Sampling controls ------------------------------------------------------
% Radars east of the Rockies, as in us_night. Nights are drawn uniformly
% from the migration seasons; the time is drawn within +/-4 h of local
% solar midnight, which covers most of the nocturnal flight.
rng(20261022)
n_volume = 300;
years = 2022:2025;
season_windows = [ ...
    "04-10" "05-31"; ...
    "08-25" "10-31"];
night_offset_h = [-4 4];
radar_bounds = [-104 -66 24 50];

%% Cleaning controls ------------------------------------------------------
% Same read_rd options as us_night; only removeRain changes between methods.
method_names = ["none" "rho" "mistnet"];
read_options = { ...
    "limitElev", 7, ...
    "minRange", 3000, ...
    "maxRange", 160000, ...
    "maxAlt", 3500, ...
    "removeOutliers", 3, ...
    "shape", true, ...
    "smoothing", 2000, ...
    "data_res", [2000 2000], ...
    "waitbar", false};
z_edges = 0:100:3000;

%% Draw the sample --------------------------------------------------------
radar = select_radars_in_bbox(fullfile(project_root, "data", "reference", "radar.csv"), radar_bounds);
radar_id = randi(height(radar), n_volume, 1);
season_id = randi(size(season_windows, 1), n_volume, 1);
year_id = years(randi(numel(years), n_volume, 1))';
night_date = NaT(n_volume, 1);
for i_v = 1:n_volume
    window = datetime(year_id(i_v) + "-" + season_windows(season_id(i_v), :));
    night_date(i_v) = window(1) + days(randi(days(diff(window)) + 1) - 1);
end
solar_midnight = night_date + hours(24 - radar.lon(radar_id) / 15);
scan_time = solar_midnight + hours(night_offset_h(1) + diff(night_offset_h) .* rand(n_volume, 1));

screening = table(radar.name(radar_id), radar.lat(radar_id), radar.lon(radar_id), night_date, scan_time, ...
    VariableNames=["radar" "lat" "lon" "night_date" "scan_time"]);
screening.season = repmat("spring", n_volume, 1);
screening.season(season_id == 2) = "fall";

%% Download volumes to the local cache ------------------------------------
screening.file = strings(n_volume, 1);
for i_v = 1:n_volume
    screening.file(i_v) = load_rd(screening.radar(i_v), screening.scan_time(i_v), download=true);
end
fprintf("[SC/A] %d/%d volumes available (%.1f min)\n", nnz(screening.file ~= ""), n_volume, toc(script_tic) / 60)

%% Clean each volume with every method ------------------------------------
volume_files = screening.file;
volume_metrics = cell(n_volume, 1);
volume_error = strings(n_volume, 1);
parfor i_v = 1:n_volume
    if volume_files(i_v) == ""
        volume_error(i_v) = "no file";
        continue
    end
    try
        scans = cell(1, numel(method_names));
        for i_m = 1:numel(method_names)
            scan_cell = read_rd(volume_files(i_v), "removeRain", method_names(i_m), read_options{:});
            scans{i_m} = scan_cell{1};
        end
        volume_metrics{i_v} = compare_cleaning_methods(scans, method_names, z_edges);
    catch err
        volume_error(i_v) = string(err.message);
    end
end

%% Assemble and save ------------------------------------------------------
screening.error = volume_error;
ok = cellfun(@isstruct, volume_metrics);
metric_table = struct2table([volume_metrics{ok}]);
for name = string(metric_table.Properties.VariableNames)
    screening.(name) = nan(n_volume, 1);
    screening.(name)(ok) = metric_table.(name);
end
screening.vid_rho_mistnet_ratio = screening.vid_rho ./ screening.vid_mistnet;

screening_controls = struct( ...
    "seed", 20261022, ...
    "n_volume", n_volume, ...
    "years", years, ...
    "season_windows", season_windows, ...
    "night_offset_h", night_offset_h, ...
    "radar_bounds", radar_bounds, ...
    "method_names", method_names, ...
    "read_options", {read_options}, ...
    "z_edges", z_edges);
save(fullfile(data_dir, "screening.mat"), "screening", "screening_controls")
writetable(screening, fullfile(data_dir, "screening.csv"))

fprintf("[SC/A] %d volumes screened, %d failed\n", nnz(ok), nnz(~ok))
disp(groupsummary(screening(ok, :), [], "median", ["removed_rho" "removed_mistnet" "disagreement"]))
fprintf("[SC/A] done in %.1f min\n", toc(script_tic) / 60)
