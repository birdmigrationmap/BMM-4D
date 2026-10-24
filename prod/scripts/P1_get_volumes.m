%% Production test 1: locate or download raw radar volumes
% Find the volume closest to each requested radar-time pair in the local raw
% caches and download the missing ones in parallel.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"))
% Case: "keywest", "midwest" or "us_east"; suffix _rho or _none for other rain
% removal. A case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
if ~isfolder(c.data_dir)
    mkdir(c.data_dir)
end

%% Find local files and download missing ones -----------------------------
stage_tic = tic;
raw_files = load_rd(c.radars, c.times, closest_accepted=c.closest_accepted_min);
% Only night volumes are used (prod_case radar_night)
raw_files(~c.radar_night) = "";
need = raw_files == "" & c.radar_night;
n_downloaded = 0;
if any(need(:))
    n_downloaded = download_volumes(c.radars, c.times, need, c.closest_accepted_min, 16);
    raw_files = load_rd(c.radars, c.times, closest_accepted=c.closest_accepted_min);
    raw_files(~c.radar_night) = "";
end
stage_log(c, "P1_get_volumes", stage_tic, sprintf("%d/%d night volumes, %d frames, %d downloaded", ...
    nnz(raw_files ~= ""), nnz(c.radar_night), numel(c.times), n_downloaded))

%% Save the manifest --------------------------------------------------------
radars = c.radars;
times = c.times;
save(fullfile(c.data_dir, "manifest.mat"), "radars", "times", "raw_files")
missing_radar = c.radars(all(raw_files == "", 2));
if ~isempty(missing_radar)
    fprintf("[prod/%s] no data: %s\n", c.name, strjoin(missing_radar, ", "))
end
