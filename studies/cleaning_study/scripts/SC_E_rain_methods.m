%% Scan cleaning study: rhoHV and depolarization ratio against MistNet
% Cleans the screening volumes (SC_A) with every gate-level rain rule and
% measures how close each one is to MistNet, taken as the reference: the
% signal it removes that MistNet keeps, the rain MistNet removes that it
% keeps, and the change in vertically integrated density.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "cleaning_study", "functions"))
data_dir = fullfile(project_root, "studies", "cleaning_study", "data");
load(fullfile(data_dir, "screening.mat"), "screening")
volumes = screening(screening.error == "" & isfile(screening.file), :);
% Rough estimate on a subset: 30 volumes where MistNet removes more than 5%
% of the signal (rain) and 30 others
rng(20261023)
rainy = find(volumes.removed_mistnet > 0.05);
other = find(volumes.removed_mistnet <= 0.05);
volumes = volumes([rainy(randperm(numel(rainy), 30)); other(randperm(numel(other), 30))], :);

%% Methods ----------------------------------------------------------------
% Same read_rd options for all (SC_A settings, 150 km, no rain dilation, no
% no-echo zeros), only the rain rule changes. MistNet is the reference.
method_names = ["none" "mistnet" "rho95" "dr12" "dr15"];
method_options = { ...
    {"removeRain", "none"}, ...
    {"removeRain", "mistnet"}, ...
    {"removeRain", "rho", "rainRhoThreshold", 0.95}, ...
    {"removeRain", "dr", "rainDrThreshold", -12}, ...
    {"removeRain", "dr", "rainDrThreshold", -15}};
reference_id = 2;
read_options = {"limitElev", 7, "minRange", 3000, "maxRange", 150000, "maxAlt", 3500, ...
    "removeOutliers", 3, "shape", true, "smoothing", 2000, "data_res", [2000 2000], ...
    "waitbar", false, "parallel", false};
z_edges = 0:100:3000;

%% Clean every volume with every method ------------------------------------
if isempty(gcp("nocreate"))
    parpool("Processes", 5);   % MistNet needs ~3 GB per worker
end
volume_files = volumes.file;
volume_metrics = cell(height(volumes), 1);
volume_error = strings(height(volumes), 1);
script_tic = tic;
parfor i_v = 1:height(volumes)
    try
        scans = cell(1, numel(method_names));
        for i_m = 1:numel(method_names)
            scan = bmm4d.radar.read_rd(volume_files(i_v), read_options{:}, method_options{i_m}{:});
            scans{i_m} = scan{1};
        end
        volume_metrics{i_v} = compare_to_reference(scans, method_names, reference_id, z_edges);
    catch err
        volume_error(i_v) = string(err.message);
    end
end
delete(gcp("nocreate"))
fprintf("[SC/E] %d volumes in %.1f min, %d failed\n", height(volumes), toc(script_tic) / 60, nnz(volume_error ~= ""))

%% Collect and save ---------------------------------------------------------
ok = volume_error == "";
rain_methods = [volumes(ok, ["radar" "scan_time" "season" "file"]), struct2table([volume_metrics{ok}])];
save(fullfile(data_dir, "rain_methods.mat"), "rain_methods", "method_names", "read_options")
writetable(rain_methods, fullfile(data_dir, "rain_methods.csv"))
