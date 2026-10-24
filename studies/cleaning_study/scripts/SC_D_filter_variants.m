%% Scan cleaning study: filter variants on problem volumes and the gallery
% Cleans each review volume with the current prod settings and with
% candidate settings, then exports one comparison figure per volume.
% Volumes: the 24 gallery cases plus the volumes behind the aberrations seen
% in the prod maps (KILX rings, near-radar discs, KFSD streaks, KBYX cells).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"))
addpath(fullfile(project_root, "studies", "cleaning_study", "functions"))
data_dir = fullfile(local_data_dir(), "prod", "cleaning_variants");
figure_dir = fullfile(project_root, "studies", "cleaning_study", "figures", "variants");
if ~isfolder(figure_dir)
    mkdir(figure_dir)
end

%% Volumes ----------------------------------------------------------------
load(fullfile(project_root, "studies", "cleaning_study", "data", "cases.mat"), "cases")
raw_dir = data_path("raw_scans", "nexrad");
problem = table( ...
    ["P01"; "P02"; "P03"; "P04"; "P05"; "P06"], ...
    ["kilx_rings"; "kfsd_streaks"; "kfsd_dawn"; "kbyx_cells"; "kbyx_cells"; "kdmx_night"], ...
    [fullfile(raw_dir, "KILX_20200516", "KILX20200516_013241_V06"); ...
     fullfile(raw_dir, "KFSD_20200516", "KFSD20200516_095726_V06"); ...
     fullfile(raw_dir, "KFSD_20200516", "KFSD20200516_113048_V06"); ...
     fullfile(raw_dir, "KBYX_20180423", "KBYX20180423_010140_V06"); ...
     fullfile(raw_dir, "KBYX_20180423", "KBYX20180423_012927_V06"); ...
     fullfile(raw_dir, "KDMX_20200516", "KDMX20200516_053352_V06")], ...
    VariableNames=["case_id" "group" "file"]);
volumes = [problem; cases(:, ["case_id" "group" "file"])];
volumes = volumes(isfile(volumes.file), :);

%% Variants ---------------------------------------------------------------
% "current": prod defaults before this study (MistNet, no-echo missing).
% "hybrid95" / "hybrid90": MistNet or median rhoHV above 0.95 / 0.90, rain
% mask grown by 2 km, no-echo gates kept as zero density.
c = prod_case("midwest");
current_options = c.read_options;
id = find(cellfun(@(x) string(x) == "noEchoRange", current_options(1:2:end)));
current_options([2*id-1 2*id]) = [];
hybrid_options = @(threshold) [current_options, {"noEchoRange", 150000, ...
    "rainRhoThreshold", threshold, "rainDilation", 2000}];
variant_names = ["current" "hybrid95" "hybrid90"];
variant_options = {current_options, hybrid_options(0.95), hybrid_options(0.90)};
id = find(cellfun(@(x) string(x) == "removeRain", current_options(1:2:end)));
variant_options{2}{2 * id} = "mistnet_rho";
variant_options{3}{2 * id} = "mistnet_rho";

%% Clean every volume with every variant ---------------------------------
[i_vol, i_var] = ndgrid(1:height(volumes), 1:numel(variant_names));
jobs_file = fullfile(data_dir, variant_names(i_var(:))', volumes.case_id(i_vol(:)) + ".mat");
for name = variant_names
    if ~isfolder(fullfile(data_dir, name))
        mkdir(fullfile(data_dir, name))
    end
end
status = strings(numel(jobs_file), 1);
todo = find(~isfile(jobs_file));
if ~isempty(todo) && isempty(gcp("nocreate"))
    parpool("Processes", 4);   % MistNet needs ~3 GB per worker
end
job_file_todo = jobs_file(todo);
job_raw_todo = volumes.file(i_vol(todo));
job_options_todo = variant_options(i_var(todo));
status_todo = strings(numel(todo), 1);
parfor k = 1:numel(todo)
    status_todo(k) = clean_volume(job_raw_todo(k), job_options_todo{k}, c.vp, c.filter, job_file_todo(k));
end
status(todo) = status_todo;
delete(gcp("nocreate"))
fprintf("[SC/D] cleaned %d volumes, %d failed\n", numel(todo), sum(status_todo ~= "ok"))

%% Comparison figures -----------------------------------------------------
for i_v = 1:height(volumes)
    files = jobs_file(i_vol(:) == i_v);
    if ~all(isfile(files))
        continue
    end
    output_file = fullfile(figure_dir, volumes.case_id(i_v) + "_" + volumes.group(i_v) + ".png");
    plot_filter_variants(volumes.file(i_v), files, variant_names, output_file);
end
fprintf("[SC/D] figures in %s\n", figure_dir)

