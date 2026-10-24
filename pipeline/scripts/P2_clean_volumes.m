%% Production test 2: clean volumes into the shared volume cache
% Clean every volume once (MistNet by default, see prod_case) and keep a
% compact copy with its velocity profile. Volumes are independent, so all
% radars and times run in parallel on a bounded pool; cleaned volumes are
% reused by every case and stage that uses the same cleaning settings.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"))
% Case: "keywest", "midwest" or "us_east"; suffix _rho or _none for other rain
% removal. A case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
% Measured worker memory (Key West): about 1.6 GB per worker with rhoHV
% cleaning and 2.6 GB with MistNet, so MistNet runs on fewer workers.
n_workers = 6;
if any(cellfun(@(x) isequal(string(x), "mistnet"), c.read_options))
    n_workers = 4;
end
load(fullfile(c.data_dir, "manifest.mat"), "raw_files")

%% Clean the volumes that are not cached yet --------------------------------
stage_tic = tic;
cache_dir = volume_cache_dir(c);
volume_files = strings(size(raw_files));
has_raw = raw_files ~= "";
raw_name = regexprep(raw_files(has_raw), ".*/", "");
volume_files(has_raw) = fullfile(cache_dir, raw_name + ".mat");
todo = find(has_raw & ~isfile(volume_files));

pool = gcp("nocreate");
if isempty(pool) || pool.NumWorkers ~= n_workers
    delete(pool);
    parpool("Processes", n_workers);
end
read_options = c.read_options;
vp = c.vp;
filter = c.filter;
status = strings(numel(todo), 1);
todo_raw = raw_files(todo);
todo_out = volume_files(todo);
parfor k = 1:numel(todo)
    status(k) = clean_volume(todo_raw(k), read_options, vp, filter, todo_out(k));
end
failed = status ~= "ok";
for k = find(failed)'
    fprintf("[prod/%s] failed %s: %s\n", c.name, todo_raw(k), status(k))
end
volume_files(todo(failed)) = "";
stage_log(c, "P2_clean_volumes", stage_tic, sprintf("%d cleaned, %d cached, %d failed, cache %s", ...
    nnz(~failed), nnz(has_raw) - numel(todo), nnz(failed), cache_dir))

save(fullfile(c.data_dir, "manifest.mat"), "volume_files", "-append")

% Idle workers keep their memory (reference maps, MistNet model); release it
% before the case build and inversion.
delete(gcp("nocreate"))
