function [n_cleaned, failures] = ar_clean_volumes(c, project_root, varargin)
%AR_CLEAN_VOLUMES bmm4d.clean_volumes for background jobs and dev filters (dev).
%   [N_CLEANED, FAILURES] = AR_CLEAN_VOLUMES(C, PROJECT_ROOT, ...) cleans
%   the volumes of the case that are not in the cache yet, as
%   bmm4d.clean_volumes, with the pool that is already open (the pool of a
%   batch job: batch(@ar_clean_volumes, 2, {c, project_root, ...}, Pool=4)),
%   and writes volume_files to the manifest. The job and its workers first
%   run setup_paths (toolbox and wsrlib). FAILURES lists the distinct error
%   messages. Options:
%   - "dev_folder": a folder put first on the path of the job and its
%     workers, whose +bmm4d functions replace the toolbox ones (e.g.
%     studies/artifacts/dev: the beam-blockage azimuth fix);
%   - "cache_dir": the cache folder (default bmm4d.radar.volume_cache_dir,
%     which only knows the read options: give a dev filter its own folder);
%   - "aspect": true cleans with ar_clean_volume (symmetric aspect
%     correction from the track; the read options must have shape false)
%     instead of bmm4d.radar.clean_volume.

p = inputParser;
addParameter(p, "dev_folder", "");
addParameter(p, "cache_dir", "");
addParameter(p, "aspect", false);
parse(p, varargin{:});
aspect = p.Results.aspect;
dev_folder = string(p.Results.dev_folder);
cache_dir = string(p.Results.cache_dir);
if cache_dir == ""
    cache_dir = bmm4d.radar.volume_cache_dir(c);
end
if ~isfolder(cache_dir)
    mkdir(cache_dir)
end

setup_paths(project_root);
if dev_folder ~= ""
    addpath(dev_folder, "-begin");
end
spmd
    setup_paths(project_root);
    if dev_folder ~= ""
        addpath(dev_folder, "-begin");
    end
end

load(fullfile(c.data_dir, "manifest.mat"), "raw_files")
volume_files = strings(size(raw_files));
has_raw = raw_files ~= "";
raw_name = regexprep(raw_files(has_raw), ".*/", "");
volume_files(has_raw) = fullfile(cache_dir, raw_name + ".mat");
todo = find(has_raw & ~isfile(volume_files));
read_options = c.read_options;
vp = c.vp;
filter = c.filter;
todo_raw = raw_files(todo);
todo_out = volume_files(todo);
status = strings(numel(todo), 1);
stage_tic = tic;
parfor k = 1:numel(todo)
    if aspect
        status(k) = ar_clean_volume(todo_raw(k), read_options, vp, filter, todo_out(k));
    else
        status(k) = bmm4d.radar.clean_volume(todo_raw(k), read_options, vp, filter, todo_out(k));
    end
end
failed = status ~= "ok";
volume_files(todo(failed)) = "";
n_cleaned = nnz(~failed);
failures = unique(status(failed));
bmm4d.util.stage_log(c, "clean_volumes", stage_tic, sprintf("%d cleaned, %d cached, %d failed, cache %s", ...
    n_cleaned, nnz(has_raw) - numel(todo), nnz(failed), cache_dir))
save(fullfile(c.data_dir, "manifest.mat"), "volume_files", "-append")
end
