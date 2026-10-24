function volume_files = clean_volumes(c)
%CLEAN_VOLUMES Clean every raw volume of a case into the shared volume cache.
%   VOLUME_FILES = BMM4D.CLEAN_VOLUMES(C) cleans each volume of the manifest
%   (get_volumes) with bmm4d.radar.clean_volume and the cleaning settings
%   of C (read_options, vp, filter), in parallel. Cleaned volumes are kept
%   in a cache folder named after those settings (volume_cache_dir), so
%   every case with the same cleaning reuses them. VOLUME_FILES
%   [radar time] is "" where missing or failed; it is appended to the
%   manifest.

load(fullfile(c.data_dir, "manifest.mat"), "raw_files")
% Measured worker memory (Key West): about 1.6 GB per worker with rhoHV
% cleaning and 2.6 GB with MistNet, so MistNet runs on fewer workers.
n_workers = 6;
if any(cellfun(@(x) isequal(string(x), "mistnet"), c.read_options))
    n_workers = 4;
end

%% Clean the volumes that are not cached yet
stage_tic = tic;
cache_dir = bmm4d.radar.volume_cache_dir(c);
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
    status(k) = bmm4d.radar.clean_volume(todo_raw(k), read_options, vp, filter, todo_out(k));
end
failed = status ~= "ok";
for k = find(failed)'
    fprintf("[bmm4d/%s] failed %s: %s\n", c.name, todo_raw(k), status(k))
end
volume_files(todo(failed)) = "";
bmm4d.util.stage_log(c, "clean_volumes", stage_tic, sprintf("%d cleaned, %d cached, %d failed, cache %s", ...
    nnz(~failed), nnz(has_raw) - numel(todo), nnz(failed), cache_dir))

save(fullfile(c.data_dir, "manifest.mat"), "volume_files", "-append")

% Idle workers keep their memory (reference maps, MistNet model); release it
% before the case build and inversion.
delete(gcp("nocreate"))
end
