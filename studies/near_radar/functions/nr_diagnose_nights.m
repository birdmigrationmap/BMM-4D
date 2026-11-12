function out_files = nr_diagnose_nights(night_names, n_workers, max_frames)
%NR_DIAGNOSE_NIGHTS Near-radar diagnostic of every volume of shared nights (dev).
%   OUT_FILES = NR_DIAGNOSE_NIGHTS(NIGHT_NAMES, N_WORKERS) runs, for each
%   night of the shared night set (studies/nights: products/prod/<name>,
%   cleaned volumes in cache 493fd773, raw files kept), nr_read_vcp and
%   nr_volume_ratio on every cleaned volume of its manifest, on a pool of
%   N_WORKERS (default 6; 0 runs serially), and saves the rows to
%   products/near_radar/NR_E/rows_<name>.mat. A night whose file exists is
%   skipped, so a cut-off job is finished by submitting it again.
%   MAX_FRAMES (default Inf) keeps the first frames of each night only (a
%   pilot); its rows go to rows_<name>_pilot.mat, so a full run never skips
%   a night because of a pilot.
%   Meant for a cluster job (hpc/submit.sh); also runs on the laptop with
%   the share mounted (nr_shared_path). ~0.3-0.5 s per volume per worker.

if nargin < 2
    n_workers = 6;
end
if nargin < 3
    max_frames = inf;
end
suffix = "";
if isfinite(max_frames)
    suffix = "_pilot";
end
products_dir = nr_shared_path("/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/BMM-4D/products");
out_dir = fullfile(products_dir, "near_radar", "NR_E");
if ~isfolder(out_dir)
    mkdir(out_dir)
end
% The headers are read through python3: fail early if it cannot decompress
night_names = string(night_names(:))';
first = load(fullfile(products_dir, "prod", night_names(1), "manifest.mat"), "raw_files");
probe = nr_shared_path(first.raw_files(find(first.raw_files ~= "", 1)));
assert(~isempty(nr_read_vcp(probe).cuts), "nr_diagnose_nights:Header", "No scan strategy read from %s.", probe)
if n_workers > 0
    % Workers need wsrlib (slant2ground) and these functions on their path
    % (parfevalOnAll, not spmd: spmd needs MPI networking, which fails off-site)
    if isempty(gcp("nocreate"))
        project_root = setup_paths();
        parpool("Processes", n_workers);
        wait(parfevalOnAll(@setup_paths, 0, project_root));
    end
    wait(parfevalOnAll(@addpath, 0, fileparts(mfilename("fullpath"))));
end
out_files = strings(size(night_names));
for i_n = 1:numel(night_names)
    name = night_names(i_n);
    out_files(i_n) = fullfile(out_dir, "rows_" + name + suffix + ".mat");
    if isfile(out_files(i_n))
        continue
    end
    night_tic = tic;
    manifest = load(fullfile(products_dir, "prod", name, "manifest.mat"), "raw_files", "volume_files");
    frames = 1:min(size(manifest.volume_files, 2), max_frames);
    volume_files = nr_shared_path(reshape(manifest.volume_files(:, frames), [], 1));
    raw_files = nr_shared_path(reshape(manifest.raw_files(:, frames), [], 1));
    todo = volume_files ~= "" & raw_files ~= "";
    volume_files = volume_files(todo);
    raw_files = raw_files(todo);
    parts = cell(numel(volume_files), 1);
    parfor (k = 1:numel(volume_files), n_workers)
        try
            vcp = nr_read_vcp(raw_files(k));
            if isempty(vcp.cuts), continue, end
            rows = nr_volume_ratio(load(volume_files(k), "volume").volume, vcp.cuts);
            if isempty(rows), continue, end
            rows.set = repmat(name, height(rows), 1);
            rows.vcp = repmat(vcp.number, height(rows), 1);
            parts{k} = rows;
        catch err
            warning("nr_diagnose_nights:Volume", "%s: %s", volume_files(k), err.message)
        end
    end
    rows = vertcat(parts{:});
    n_volume = nnz(~cellfun(@isempty, parts));
    save(out_files(i_n), "rows", "n_volume")
    fprintf("[near_radar] %s: %d of %d volumes, %d rows, %.0f s\n", name, n_volume, numel(volume_files), height(rows), toc(night_tic))
end
end
