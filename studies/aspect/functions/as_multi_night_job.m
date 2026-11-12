function run_log = as_multi_night_job(nights, project_root, multi_dir, delete_raw)
%AS_MULTI_NIGHT_JOB Download, clean and bin several nights for the aspect fit (dev).
%   RUN_LOG = AS_MULTI_NIGHT_JOB(NIGHTS, PROJECT_ROOT, MULTI_DIR, DELETE_RAW)
%   is meant for a batch job with a pool:
%       batch(@as_multi_night_job, 1, {nights, project_root, multi_dir, false}, ...
%           Pool=4, AutoAddClientPath=true, CurrentFolder=project_root)
%   For each row of NIGHTS (AS_E_select_nights.m: name, first_frame,
%   n_frames), on the 107 radars of the us_east case:
%   1. bmm4d.config of the night (hourly frames from first_frame, cleaning
%      as us_east but without the shape correction: read option "shape"
%      false, so the volumes carry no aspect correction);
%   2. bmm4d.get_volumes: local raw files, the missing ones downloaded;
%   3. ar_clean_volumes (studies/artifacts): bmm4d.radar.clean_volume on
%      the job's pool, into the shared volume cache of these settings;
%   4. as_night_bins of every frame: MULTI_DIR/<name>/frame_XX.mat;
%   5. DELETE_RAW true: deletes the raw files this job downloaded for the
%      night (never files that were there before). Only with the user's OK.
%   A night is skipped when MULTI_DIR/<name>/done.txt exists, so the job
%   can be stopped and restarted. The job stops before a night when less
%   than 30 GB are free on the data disk. One row per night is appended to
%   MULTI_DIR/run_log.csv (times in s, volumes, downloads, failures).

setup_paths(project_root);
aspect_functions = fullfile(project_root, "studies", "aspect", "functions");
addpath(aspect_functions, fullfile(project_root, "studies", "artifacts", "functions"), ...
    fullfile(project_root, "examples"))
spmd
    setup_paths(project_root);
    addpath(aspect_functions)
end
us = case_settings("us_east");
log_file = fullfile(multi_dir, "run_log.csv");
run_log = table();
for k = 1:height(nights)
    name = string(nights.name(k));
    night_dir = fullfile(multi_dir, name);
    if isfile(fullfile(night_dir, "done.txt"))
        continue
    end
    if ~isfolder(night_dir)
        mkdir(night_dir)
    end
    free_gb = java.io.File(char(multi_dir)).getUsableSpace() / 1e9;
    if free_gb < 30
        warning("as_multi_night_job:Disk", "%.0f GB free: stopped before %s.", free_gb, name)
        break
    end

    % 1. Settings of the night: us_east radars and cleaning, no shape correction
    times = nights.first_frame(k) + hours(0:nights.n_frames(k) - 1);
    c = bmm4d.config(name, "radars", us.radars, "times", times, "grid_bounds", us.grid_bounds, ...
        "filter", us.filter, "solve_aggregate", us.solve_aggregate, "domain", us.domain, ...
        "read_options", {"shape", false});

    % 2. Raw volumes: what is local, then the download
    before = bmm4d.radar.load_rd(c.radars, c.times, closest_accepted=c.closest_accepted_min);
    stage_tic = tic;
    raw_files = bmm4d.get_volumes(c);
    download_s = toc(stage_tic);
    downloaded = raw_files ~= "" & before == "";

    % 3. Cleaning on the job's pool
    stage_tic = tic;
    [n_cleaned, failures] = ar_clean_volumes(c, project_root);
    clean_s = toc(stage_tic);
    load(fullfile(c.data_dir, "manifest.mat"), "volume_files")

    % 4. Overlap bins of every frame
    stage_tic = tic;
    frame_times = c.times;
    n_cells = zeros(numel(frame_times), 1);
    parfor i_frame = 1:numel(frame_times)
        bins = as_night_bins(volume_files(:, i_frame), frame_times(i_frame), i_frame, ...
            fullfile(night_dir, sprintf("frame_%02d.mat", i_frame)));
        n_cells(i_frame) = height(bins);
    end
    bin_s = toc(stage_tic);
    save(fullfile(night_dir, "night.mat"), "frame_times", "volume_files")

    % 5. Raw files downloaded by this job (only with the user's OK)
    if delete_raw
        delete(raw_files(downloaded))
    end

    row = table(name, numel(frame_times), nnz(raw_files ~= ""), nnz(downloaded), n_cleaned, ...
        nnz(volume_files ~= ""), sum(n_cells), download_s, clean_s, bin_s, "", ...
        VariableNames=["name" "n_frames" "n_raw" "n_downloaded" "n_cleaned" "n_volumes" "n_bin_rows" ...
        "download_s" "clean_s" "bin_s" "failures"]);
    row.failures = strjoin(string(failures(:))', " | ");
    run_log = [run_log; row]; %#ok<AGROW>
    writetable(row, log_file, WriteMode="append")
    writelines(string(datetime("now")), fullfile(night_dir, "done.txt"))
end
end
