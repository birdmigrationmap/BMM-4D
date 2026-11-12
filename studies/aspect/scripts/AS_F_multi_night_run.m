%% Aspect F: download, clean and bin the chosen nights in a background job (dev)
% Runs as_multi_night_job on the nights of AS_E_select_nights.m: per night,
% the us_east radars every hour for 6 h, cleaned as production but without
% the shape correction, then binned on the overlap cells (as_night_bins).
% A pilot night first (about 1.5 h), checked by hand, then the others
% (about 20 h for 15 nights). Both are batch jobs with 4 workers: MistNet
% workers hold 3.5-6 GB each on this 36 GB machine, so nothing else large
% (no other pool, no full inversion) should run beside them.
% Before running:
% - the user's OK for the download: the pilot about 5.7 GB, the main run
%   about 85 GB (AS_E prints the total; volumes already local are not
%   downloaded again);
% - the disk: the raw volumes stay in nexrad_raw_dir unless the user
%   agrees to delete_raw = true (only the files this job downloaded are
%   deleted, each night after its bins are saved); the cleaned volumes
%   take about 2 MB each (~21 GB for 16 nights) in the shared cache.
% The job stops before a night when less than 30 GB are free, and skips
% the nights already done (done.txt), so it can be stopped and restarted.
% Outputs: <data>/aspect/multi_night/<name>/frame_XX.mat, night.mat,
% run_log.csv.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"), ...
    fullfile(project_root, "studies", "artifacts", "functions"), fullfile(project_root, "examples"))
data_root = bmm4d.util.local_data_dir();
multi_dir = fullfile(data_root, "aspect", "multi_night");
nights = readtable(fullfile(multi_dir, "nights.csv"), TextType="string");
nights.first_frame = datetime(nights.first_frame);
delete_raw = false;     % true only with the user's OK

%% 1. Pilot: the strongest night
% Measures the download rate, the cleaning time per volume, the worker
% memory and the number of overlap pairs, before the long run.
pilot = nights(nights.median_refl == max(nights.median_refl), :);
fprintf("pilot %s: %d volumes, about %.1f GB to download\n", pilot.name, pilot.volumes, pilot.download_gb)
pilot_job = batch(@as_multi_night_job, 1, {pilot, project_root, multi_dir, delete_raw}, Pool=4, ...
    AutoAddClientPath=true, CurrentFolder=project_root);

%% 2. Check the pilot (once pilot_job.State is "finished")
% Times per stage and failures; then a first fit of the night alone,
% against the us_east night (b2_per_km 0.092 +- 0.005, b1 -0.21 +- 0.02).
disp(pilot_job.State)
run_log = readtable(fullfile(multi_dir, "run_log.csv"), TextType="string");
run_log.clean_s_per_volume = run_log.clean_s ./ max(run_log.n_cleaned, 1);
disp(run_log)
pilot_dir = fullfile(multi_dir, pilot.name);
frame_files = dir(fullfile(pilot_dir, "frame_*.mat"));
pair_parts = cell(numel(frame_files), 1);
volume_parts = cell(numel(frame_files), 1);
for i_frame = 1:numel(frame_files)
    loaded = load(fullfile(pilot_dir, frame_files(i_frame).name));
    pair_parts{i_frame} = as_overlap_pairs(loaded.bins);
    volume_parts{i_frame} = loaded.volumes;
end
pairs = vertcat(pair_parts{:});
volumes = [volume_parts{:}];
radar_height = nan(max([volumes.radar]), 1);
radar_height([volumes.radar]) = [volumes.height];
D = as_overlap_design(pairs, "raw");
pilot_fit = as_fit_aspect_model(D, true(size(D.y)), radar_height, "height");
fprintf("%s: %d cell pairs, %d radar pairs; b2_per_km %.3f (%.3f), b1 %.3f (%.3f)\n", pilot.name, ...
    numel(D.y), max(D.pair), pilot_fit.b2_per_km, pilot_fit.se(1), pilot_fit.b1, pilot_fit.se(2))

%% 3. Main run: the other nights (about 20 h)
% Restartable: nights with done.txt are skipped. Check progress with
% run_log.csv and main_job.State; the job's own log: main_job.Tasks(1).Diary
% (or diary(main_job)).
main_job = batch(@as_multi_night_job, 1, {nights, project_root, multi_dir, delete_raw}, Pool=4, ...
    AutoAddClientPath=true, CurrentFolder=project_root);

%% 4. Progress
disp(main_job.State)
run_log = readtable(fullfile(multi_dir, "run_log.csv"), TextType="string");
disp(run_log(:, ["name" "n_volumes" "n_downloaded" "download_s" "clean_s" "bin_s" "failures"]))
