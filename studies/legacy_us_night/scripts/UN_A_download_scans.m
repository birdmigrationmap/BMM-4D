%% US night reconstruction: download radar volumes
% Select the NEXRAD radars east of the Rocky Mountains, pick the volume
% closest to each requested time, and download them from the public
% Unidata S3 bucket over HTTPS.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "us_night", "functions"))
script_tic = tic;

%% Night controls ---------------------------------------------------------
% night_date is the local evening date. Times are UTC; 23:30-12:00 UTC
% covers dusk on the East Coast to dawn in the Great Plains.
night_date = datetime("2026-10-21");
time_step = minutes(30);
scan_time_requested = (night_date + hours(23.5)):time_step:(night_date + hours(36));
closest_accepted = minutes(10);

% Radars east of the Rockies. The western limit keeps the Front Range
% stations out to save memory; KFDX and KUDX are the westernmost kept.
radar_bounds = [-104 -66 24 50];
bucket_url = "https://unidata-nexrad-level2.s3.amazonaws.com";
n_parallel_download = 16;

night_tag = string(night_date, "yyyyMMdd");
night_dir = fullfile(local_data_dir(project_root), "us_night", "data", night_tag);
if ~exist(night_dir, "dir")
    mkdir(night_dir)
end
% Raw volumes go to the local cache outside OneDrive (local_settings.json).
raw_dir = nexrad_raw_dir();

%% Select radars ----------------------------------------------------------
radar_selected = select_radars_in_bbox(fullfile(project_root, "data", "reference", "radar.csv"), radar_bounds);
n_radar = height(radar_selected);
n_time = numel(scan_time_requested);
fprintf("[UN/A] night %s | %d radars, %d times (%s to %s UTC)\n", night_tag, n_radar, n_time, ...
    string(scan_time_requested(1), "dd-MMM HH:mm"), string(scan_time_requested(end), "dd-MMM HH:mm"))

%% List available files and pick the closest volume per time -------------
% One listing per radar and UTC day; the night spans two UTC days.
scan_days = unique(dateshift(scan_time_requested, "start", "day"));
scan_key = strings(n_radar, n_time);
scan_offset = nan(n_radar, n_time);
scan_bytes = zeros(n_radar, n_time);
for i_r = 1:n_radar
    station = radar_selected.name(i_r);
    files = cell(numel(scan_days), 1);
    for i_d = 1:numel(scan_days)
        try
            files{i_d} = list_nexrad_day(station, scan_days(i_d), bucket_url);
        catch err
            warning("UN_A:ListFailed", "Listing failed for %s on %s: %s", station, string(scan_days(i_d)), err.message)
        end
    end
    files = vertcat(files{:});
    if isempty(files) || height(files) == 0
        fprintf("[UN/A] %s | no files\n", station)
        continue
    end
    for i_t = 1:n_time
        [offset, id] = min(abs(files.scan_time - scan_time_requested(i_t)));
        if offset <= closest_accepted
            scan_key(i_r, i_t) = files.key(id);
            scan_offset(i_r, i_t) = minutes(offset);
            scan_bytes(i_r, i_t) = files.bytes(id);
        end
    end
    fprintf("[UN/A] %s | %d/%d times available\n", station, nnz(scan_key(i_r, :) ~= ""), n_time)
end

%% Download missing volumes in parallel -----------------------------------
% Files follow the cache layout <raw_dir>/<RADAR>_<yyyymmdd>/<file>, shared
% with load_rd. curl --parallel runs the downloads concurrently.
has_key = scan_key ~= "";
scan_name = regexprep(scan_key(has_key), ".*/", "");
scan_folder = extractBefore(scan_name, 5) + "_" + extractBetween(scan_name, 5, 12);
scan_file = strings(n_radar, n_time);
scan_file(has_key) = string(raw_dir) + filesep + scan_folder + filesep + scan_name;
to_download = has_key & ~isfile(scan_file);
fprintf("[UN/A] downloading %d files (%.1f GB), %d already present\n", ...
    nnz(to_download), sum(scan_bytes(to_download)) / 1e9, nnz(has_key & ~to_download))

if any(to_download(:))
    curl_config_file = fullfile(night_dir, "download_curl.txt");
    writelines(compose("url = ""%s/%s""\noutput = ""%s""", bucket_url, scan_key(to_download), scan_file(to_download)), ...
        curl_config_file)
    download_cmd = sprintf("curl --parallel --parallel-max %d --create-dirs -s -S -f --retry 3 -K '%s'", ...
        n_parallel_download, curl_config_file);
    [status, cmd_out] = system(download_cmd);
    if status ~= 0
        warning("UN_A:DownloadErrors", "Some downloads failed:\n%s", cmd_out)
    end
end

% Older archive files are gzipped; unpack them so read_rd can open them.
gz_id = find(endsWith(scan_file, ".gz") & isfile(scan_file));
for i = gz_id(:)'
    gunzip(scan_file(i), fileparts(scan_file(i)))
    delete(scan_file(i))
end
scan_file(endsWith(scan_file, ".gz")) = extractBefore(scan_file(endsWith(scan_file, ".gz")), ".gz");
scan_file(~isfile(scan_file)) = "";

%% Save manifest ----------------------------------------------------------
night_controls = struct( ...
    "night_date", night_date, ...
    "night_tag", night_tag, ...
    "scan_time_requested", scan_time_requested, ...
    "time_step", time_step, ...
    "closest_accepted", closest_accepted, ...
    "radar_bounds", radar_bounds, ...
    "bucket_url", bucket_url);
manifest_file = fullfile(night_dir, "manifest.mat");
save(manifest_file, "radar_selected", "scan_key", "scan_file", "scan_offset", "night_controls")

availability = mean(scan_file ~= "", 2);
fprintf("[UN/A] %d radars with all times, %d with none\n", nnz(availability == 1), nnz(availability == 0))
fprintf("[UN/A] saved %s in %.1f min\n", manifest_file, toc(script_tic) / 60)
