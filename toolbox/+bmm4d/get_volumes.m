function raw_files = get_volumes(c)
%GET_VOLUMES Locate or download the raw radar volumes of a case.
%   RAW_FILES = BMM4D.GET_VOLUMES(C) finds the volume closest to each
%   radar-time pair of C (bmm4d.config) in the local raw caches and
%   downloads the missing ones in parallel. Only night volumes
%   (C.radar_night) are kept; RAW_FILES is [radar time], "" where missing.
%   Saves radars, times and raw_files to <C.data_dir>/manifest.mat.

if ~isfolder(c.data_dir)
    mkdir(c.data_dir)
end

%% Find local files and download missing ones
stage_tic = tic;
raw_files = bmm4d.radar.load_rd(c.radars, c.times, closest_accepted=c.closest_accepted_min);
raw_files(~c.radar_night) = "";
need = raw_files == "" & c.radar_night;
n_downloaded = 0;
if any(need(:))
    n_downloaded = bmm4d.radar.download_volumes(c.radars, c.times, need, c.closest_accepted_min, 16);
    raw_files = bmm4d.radar.load_rd(c.radars, c.times, closest_accepted=c.closest_accepted_min);
    raw_files(~c.radar_night) = "";
end
bmm4d.util.stage_log(c, "get_volumes", stage_tic, sprintf("%d/%d night volumes, %d frames, %d downloaded", ...
    nnz(raw_files ~= ""), nnz(c.radar_night), numel(c.times), n_downloaded))

%% Save the manifest
radars = c.radars;
times = c.times;
save(fullfile(c.data_dir, "manifest.mat"), "radars", "times", "raw_files")
missing_radar = c.radars(all(raw_files == "", 2));
if ~isempty(missing_radar)
    fprintf("[bmm4d/%s] no data: %s\n", c.name, strjoin(missing_radar, ", "))
end
end
