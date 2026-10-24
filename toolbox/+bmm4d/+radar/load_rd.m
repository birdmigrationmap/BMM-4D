function scans_filepath = load_rd(scan_radar, scan_time, varargin)
%LOAD_RD Locate or download raw radar files for requested radar-time pairs.
%   SCANS_FILEPATH = LOAD_RD(SCAN_RADAR, SCAN_TIME, ...) returns a string
%   array with one filepath per requested radar-time combination, or ""
%   when no file lies within the accepted offset.
%
%   Files are searched in the local cache (nexrad_raw_dir) and in
%   raw_scans/nexrad of the data folder (data_path) used by earlier
%   workflows. Downloads always go to the local cache.
%
%   Options:
%   - 'closest_accepted': maximum time offset in minutes.
%   - 'download': download missing files from the Unidata S3 bucket when true.
%   - 'search_dirs': folders to search instead of the two above.

p = inputParser;
addParameter(p,'closest_accepted',10);
addParameter(p,'download',false);
addParameter(p,'search_dirs',strings(0));
parse(p,varargin{:});

% A char radar name such as 'KMPX' would otherwise loop over its letters.
scan_radar = string(scan_radar);
bucket_url = "https://unidata-nexrad-level2.s3.amazonaws.com";
cache_dir = bmm4d.radar.nexrad_raw_dir();
search_dirs = unique([string(cache_dir) string(bmm4d.util.data_path("raw_scans", "nexrad"))], "stable");
if ~isempty(p.Results.search_dirs)
    search_dirs = string(p.Results.search_dirs);
end

% List of unique date
[scan_day,~,scan_day_id] = unique(dateshift(scan_time,"start","day"));

scans_filepath=strings(numel(scan_radar),numel(scan_time));

% Also search the neighbouring days so requests near midnight can match a
% file stored in the previous or next day's folder.
search_days = unique([scan_day(:) - 1; scan_day(:); scan_day(:) + 1]);

for i_s=1:numel(scan_radar)
    % Local files from every search folder, with time parsed from the name
    local_files = strings(0, 1);
    for i_day = 1:numel(search_days)
        folder_name = scan_radar(i_s) + "_" + string(search_days(i_day), "yyyyMMdd");
        for i_dir = 1:numel(search_dirs)
            fl = dir(fullfile(search_dirs(i_dir), folder_name, scan_radar(i_s) + "*"));
            fl = fl(~contains({fl.name}, "_MDM") & ~endsWith({fl.name}, ".gz"));
            local_files = [local_files; string(fullfile({fl.folder}, {fl.name}))']; %#ok<AGROW>
        end
    end
    local_time = parse_scan_time(local_files);

    for i_t = 1:numel(scan_time)
        [offset, id] = min(abs(local_time - scan_time(i_t)));
        if ~isempty(id) && minutes(offset) < p.Results.closest_accepted
            scans_filepath(i_s, i_t) = local_files(id);
        end
    end

    for i_d=1:numel(scan_day)
        folder_name = scan_radar(i_s) + "_" + string(scan_day(i_d), "yyyyMMdd");
        id_t = find(scan_day_id==i_d);

        missing_t = id_t(scans_filepath(i_s, id_t) == "");
        if ~p.Results.download || isempty(missing_t)
            continue
        end

        try
            remote = bmm4d.radar.list_nexrad_day(scan_radar(i_s), scan_day(i_d), bucket_url);
        catch err
            warning("load_rd:ListFailed", "Could not list files for %s on %s: %s", ...
                scan_radar(i_s), string(scan_day(i_d), "yyyy-MM-dd"), err.message)
            continue
        end
        if height(remote) == 0
            warning("load_rd:NoData", "No data for %s on %s.", scan_radar(i_s), string(scan_day(i_d), "yyyy-MM-dd"))
            continue
        end

        target_dir = fullfile(cache_dir, folder_name);
        if ~isfolder(target_dir)
            mkdir(target_dir)
        end
        for i_t = missing_t(:)'
            [offset, id] = min(abs(remote.scan_time - scan_time(i_t)));
            if minutes(offset) >= p.Results.closest_accepted
                warning("load_rd:NoCloseFile", "No file within %g minutes for %s at %s.", ...
                    p.Results.closest_accepted, scan_radar(i_s), string(scan_time(i_t)))
                continue
            end
            file = fullfile(target_dir, regexprep(remote.key(id), ".*/", ""));
            try
                websave(file, bucket_url + "/" + remote.key(id));
            catch err
                % Do not leave a partial file that later looks like a valid scan.
                if isfile(file), delete(file), end
                warning("load_rd:DownloadFailed", "Could not download %s: %s", remote.key(id), err.message)
                continue
            end
            if endsWith(file, ".gz")
                gunzip(file, target_dir)
                delete(file)
                file = extractBefore(file, ".gz");
            end
            scans_filepath(i_s, i_t) = file;
        end
    end
end
end

function scan_time = parse_scan_time(files)
% Parse scan start time from names such as KLOT20261022_020508_V06.
time_str = string(regexp(files, "\d{8}_\d{6}", "match", "once"));
scan_time = NaT(numel(files), 1);
valid = ~ismissing(time_str) & time_str ~= "";
scan_time(valid) = datetime(time_str(valid), InputFormat="yyyyMMdd_HHmmss");
end
