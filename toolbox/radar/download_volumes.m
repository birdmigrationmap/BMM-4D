function n_downloaded = download_volumes(radars, times, need, closest_accepted_min, n_parallel)
%DOWNLOAD_VOLUMES Download the volumes closest to the requested times.
%   NEED(i_r, i_t) marks radar-time pairs without a local file. Files are
%   listed over HTTPS from the public Unidata bucket (one listing per radar
%   and UTC day) and fetched concurrently with curl into the raw cache
%   (nexrad_raw_dir, layout <RADAR>_<yyyymmdd>/<file>), where load_rd finds
%   them.

bucket_url = "https://unidata-nexrad-level2.s3.amazonaws.com";
raw_dir = nexrad_raw_dir();
urls = strings(0, 1);
outputs = strings(0, 1);
for i_r = find(any(need, 2))'
    % UTC days to list, including the neighbouring day for times near midnight.
    window = [times(need(i_r, :)) - minutes(closest_accepted_min), times(need(i_r, :)) + minutes(closest_accepted_min)];
    days = unique(dateshift(window, "start", "day"));
    listing = cell(numel(days), 1);
    for i_d = 1:numel(days)
        try
            listing{i_d} = list_nexrad_day(radars(i_r), days(i_d), bucket_url);
        catch
            listing{i_d} = [];
        end
    end
    files = vertcat(listing{:});
    if isempty(files) || height(files) == 0
        continue
    end
    for i_t = find(need(i_r, :))
        [offset, id] = min(abs(files.scan_time - times(i_t)));
        if minutes(offset) < closest_accepted_min
            name = regexprep(files.key(id), ".*/", "");
            urls(end + 1, 1) = bucket_url + "/" + files.key(id); %#ok<AGROW>
            outputs(end + 1, 1) = fullfile(raw_dir, radars(i_r) + "_" + string(files.scan_time(id), "yyyyMMdd"), name); %#ok<AGROW>
        end
    end
end

[outputs, keep] = unique(outputs);
urls = urls(keep);
n_downloaded = numel(urls);
if n_downloaded == 0
    return
end
config_file = string(tempname) + ".txt";
writelines(compose("url = ""%s""\noutput = ""%s""", urls, outputs), config_file)
[status, out] = system(sprintf("curl --parallel --parallel-max %d --create-dirs -s -S -f --retry 3 -K '%s'", ...
    n_parallel, config_file));
if status ~= 0
    warning("download_volumes:CurlErrors", "Some downloads failed:\n%s", out)
end
delete(config_file)
end
