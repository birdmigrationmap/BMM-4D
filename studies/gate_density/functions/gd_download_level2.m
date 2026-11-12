function raw_file = gd_download_level2(station, target_time, out_dir, varargin)
%GD_DOWNLOAD_LEVEL2 Download the NEXRAD Level II volume closest to a time.
%   RAW_FILE = GD_DOWNLOAD_LEVEL2(STATION, TARGET_TIME, OUT_DIR) lists the
%   public archive (unidata-nexrad-level2, over https) for that station and
%   UTC day, picks the volume closest to TARGET_TIME (within MAX_OFFSET_MIN),
%   downloads it to OUT_DIR (gunzipped if needed) and returns its path, or
%   "" when there is none. Network calls are retried N_TRY times.
%   Options: max_offset_min (10), n_try (3).

p = inputParser;
addParameter(p, "max_offset_min", 10);
addParameter(p, "n_try", 3);
parse(p, varargin{:});
opt = p.Results;

bucket = "https://unidata-nexrad-level2.s3.amazonaws.com/";
target_time.TimeZone = "";
prefix = sprintf("%s/%s/", string(target_time, "yyyy/MM/dd"), station);
listing = retry(@() webread(bucket + "?list-type=2&prefix=" + prefix, weboptions(Timeout=60)), opt.n_try);
keys = string(regexp(char(listing), "<Key>([^<]+)</Key>", "tokens"));
keys = keys(~endsWith(keys, "_MDM") & contains(keys, "_V0"));
raw_file = "";
if isempty(keys)
    return
end
% Volume start time from the file name, e.g. KDVN20230510_030436_V06
names = extractAfter(keys, prefix);
stamp = extractBetween(names, strlength(station) + 1, strlength(station) + 15);
key_time = datetime(stamp, InputFormat="yyyyMMdd_HHmmss");
[offset, i_k] = min(abs(key_time - target_time));
if minutes(offset) > opt.max_offset_min
    return
end
[~, name, ext] = fileparts(keys(i_k));
raw_file = fullfile(out_dir, name + ext);
if ~isfile(raw_file)
    retry(@() websave(raw_file, bucket + keys(i_k), weboptions(Timeout=120)), opt.n_try);
end
if endsWith(raw_file, ".gz")
    gunzip(raw_file, out_dir);
    delete(raw_file)
    raw_file = extractBefore(raw_file, strlength(raw_file) - 2);
end
end

function out = retry(fn, n_try)
% Call FN up to N_TRY times, waiting a little longer after each failure
for k = 1:n_try
    try
        out = fn();
        return
    catch err
        if k == n_try
            rethrow(err)
        end
        pause(10 * k)
    end
end
end
