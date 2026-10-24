function files = list_nexrad_day(station, day, bucket_url)
%LIST_NEXRAD_DAY List Level II volume files for one station and UTC day.
%   FILES = LIST_NEXRAD_DAY(STATION, DAY, BUCKET_URL) queries the public S3
%   bucket over HTTPS (no AWS CLI needed) and returns a table with the
%   object key, size in bytes, and scan start time (UTC, no time zone).
%   MDM metadata files are skipped.

prefix = sprintf("%s/%s/", string(day, "yyyy/MM/dd"), station);
keys = strings(0, 1);
sizes = zeros(0, 1);
continuation = "";
while true
    url = bucket_url + "/?list-type=2&prefix=" + prefix;
    if continuation ~= ""
        url = url + "&continuation-token=" + urlencode(continuation);
    end
    xml = string(webread(url, weboptions(ContentType="text", Timeout=60)));
    key_tokens = regexp(xml, "<Key>([^<]+)</Key>", "tokens");
    size_tokens = regexp(xml, "<Size>(\d+)</Size>", "tokens");
    keys = [keys; string(cellfun(@(c) c{1}, key_tokens, UniformOutput=false))']; %#ok<AGROW>
    sizes = [sizes; cellfun(@(c) str2double(c{1}), size_tokens)']; %#ok<AGROW>
    token = regexp(xml, "<NextContinuationToken>([^<]+)</NextContinuationToken>", "tokens", "once");
    if isempty(token)
        break
    end
    continuation = string(token{1});
end

% Parse scan start time from names such as KLOT20261022_020508_V06.
names = extractAfter(keys, prefix);
time_str = string(regexp(names, "\d{8}_\d{6}", "match", "once"));
keep = ~ismissing(time_str) & time_str ~= "" & ~contains(names, "_MDM");
keys = keys(keep);
sizes = sizes(keep);
scan_time = datetime(time_str(keep), InputFormat="yyyyMMdd_HHmmss");

files = table(keys(:), sizes(:), scan_time(:), VariableNames=["key" "bytes" "scan_time"]);
end
