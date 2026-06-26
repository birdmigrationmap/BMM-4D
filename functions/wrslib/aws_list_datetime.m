function fileinfo = aws_list_datetime(station, time)
%AWS_LIST_DATETIME List archive files for one station-day or station-hour.

y = year(time);
m = month(time);
d = day(time);
h = hour(time);

s3path = sprintf('%04d/%02d/%02d/%04s', y, m, d, station);
if h > 0
    s3path = sprintf('%s/%04s%04d%02d%02d_%02d', ...
        s3path, station, y, m, d, h);
end

aws_bin = "/usr/local/bin/aws";
cmd = sprintf([ ...
    'AWS_PAGER="" %s s3api list-objects ' ...
    '--bucket unidata-nexrad-level2 --prefix %s --max-items 10000 ' ...
    '--query ''Contents[].{Key: Key, Size: Size}'' --output json --no-sign-request'], ...
    aws_bin, s3path);

[status, result] = system(cmd);
if status ~= 0
    error("aws_list_datetime:AwsCliFailed", ...
        "AWS list failed for %s on %s. Command: %s. Output: %s", ...
        station, string(datetime(y, m, d)), cmd, strtrim(result))
end

result = strtrim(result);
if result == "" || result == "null"
    fileinfo = [];
    return
end

filedata = jsondecode(result);
if isempty(filedata)
    fileinfo = [];
    return
end
if ~isstruct(filedata)
    error("aws_list_datetime:UnexpectedResponse", ...
        "Unexpected AWS response type for %s on %s.", ...
        station, string(datetime(y, m, d)))
end

fileinfo = arrayfun(@(entry) aws_parse(entry.Key), filedata, UniformOutput=true);
end
