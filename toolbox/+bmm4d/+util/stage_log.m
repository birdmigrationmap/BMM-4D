function stage_log(c, stage, stage_tic, note)
%STAGE_LOG Append one stage timing line to the case log and print it.
%   Writes time stamp, stage, duration and the MATLAB client resident memory
%   to <data_dir>/stage_log.csv. Worker memory is sampled separately (see
%   examples/README.md) because a client snapshot misses parfor peaks.

if nargin < 4
    note = "";
end
[~, rss_kb] = system(sprintf("ps -o rss= -p %d", feature("getpid")));
rss_gb = str2double(strtrim(rss_kb)) / 1024 ^ 2;
duration_s = toc(stage_tic);
log_file = fullfile(c.data_dir, "stage_log.csv");
if ~isfile(log_file)
    writelines("time,case,stage,seconds,client_rss_gb,note", log_file)
end
writelines(sprintf("%s,%s,%s,%.1f,%.2f,%s", string(datetime("now"), "yyyy-MM-dd HH:mm:ss"), ...
    c.name, stage, duration_s, rss_gb, note), log_file, WriteMode="append")
fprintf("[bmm4d/%s] %s | %.1f s | client %.2f GB | %s\n", c.name, stage, duration_s, rss_gb, note)
end
