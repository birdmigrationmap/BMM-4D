%% Near radar A: scan strategies and cut types of the volumes on disk (dev)
% Which VCPs occur, and where the split cuts (a contiguous surveillance
% cut for reflectivity) end and the batch or Doppler cuts begin. The
% near-radar deficit follows that boundary (studies/artifacts, AR_J), so
% the cut type from the Level II header may be a better criterion than a
% fixed elevation. Reads only the metadata record of each raw file
% (nr_read_vcp, ~0.1 s per file).
% Sets: the us_east night (21-22 Oct 2026, one file per radar every 3 h),
% the spring night of 12-13 May 2021 (28 Midwest radars, hourly) and the
% 9 sandbox nights (one volume per radar).
% Output: figures/NR_A_vcp_volumes.csv (one row per volume),
% figures/NR_A_vcp_summary.csv (one row per VCP and boundary).

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "near_radar", "functions"))
figure_dir = fullfile(project_root, "studies", "near_radar", "figures");
raw_dir = bmm4d.radar.nexrad_raw_dir();
sandbox_raw_dir = fullfile(bmm4d.util.local_data_dir(), "sandbox_nights", "raw");

%% List the raw volumes
sandbox_files = string({dir(fullfile(sandbox_raw_dir, "K*_V06")).name})';
files = table(fullfile(sandbox_raw_dir, sandbox_files), repmat("sandbox_nights", numel(sandbox_files), 1), ...
    VariableNames=["file" "set"]);
night_sets = struct("set", {"us_east_20261021", "night_20210512"}, "pattern", {"K*_2026102*", "K*_20210513"}, ...
    "every_h", {3, 1});
for night = night_sets
    radar_days = dir(fullfile(raw_dir, night.pattern));
    for radar_day = radar_days(:)'
        volumes = dir(fullfile(radar_day.folder, radar_day.name, "K*V06"));
        if isempty(volumes), continue, end
        stamp = datetime(extractBetween(string({volumes.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
        % First volume of every every_h-hour block
        [~, first] = unique(floor((hour(stamp) + day(stamp) * 24) / night.every_h), "stable");
        files = [files; table(string(fullfile(radar_day.folder, radar_day.name, {volumes(first).name}))', ...
            repmat(string(night.set), numel(first), 1), VariableNames=["file" "set"])]; %#ok<AGROW>
    end
end
fprintf("%d raw volumes\n", height(files))

%% Read the scan strategy of each volume
n_file = height(files);
[vcp_number, n_split, top_split_deg, first_other_deg, other_types] = deal(nan(n_file, 1), nan(n_file, 1), ...
    nan(n_file, 1), nan(n_file, 1), strings(n_file, 1));
read_tic = tic;
for i_f = 1:n_file
    try
        vcp = nr_read_vcp(files.file(i_f));
    catch err
        warning("%s: %s", files.file(i_f), err.message)
        continue
    end
    if isempty(vcp.cuts), continue, end
    cuts = vcp.cuts(vcp.cuts.elev < 7, :);
    vcp_number(i_f) = vcp.number;
    split_elev = unique(cuts.elev(cuts.split));
    n_split(i_f) = numel(split_elev);
    top_split_deg(i_f) = max(split_elev);
    other = cuts(~cuts.split, :);
    if ~isempty(other)
        first_other_deg(i_f) = min(other.elev);
        other_types(i_f) = strjoin(unique(other.waveform_name)', "+");
    end
end
fprintf("read in %.0f s\n", toc(read_tic))
[~, name] = fileparts(files.file);
files.radar = extractBefore(name, 5);
files.time = datetime(extractBetween(name, 5, 19), InputFormat="yyyyMMdd_HHmmss");
files.vcp = vcp_number;
files.n_split = n_split;
files.top_split_deg = top_split_deg;
files.first_other_deg = first_other_deg;
files.other_types = other_types;
writetable(files(:, 2:end), fullfile(figure_dir, "NR_A_vcp_volumes.csv"))

%% Summary: VCPs, boundary and cut types of the non-split cuts below 7 deg
summary = groupsummary(files(~isnan(files.vcp), :), ["vcp" "top_split_deg" "first_other_deg" "other_types"]);
summary.radars = arrayfun(@(i) numel(unique(files.radar(files.vcp == summary.vcp(i) & ...
    files.top_split_deg == summary.top_split_deg(i)))), (1:height(summary))');
summary.sets = arrayfun(@(i) strjoin(unique(files.set(files.vcp == summary.vcp(i) & ...
    files.top_split_deg == summary.top_split_deg(i)))', " "), (1:height(summary))');
summary = sortrows(summary, "GroupCount", "descend");
disp(summary)
writetable(summary, fullfile(figure_dir, "NR_A_vcp_summary.csv"))

% VCP changes during a night (a radar switching scan strategy)
per_radar_night = groupsummary(files(~isnan(files.vcp) & files.set ~= "sandbox_nights", :), ["set" "radar"], ...
    @(v) numel(unique(v)), "vcp");
fprintf("radar-nights with more than one VCP: %d of %d\n", nnz(per_radar_night{:, end} > 1), height(per_radar_night))
