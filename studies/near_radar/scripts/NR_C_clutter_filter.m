%% Near radar C: the clutter filter at the gates that read low (dev)
% Is the near-radar deficit carried by the gates the clutter filter took
% power from? If so, the filter's own per-gate flag (Level II CFP) is a
% rule without a range threshold. Per gate of the cleaned volume:
% density / predicted (nr_volume_ratio: the split-cut profile, beam
% weighted) and the CFP code of the raw cut the gate came from
% (bmm4d.radar.read_level2_moments; codes: 0 below threshold, 1 range
% folded, 2 filter not applied, 3 point clutter, >= 8: removed dB = code - 8).
% Classes: not applied (2), applied with 0 dB removed (8), 1-3 dB (9-11),
% >= 4 dB (>= 12). The split cuts are the control.
% Caveat: where the filter acts is set by the bypass map and the clutter
% decision; on cuts where it depends on the echo (the split cuts' CMD),
% gates with removed power are a selection, not a random subset.
% Volumes: one per radar near 05 UTC of the us_east night (cache
% 9c2b4ef2), of the 13 May 2021 night and of two sandbox nights (cache
% fad5e1e1); ~4 s per volume to read REF and CFP (~12 min).
% Outputs: figures/NR_C_cfp.csv (waveform, ring, class: gates, ratio),
% figures/NR_C_cfp_share.csv (share of gates by class per ring).

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "near_radar", "functions"))
figure_dir = fullfile(project_root, "studies", "near_radar", "figures");
volume_root = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes");
sandbox_raw_dir = fullfile(bmm4d.util.local_data_dir(), "sandbox_nights", "raw");
class_edges = [2 3 8 9 12 256];   % 2 | 3-7 | 8 | 9-11 | >= 12
class_names = ["not applied" "point/other" "0 dB" "1-3 dB" ">= 4 dB"];

%% Volumes: one per radar near 05 UTC
patterns = ["9c2b4ef2/K*20261022_05*_V06.mat" "fad5e1e1/K*20210513_05*_V06.mat" ...
    "fad5e1e1/K*20230517_*_V06.mat" "fad5e1e1/K*20211022_*_V06.mat"];
files = strings(0, 1);
for pattern = patterns
    listed = dir(fullfile(volume_root, pattern));
    names = string({listed.name})';
    % First volume of each radar
    [~, first] = unique(extractBefore(names, 5), "first");
    files = [files; fullfile(string({listed(first).folder})', names(first))]; %#ok<AGROW>
end
fprintf("%d volumes\n", numel(files))

%% Gates with their CFP code
parts = cell(numel(files), 1);
run_tic = tic;
for i_f = 1:numel(files)
    volume = load(files(i_f), "volume").volume;
    raw_file = volume.source_file;
    if ~isfile(raw_file)
        [~, raw_name] = fileparts(files(i_f));
        raw_file = fullfile(sandbox_raw_dir, raw_name);
    end
    if ~isfile(raw_file), continue, end
    vcp = nr_read_vcp(raw_file);
    [~, ~, gates] = nr_volume_ratio(volume, vcp.cuts);
    if isempty(gates), continue, end
    cuts = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], 7);
    % Cut type of each sweep, as in nr_volume_ratio
    sweep_waveform = strings(numel(volume.dens.sweeps), 1);
    for k = 1:numel(volume.dens.sweeps)
        [gap, i_cut] = min(abs(vcp.cuts.elev - volume.dens.sweeps(k).elev));
        if gap > 0.1, continue, end
        sweep_waveform(k) = vcp.cuts.waveform_name(i_cut);
        if vcp.cuts.split(i_cut), sweep_waveform(k) = "split"; end
    end
    gates.code = nan(height(gates), 1);
    for k = unique(gates.sweep)'
        sweep = volume.dens.sweeps(k);
        % Raw cut of this sweep: same radials, closest elevation, with CFP (the surveillance cut of a split pair)
        candidates = find(arrayfun(@(cut) numel(cut.az) == numel(sweep.az) && ~isempty(cut.CFP_code), cuts));
        if isempty(candidates), continue, end
        [gap, best] = min(abs([cuts(candidates).elev] - sweep.elev));
        if gap > 0.2, continue, end
        cut = cuts(candidates(best));
        mine = find(gates.sweep == k);
        [~, i_az] = min(abs(mod(gates.az_deg(mine)' - cut.az(:) + 180, 360) - 180), [], 1);
        i_rg = interp1(double(cut.CFP_range), 1:numel(cut.CFP_range), gates.range_m(mine), "nearest", NaN);
        ok = ~isnan(i_rg);
        gates.code(mine(ok)) = double(cut.CFP_code(sub2ind(size(cut.CFP_code), i_rg(ok), i_az(ok)')));
    end
    gates.waveform = sweep_waveform(gates.sweep);
    gates.station = repmat(string(volume.station), height(gates), 1);
    parts{i_f} = gates;
    fprintf("%d / %d %s, %.0f s\n", i_f, numel(files), volume.station, toc(run_tic))
end
gates = vertcat(parts{:});

%% Ratio and share of gates by cut type, ring and filter class
ring_edges = [3 5 7.5 10 15 20 30 40 60 100];
gates.ring = discretize(gates.range_m / 1000, ring_edges);
gates.class = discretize(gates.code, class_edges);
use = ~isnan(gates.ring) & ~isnan(gates.class) & gates.waveform ~= "";
summary = groupsummary(gates(use, :), ["waveform" "ring" "class"], "sum", ["density" "predicted"]);
summary.lo_km = ring_edges(summary.ring)';
summary.class_name = class_names(summary.class)';
summary.ratio = summary.sum_density ./ summary.sum_predicted;
writetable(summary(:, ["waveform" "lo_km" "class_name" "GroupCount" "sum_density" "sum_predicted" "ratio"]), ...
    fullfile(figure_dir, "NR_C_cfp.csv"))
% Share of the gates of each cut type and ring in each class
group = findgroups(summary(:, ["waveform" "ring"]));
group_total = splitapply(@sum, summary.GroupCount, group);
summary.share = summary.GroupCount ./ group_total(group);
share = unstack(summary(:, ["waveform" "lo_km" "class_name" "share"]), "share", "class_name", VariableNamingRule="preserve");
disp(share)
writetable(share, fullfile(figure_dir, "NR_C_cfp_share.csv"))
for waveform = ["split" "batch" "staggered"]
    t = summary(summary.waveform == waveform, :);
    if isempty(t), continue, end
    t.ratio = round(t.ratio, 2);
    fprintf("%s cuts: density / predicted by ring (rows) and filter class (columns)\n", waveform)
    disp(unstack(t(:, ["lo_km" "class_name" "ratio"]), "ratio", "class_name", VariableNamingRule="preserve"))
end
