%% Sandbox: why the batch cuts read low near the radar
% In the sandbox reconstruction the cuts at 1.8 deg and above read 0.6-0.8
% of a solve without near rows at 10-20 km, the split cuts (0.5-1.3 deg)
% 1.1-1.3. This compares cuts at the same height: the mean linear Z of an
% upper cut at range r against that of a lower cut over the ranges where
% its beam centre spans the same heights (about 1.4 r for the pairs below).
% Pairs: 1.8/1.3 (batch/split), 1.3/0.9 (split/split, control), 2.4/1.8
% (batch/batch, control). Split cuts use the surveillance pass (CS).
% Versions: the radar output (REF, below threshold = 0, range folded
% missing); the output with the clutter filter's removed power put back
% (CFP codes 8-18: up to 10 dB; gates with more removed, likely ground
% clutter, and point-clutter gates (code 1, value replaced) left out, on
% both cuts); and the estimator (cleaned volumes, density).
% Also the share of gates the clutter filter leaves alone (CFP code 0)
% and the share with >= 3 dB removed, by cut type and range.
% Outputs: figures/cut_type_deficit.png, cut_type_deficit.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
load(fullfile(c.data_dir, "case.mat"), "radar_selected")
volume_dir = fullfile(fileparts(c.data_dir), "volumes", "9c2b4ef2");
effective_radius = 6371200 * 4 / 3;
beam_height = @(r, e) sqrt(r .^ 2 + effective_radius ^ 2 + 2 * r * effective_radius .* sind(e)) - effective_radius;
time_ids = [1 3 6];
range_edges_km = [6 9 12 16 20 25 30 40 55];
share_edges_km = [2.5 5 10 15 20 30 40 60 80];
pairs = table(["batch/split (1.8/1.3)"; "split/split (1.3/0.9)"; "batch/batch (2.4/1.8)"], ...
    [1.6 2.0; 1.1 1.4; 2.2 2.6], [1.1 1.4; 0.7 1.0; 1.6 2.0], VariableNames=["name" "upper" "lower"]);
versions = ["radar output" "removed power put back" "estimator"];

%% Same-height ratios, every radar and time -----------------------------------
n_radar = height(radar_selected);
ratio = nan(n_radar, numel(time_ids), numel(range_edges_km) - 1, height(pairs), numel(versions));
share_unfiltered = nan(n_radar, numel(time_ids), numel(share_edges_km) - 1, 2);   % split, batch
share_removed = nan(size(share_unfiltered));
for i_r = 1:n_radar
    files = dir(fullfile(volume_dir, radar_selected.name(i_r) + "*.mat"));
    for i_t = 1:numel(time_ids)
        if numel(files) < time_ids(i_t), continue, end
        volume = load(fullfile(files(time_ids(i_t)).folder, files(time_ids(i_t)).name)).volume;
        cuts = bmm4d.radar.read_level2_moments(volume.source_file, ["REF" "CFP"], 2.7);
        cuts = cuts(arrayfun(@(x) ~isempty(x.CFP), cuts));   % CS of split cuts, and batch cuts
        elev = [cuts.elev];
        % Radar output and output with the removed power, as mean linear Z
        z_out = cell(1, numel(cuts));
        z_back = cell(1, numel(cuts));
        for k = 1:numel(cuts)
            code = double(cuts(k).CFP_code);
            z = 10 .^ (cuts(k).REF / 10);
            z(cuts(k).REF_code == 0) = 0;
            z(cuts(k).REF_code == 1) = nan;
            removed_db = max(code - 8, 0);
            z_out{k} = z;
            z_back{k} = z .* 10 .^ (removed_db / 10);
            z_back{k}(removed_db > 10 | code == 1) = nan;
        end
        for i_p = 1:height(pairs)
            upper = find(elev > pairs.upper(i_p, 1) & elev < pairs.upper(i_p, 2), 1);
            lower = find(elev > pairs.lower(i_p, 1) & elev < pairs.lower(i_p, 2), 1);
            if isempty(upper) || isempty(lower), continue, end
            ratio(i_r, i_t, :, i_p, 1) = same_height_ratio(z_out{upper}, z_out{lower}, cuts(upper), cuts(lower), range_edges_km, beam_height);
            ratio(i_r, i_t, :, i_p, 2) = same_height_ratio(z_back{upper}, z_back{lower}, cuts(upper), cuts(lower), range_edges_km, beam_height);
            % Estimator: the cleaned sweeps of the same elevations
            sweeps = volume.dens.sweeps;
            sweep_elev = [sweeps.elev];
            upper = find(sweep_elev > pairs.upper(i_p, 1) & sweep_elev < pairs.upper(i_p, 2), 1);
            lower = find(sweep_elev > pairs.lower(i_p, 1) & sweep_elev < pairs.lower(i_p, 2), 1);
            ratio(i_r, i_t, :, i_p, 3) = same_height_ratio(sweeps(upper).data, sweeps(lower).data, ...
                struct("REF_range", sweeps(upper).range, "elev", sweeps(upper).elev), ...
                struct("REF_range", sweeps(lower).range, "elev", sweeps(lower).elev), range_edges_km, beam_height);
        end
        % Clutter filter: share of gates left alone and with >= 3 dB removed
        for i_c = 1:2
            if i_c == 1
                k = find(elev > 1.1 & elev < 1.4, 1);
            else
                k = find(elev > 1.6 & elev < 2.0, 1);
            end
            if isempty(k), continue, end
            code = double(cuts(k).CFP_code);
            range_km = double(cuts(k).REF_range(:)) / 1000;
            for i_b = 1:numel(share_edges_km) - 1
                in_bin = range_km >= share_edges_km(i_b) & range_km < share_edges_km(i_b + 1);
                share_unfiltered(i_r, i_t, i_b, i_c) = mean(code(in_bin, :) == 0, "all");
                share_removed(i_r, i_t, i_b, i_c) = mean(code(in_bin, :) >= 11, "all");
            end
        end
    end
    fprintf("%s done\n", radar_selected.name(i_r))
end

%% Summary table ------------------------------------------------------------------
range_mid_km = (range_edges_km(1:end - 1) + range_edges_km(2:end)) / 2;
summary = table();
for i_v = 1:numel(versions)
    for i_p = 1:height(pairs)
        values = reshape(ratio(:, :, :, i_p, i_v), [], numel(range_mid_km));
        summary = [summary; table(repmat(versions(i_v), numel(range_mid_km), 1), repmat(pairs.name(i_p), numel(range_mid_km), 1), ...
            range_mid_km', median(values, 1, "omitnan")', prctile(values, 25, 1)', prctile(values, 75, 1)', ...
            VariableNames=["version" "pair" "upper_range_km" "median" "q25" "q75"])]; %#ok<AGROW>
    end
end
writetable(summary, fullfile(figure_dir, "cut_type_deficit.csv"))
disp(unstack(summary(:, ["version" "pair" "upper_range_km" "median"]), "median", "upper_range_km"))

%% Plot results -------------------------------------------------------------------
pair_colors = [0.80 0.33 0.10; 0.25 0.25 0.25; 0.20 0.45 0.75];
fig = bmm4d.plot.create_export_figure(3000, 1000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(versions) + 1, TileSpacing="tight", Padding="tight");
for i_v = 1:numel(versions)
    ax = nexttile(t);
    hold(ax, "on")
    for i_p = 1:height(pairs)
        rows = summary.version == versions(i_v) & summary.pair == pairs.name(i_p);
        fill(ax, [range_mid_km fliplr(range_mid_km)], [summary.q25(rows)' fliplr(summary.q75(rows)')], pair_colors(i_p, :), ...
            FaceAlpha=0.15, EdgeColor="none", HandleVisibility="off")
        plot(ax, range_mid_km, summary.median(rows), "-o", Color=pair_colors(i_p, :), LineWidth=2, MarkerFaceColor=pair_colors(i_p, :), DisplayName=pairs.name(i_p))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, FontSize=12, YLim=[0.3 1.2], XLim=[5 50])
    grid(ax, "on")
    title(ax, versions(i_v), FontSize=14)
    xlabel(ax, "Range of the upper cut (km)")
    if i_v == 1
        ylabel(ax, "Upper / lower cut at the same height")
        legend(ax, Location="northwest", FontSize=11)
    end
end
ax = nexttile(t);
hold(ax, "on")
share_mid_km = (share_edges_km(1:end - 1) + share_edges_km(2:end)) / 2;
cut_names = ["split 1.3°" "batch 1.8°"];
cut_colors = [0.25 0.25 0.25; 0.80 0.33 0.10];
for i_c = 1:2
    plot(ax, share_mid_km, squeeze(median(reshape(share_unfiltered(:, :, :, i_c), [], numel(share_mid_km)), 1, "omitnan")), "-o", ...
        Color=cut_colors(i_c, :), LineWidth=2, MarkerFaceColor=cut_colors(i_c, :), DisplayName=cut_names(i_c) + ", not filtered")
    plot(ax, share_mid_km, squeeze(median(reshape(share_removed(:, :, :, i_c), [], numel(share_mid_km)), 1, "omitnan")), "--s", ...
        Color=cut_colors(i_c, :), LineWidth=2, DisplayName=cut_names(i_c) + ", \geq3 dB removed")
end
set(ax, FontSize=12, YLim=[0 1], XLim=[0 75])
grid(ax, "on")
title(ax, "Clutter filter", FontSize=14)
xlabel(ax, "Range (km)")
% Labels on the lines instead of a legend, which does not fit beside the panel
text(ax, 72, 0.88, "split 1.3°, not filtered", Color=cut_colors(1, :), FontSize=11, HorizontalAlignment="right")
text(ax, 72, 0.47, "batch 1.8°, not filtered", Color=cut_colors(2, :), FontSize=11, HorizontalAlignment="right")
text(ax, 3, 0.97, "dashed: \geq3 dB removed", Color=[0.3 0.3 0.3], FontSize=11)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "cut_type_deficit.png"), Resolution=300)
close(fig)

%% Local functions ----------------------------------------------------------
function ratio = same_height_ratio(z_upper, z_lower, upper, lower, range_edges_km, beam_height)
% Mean of the upper cut per range bin over the mean of the lower cut where
% its beam centre spans the same heights
range_upper = double(upper.REF_range(:));
range_lower = double(lower.REF_range(:));
height_upper = beam_height(range_upper, upper.elev);
height_lower = beam_height(range_lower, lower.elev);
ratio = nan(1, numel(range_edges_km) - 1);
for i_b = 1:numel(range_edges_km) - 1
    in_bin = range_upper >= range_edges_km(i_b) * 1000 & range_upper < range_edges_km(i_b + 1) * 1000;
    same_height = height_lower >= min(height_upper(in_bin)) & height_lower < max(height_upper(in_bin));
    ratio(i_b) = mean(z_upper(in_bin, :), "all", "omitnan") / mean(z_lower(same_height, :), "all", "omitnan");
end
end
