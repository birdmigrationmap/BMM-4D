function [radar_night, deficit, deficit_se, ring_edges] = nr_radar_night(rows, cuts, tolerance)
%NR_RADAR_NIGHT Near-radar deficit and adaptive R_min per radar-night (dev).
%   [RADAR_NIGHT, DEFICIT, DEFICIT_SE, RING_EDGES] = NR_RADAR_NIGHT(ROWS)
%   pools the nr_volume_ratio rows of many volumes (with columns set and
%   night added) by radar-night, for the cuts CUTS (default "batch"; a
%   list of waveform names). Two levels are kept apart:
%   - far_level: these cuts at 25-50 km over the split cuts at 10-40 km of
%     the same radar-night, a cut-type offset that does not depend on range
%     (not what dropping near rows fixes);
%   - DEFICIT [n_radar_night n_ring]: each ring's ratio over that
%     radar-night's own level of these cuts at 25-50 km (ratio of sums over
%     volumes), DEFICIT_SE the jackknife (over volumes) standard error of
%     its log.
%   A ring within 25 km is low when its deficit is below 1 - TOLERANCE
%   (default 0.05) by more than two standard errors; r_min_km is the outer
%   edge of the outermost low ring (0: none; NaN with fewer than 3
%   volumes). r_min_significant_km is the same with no tolerance: it grows
%   with the amount of data, since any small deficit becomes significant.
%   RADAR_NIGHT has set, station, night, vcp, n_volume, has_cuts,
%   far_level, deficit_<ring>, r_min_km, r_min_significant_km.

arguments
    rows table
    cuts (1, :) string = "batch"
    tolerance (1, 1) double = 0.05
end

volume_key = rows.station + "_" + string(rows.datetime, "yyyyMMdd_HHmmss");
ring_edges = unique([rows.range_lo_km; rows.range_hi_km])';
[radar_night_id, radar_night] = findgroups(rows(:, ["set" "station" "night"]));
n_rn = height(radar_night);
n_ring = numel(ring_edges) - 1;
far_rings = ring_edges(1:end-1) >= 25 & ring_edges(2:end) <= 50;
[deficit, deficit_se] = deal(nan(n_rn, n_ring));
[far_level, r_min, r_min_significant] = deal(nan(n_rn, 1));
[n_volume, has_cuts] = deal(zeros(n_rn, 1));
vcps = strings(n_rn, 1);
for i_rn = 1:n_rn
    mine = rows(radar_night_id == i_rn, :);
    vcps(i_rn) = strjoin(string(unique(mine.vcp))', "+");
    volume_id = findgroups(volume_key(radar_night_id == i_rn));
    n_volume(i_rn) = max(volume_id);
    is_cut = ismember(mine.waveform, cuts);
    has_cuts(i_rn) = any(is_cut);
    if ~any(is_cut), continue, end
    % Sums per volume and ring [n_volume n_ring]
    ring = discretize(mine.range_lo_km, ring_edges);
    D = accumarray([volume_id(is_cut) ring(is_cut)], mine.density(is_cut), [n_volume(i_rn) n_ring]);
    P = accumarray([volume_id(is_cut) ring(is_cut)], mine.predicted(is_cut), [n_volume(i_rn) n_ring]);
    split_ctl = mine.split & mine.range_lo_km >= 10 & mine.range_hi_km <= 40;
    far_cut = sum(D(:, far_rings), "all") / sum(P(:, far_rings), "all");
    far_level(i_rn) = far_cut / (sum(mine.density(split_ctl)) / sum(mine.predicted(split_ctl)));
    deficit(i_rn, :) = sum(D, 1) ./ sum(P, 1) / far_cut;
    if n_volume(i_rn) < 3 || ~isfinite(far_cut) || far_cut <= 0, continue, end
    % Jackknife over volumes of the log deficit (far level recomputed too)
    leave_out = nan(n_volume(i_rn), n_ring);
    for v = 1:n_volume(i_rn)
        D_v = sum(D, 1) - D(v, :);
        P_v = sum(P, 1) - P(v, :);
        leave_out(v, :) = log(D_v ./ P_v) - log(sum(D_v(far_rings)) / sum(P_v(far_rings)));
    end
    n_v = n_volume(i_rn);
    deficit_se(i_rn, :) = sqrt((n_v - 1) / n_v * sum((leave_out - mean(leave_out, 1, "omitnan")) .^ 2, 1, "omitnan"));
    within = ring_edges(2:end) <= 25;
    low = log(deficit(i_rn, :)) + 2 * deficit_se(i_rn, :) < log(1 - tolerance) & within;
    significant = log(deficit(i_rn, :)) + 2 * deficit_se(i_rn, :) < 0 & within;
    r_min(i_rn) = max([0 ring_edges(find(low) + 1)]);
    r_min_significant(i_rn) = max([0 ring_edges(find(significant) + 1)]);
end
radar_night.vcp = vcps;
radar_night.n_volume = n_volume;
radar_night.has_cuts = has_cuts;
radar_night.far_level = far_level;
ring_names = compose("deficit_%g_%g_km", ring_edges(1:end-1)', ring_edges(2:end)');
radar_night = [radar_night array2table(deficit, VariableNames=ring_names)];
radar_night.r_min_km = r_min;
radar_night.r_min_significant_km = r_min_significant;
end
