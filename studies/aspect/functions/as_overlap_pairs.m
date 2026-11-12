function pairs = as_overlap_pairs(bins)
%AS_OVERLAP_PAIRS Pairs of radars seeing the same cell at the same time (dev).
%   PAIRS = AS_OVERLAP_PAIRS(BINS) takes the cells of as_frame_bins (one
%   row per radar and cell) and returns one row per pair of radars in the
%   same cell and frame, radar_a < radar_b, with the a / b values of each
%   variable (raw, ell, fix, az, range_m, elev, u, v, sd2, n) and the cell
%   height. Cells seen by m radars give m (m - 1) / 2 pairs.

bins = sortrows(bins, ["frame" "key" "radar"]);
group = [true; diff(double(bins.key)) ~= 0 | diff(double(bins.frame)) ~= 0];
group_id = cumsum(group);
max_m = max(accumarray(group_id, 1));
ia = [];
ib = [];
for offset = 1:max_m - 1
    i = (1:height(bins) - offset)';
    same = group_id(i) == group_id(i + offset);
    ia = [ia; i(same)]; %#ok<AGROW>
    ib = [ib; i(same) + offset]; %#ok<AGROW>
end
pairs = table(bins.frame(ia), bins.key(ia), bins.radar(ia), bins.radar(ib), ...
    (bins.height_asl_m(ia) + bins.height_asl_m(ib)) / 2, VariableNames=["frame" "key" "radar_a" "radar_b" "height_asl_m"]);
for name = setdiff(string(bins.Properties.VariableNames), ["frame" "key" "radar" "height_asl_m"], "stable")
    pairs.(name + "_a") = bins.(name)(ia);
    pairs.(name + "_b") = bins.(name)(ib);
end
end
