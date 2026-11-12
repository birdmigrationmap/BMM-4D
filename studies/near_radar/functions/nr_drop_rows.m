function drop = nr_drop_rows(observations, cut_code, rule)
%NR_DROP_ROWS Observation entries a near-radar rule removes (dev).
%   DROP = NR_DROP_ROWS(OBSERVATIONS, CUT_CODE, RULE) is logical
%   [n_row n_time], true for the observed entries (values ~= 0) the rule
%   removes. CUT_CODE is nr_cut_types. RULE has:
%     cuts      - which cuts: "none" (drop nothing), "elevation"
%                 (row elevation >= 1.6 deg, i.e. 1.8 deg and above),
%                 "batch" (batch waveform from the header),
%                 "batch_staggered" (batch and staggered pulse pair), "not_split"
%                 (every cut whose reflectivity is not from a contiguous
%                 surveillance scan: batch, Doppler, staggered);
%     max_range_km - the cuts are dropped closer than this to their radar:
%                 a scalar, Inf (the whole cut), or a vector with one
%                 value per radar (adaptive, indexed by row_radar_id).

drop = false(size(observations.values));
if rule.cuts == "none"
    return
end
switch rule.cuts
    case "elevation"
        selected = repmat(observations.row_elev_deg(:) >= 1.6, 1, size(drop, 2));
    case "batch"
        selected = cut_code == 4;
    case "batch_staggered"
        selected = cut_code == 4 | cut_code == 5;
    case "not_split"
        selected = cut_code >= 2;
    otherwise
        error("nr_drop_rows:UnknownCuts", "Unknown cut set: %s", rule.cuts)
end
max_range_m = 1000 * rule.max_range_km(:);
if isscalar(max_range_m)
    near = observations.row_range_m(:) < max_range_m;
else
    near = observations.row_range_m(:) < max_range_m(observations.row_radar_id);
end
drop = selected & near & (observations.values ~= 0);
end
