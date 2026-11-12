function out_files = nr_gap_fill(night_names, region_name, cut_set, fixed_range_km, keep_cases, max_frames)
%NR_GAP_FILL Drop rules solved on regional cases of shared nights (dev).
%   OUT_FILES = NR_GAP_FILL(NIGHT_NAMES, REGION_NAME, CUT_SET,
%   FIXED_RANGE_KM, KEEP_CASES) builds, for each night of the shared night
%   set and the region REGION_NAME ("midwest": the 28 sandbox radars,
%   [-100 -85 36 46]; "east": Appalachians to the coast, [-85 -74 34 42]),
%   a regional case from the cleaned volumes (nr_night_case,
%   bmm4d.build_observations with its 6-worker pool), then solves it with
%   four rules and scores each as NR_D:
%     all rows; CUT_SET within FIXED_RANGE_KM; CUT_SET within this night's
%     adaptive R_min per radar (nr_radar_night on
%     products/near_radar/NR_E/rows_<night>.mat, from nr_diagnose_nights);
%     the whole CUT_SET.
%   CUT_SET is "batch" or ["batch" "staggered"]. The split entries within
%   20 km of every second radar are held out of every solve. Results:
%   products/near_radar/NR_E/gap_<region>_<night>.mat (summary, distance);
%   a night whose file exists is skipped. KEEP_CASES false (default)
%   deletes each case.mat (~1-3 GB) once scored. MAX_FRAMES (default Inf)
%   keeps the first frames only (a pilot: case nr_..._pilot, results
%   gap_<region>_<night>_pilot.mat, adaptive R_min from the pilot rows if
%   there are no full ones). Meant for a cluster job; ~10 min per night
%   and region.

arguments
    night_names (1, :) string
    region_name (1, 1) string
    cut_set (1, :) string = "batch"
    fixed_range_km (1, 1) double = 20
    keep_cases (1, 1) logical = false
    max_frames (1, 1) double = inf
end
suffix = "";
if isfinite(max_frames)
    suffix = "_pilot";
end
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"), ...
    fileparts(mfilename("fullpath")))
regions = struct("midwest", [-100 -85 36 46], "east", [-85 -74 34 42]);
products_dir = nr_shared_path("/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/BMM-4D/products");
out_dir = fullfile(products_dir, "near_radar", "NR_E");
rule_cuts = "batch";
if isequal(sort(cut_set), ["batch" "staggered"])
    rule_cuts = "batch_staggered";
end
rules = struct("name", {"all rows", "fixed", "adaptive", "whole cut set"}, ...
    "cuts", {"none", rule_cuts, rule_cuts, rule_cuts}, "max_range_km", {0, fixed_range_km, NaN, Inf});

out_files = strings(size(night_names));
for i_n = 1:numel(night_names)
    night = night_names(i_n);
    out_files(i_n) = fullfile(out_dir, "gap_" + region_name + "_" + night + suffix + ".mat");
    if isfile(out_files(i_n))
        continue
    end
    case_tic = tic;

    % Adaptive R_min of this night's radars, from the per-volume diagnostic
    rows_file = fullfile(out_dir, "rows_" + night + ".mat");
    if ~isfile(rows_file)
        rows_file = fullfile(out_dir, "rows_" + night + "_pilot.mat");
    end
    diagnostic = load(rows_file, "rows");
    diagnostic.rows.night = dateshift(diagnostic.rows.datetime - hours(16), "start", "day");
    radar_night = nr_radar_night(diagnostic.rows, cut_set);
    radar_night = radar_night(isfinite(radar_night.r_min_km), :);

    % Regional case
    c = nr_night_case("nr_" + region_name + "_" + extractAfter(night, "night_us_") + suffix, ...
        fullfile(products_dir, "prod", night), regions.(region_name), max_frames);
    if ~isfile(fullfile(c.data_dir, "case.mat"))
        bmm4d.build_observations(c);
    end
    S = ar_load(c);
    manifest = load(fullfile(c.data_dir, "manifest.mat"), "radars", "raw_files");
    cut_code = nr_cut_types(S.observations, S.radar_selected, manifest);
    observed = S.observations.values ~= 0;
    range_km = double(S.observations.row_range_m) / 1000;
    holdout = observed & cut_code == 1 & range_km < 20 & mod(double(S.observations.row_radar_id), 2) == 0;
    [found, id] = ismember(S.radar_selected.name, radar_night.station);
    adaptive_km = repmat(median(radar_night.r_min_km), height(S.radar_selected), 1);
    adaptive_km(found) = radar_night.r_min_km(id(found));

    summary = table();
    distance = table();
    for rule = rules
        if isnan(rule.max_range_km)
            rule.max_range_km = adaptive_km;
        end
        drop = nr_drop_rows(S.observations, cut_code, rule);
        result = ar_solve(nr_apply_drop(S, drop | holdout), night + " " + region_name + ": " + rule.name);
        score = nr_score(S, result.field, cut_code, drop);
        predicted = S.observations.F(:, S.active_id) * double(result.field);
        values = full(S.observations.values);
        ok = holdout & values > 0.5 & predicted > 0.5;
        dropped = score.residuals(score.residuals.dropped, :);
        kept_batch = score.residuals(~score.residuals.dropped & score.residuals.cuts == "batch" & ...
            score.residuals.lo_km >= 20 & score.residuals.hi_km <= 40, :);
        summary = [summary; table(region_name, night, string(rule.name), nnz(drop) / nnz(observed), ...
            sum(dropped.observed) / sum(dropped.predicted), sum(values(ok)) / sum(predicted(ok)), ...
            sum(kept_batch.observed) / sum(kept_batch.predicted), score.misfit_trusted, ...
            score.distance.column_ratio(1), score.distance.column_ratio(2), score.distance.column_ratio(3), ...
            result.seconds, VariableNames=["region" "night" "rule" "dropped_share" "dropped_ratio" ...
            "holdout_ratio" "kept_batch_20_40" "misfit_trusted" "column_0_10" "column_10_20" "column_20_30" ...
            "seconds"])]; %#ok<AGROW>
        distance = [distance; [table(repmat(region_name, height(score.distance), 1), ...
            repmat(night, height(score.distance), 1), repmat(string(rule.name), height(score.distance), 1), ...
            VariableNames=["region" "night" "rule"]) score.distance]]; %#ok<AGROW>
        clear result predicted values
    end
    clear S cut_code observed holdout
    save(out_files(i_n), "summary", "distance", "adaptive_km")
    if ~keep_cases
        delete(fullfile(c.data_dir, "case.mat"))   % rebuilt by build_observations if needed again
    end
    disp(summary)
    fprintf("[near_radar] %s %s: %.0f s\n", night, region_name, toc(case_tic))
end
end
