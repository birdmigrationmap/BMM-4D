%% Artefacts I: the shared aspect function in the inversion (dev, fix 1 follow-up)
% studies/aspect measured the aspect effect where two radars see the same
% cell (gradient-free): side-on 1.29, head-on 1.01, tail-on 0.70 of the
% mean, growing with height, and proposed one shared function instead of
% a per-volume fit (as_apply_aspect). Its scores were on polar volumes;
% here the four cleanings of the Midwest case are solved (production
% settings, white prior, sea-level links) and compared on the
% reconstruction:
% - ar_midwest: elliptical shape correction (round 0, the paper's);
% - ar_midwest_none: no aspect correction;
% - ar_midwest_fix: per-volume symmetric fit (ar_clean_volume);
% - ar_midwest_shared: shared aspect function (studies/aspect).
% The last three share the cleaned volumes (beam-blockage fix) and differ
% only by the correction (AR_A_cases.m, section 5).
% Scores: wedges (ar_wedge_score), harmonics of the observations at
% 10-60 km (ar_obs_harmonics: first, side-on, front/back), misfit, and the
% prediction of held-out radars (ar_holdout_score, prior rebuilt without
% the radar): each cleaning predicts its own data, so a correction that
% removes a bias common to all radars makes the radars agree better.
% Outputs: figures/AR_I_aspect.csv, AR_I_holdout.csv, AR_I_aspect.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
box = [-100 -85 36 46];
no_shape = {"read_options", {"shape", false}};
cleanings = struct("case_name", {"ar_midwest", "ar_midwest_none", "ar_midwest_fix", "ar_midwest_shared"}, ...
    "name", {"ellipse (round 0)", "no correction", "per-volume symmetric", "shared function"}, ...
    "shape", {true, false, false, false});

%% Solves and scores --------------------------------------------------------------
results = cell(1, numel(cleanings));
scores = table();
holdout = table();
for i_c = 1:numel(cleanings)
    cl = cleanings(i_c);
    if cl.shape
        c = ar_test_case(cl.case_name, box, times);
    else
        c = ar_test_case(cl.case_name, box, times, no_shape{:});
    end
    S = ar_load(c);
    results{i_c} = ar_solve(S, cl.name);
    harmonics = ar_obs_harmonics(S);
    wedge = ar_wedge_score(S, results{i_c}.column);
    kdmx = S.radar_selected.name == "KDMX";
    scores = [scores; table(string(cl.name), wedge.median, wedge.share_below_08, min(wedge.per_radar(kdmx, :)), ...
        median(harmonics.first, "all", "omitnan"), median(harmonics.side_on, "all", "omitnan"), ...
        median(harmonics.front_back, "all", "omitnan"), results{i_c}.misfit, results{i_c}.misfit_trusted, ...
        VariableNames=["cleaning" "wedge_median" "wedge_share_below_08" "kdmx_deepest" "first_harmonic" ...
        "side_on_harmonic" "front_back_harmonic" "misfit" "misfit_30_120km"])]; %#ok<AGROW>
    % Held-out radars with neighbours on all sides
    held_out = intersect(["KDMX" "KDVN" "KOAX" "KEAX" "KARX"], string(S.radar_selected.name), "stable");
    for radar_name = held_out(:)'
        score = ar_holdout_score(S, radar_name, @(S_out, keep) ar_solve(S_out, cl.name + " without " + radar_name, ...
            prior=ar_profile_prior(S.g, S.radar_scans(keep, :), S.radar_profiles(:, keep, :), S.active_mask_3d, ...
            S.c.inversion.prior_length_scale_m, "asl")));
        holdout = [holdout; table(radar_name, string(cl.name), score.relative_error, score.median_abs_log_ratio, ...
            VariableNames=["radar" "cleaning" "relative_error" "median_abs_log_ratio"])]; %#ok<AGROW>
    end
    if i_c == 1
        S_map = S;
    end
end
disp(scores)
disp(groupsummary(holdout, "cleaning", "mean", ["relative_error" "median_abs_log_ratio"]))
writetable(scores, fullfile(figure_dir, "AR_I_aspect.csv"))
writetable(holdout, fullfile(figure_dir, "AR_I_holdout.csv"))

%% Maps around KDMX --------------------------------------------------------------
ar_plot_maps(S_map, results, [3 5 7], fullfile(figure_dir, "AR_I_aspect.png"), bounds=[-96.5 -91 39.8 43.8]);
