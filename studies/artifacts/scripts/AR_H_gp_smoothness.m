%% Artefacts H: does the correlated prior need the smoothness term? (dev)
% The correlated (GP) prior of fix 3 penalizes (x - prior)' Q (x - prior)
% with Q = I + (l/dx)^2 L_h: a horizontal Laplacian of the deviation from
% the prior. The production smoothness term is a Laplacian of the field
% itself (horizontal and vertical). As the prior mean is smooth, the two
% overlap horizontally; the GP has no vertical part. Solves of
% ar_plains_bb (12 radars, 04:30-07:30 UTC, beam-blockage fix):
% - white prior, smoothness 0.3 (production);
% - GP 100 km with smoothness 0.3 (fix 3) and 0;
% - GP 100 km horizontally and 300 m vertically (3D), smoothness 0.
% Scores: edge (ar_edge_score), misfit, iterations, roughness (mean
% |second difference| along z and along x over the mean density), and the
% prediction error of held-out radars (ar_holdout_score), the test that the
% data of the other radars decide.
% Outputs: figures/AR_H_gp_smoothness.csv, AR_H_holdout.csv,
% AR_H_gp_smoothness.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_plains_bb", [-104 -94 37.5 47.5], times);
S = ar_load(c);
variants = struct("name", {"white, smoothness 0.3 (production)", "GP 100 km, smoothness 0.3", ...
    "GP 100 km, smoothness 0", "GP 100 km / 300 m, smoothness 0"}, ...
    "correlation_m", {0, 100000, 100000, 100000}, "vertical_m", {0, 0, 0, 300}, "smoothness", {0.3, 0.3, 0, 0});

%% Solves ------------------------------------------------------------------------
results = cell(numel(variants), 1);
scores = table();
for i_v = 1:numel(variants)
    v = variants(i_v);
    results{i_v} = ar_solve(S, v.name, prior_correlation_m=v.correlation_m, prior_correlation_vertical_m=v.vertical_m, ...
        smoothness=v.smoothness);
    score = ar_edge_score(S, results{i_v}.column);
    % Roughness along z and x of the field (frames stacked), over its mean
    field = nan([S.g.sz numel(times)]);
    for i_t = 1:numel(times)
        frame = nan(S.g.sz);
        frame(S.active_id) = results{i_v}.field(:, i_t);
        field(:, :, :, i_t) = frame;
    end
    level = mean(field, "all", "omitnan");
    rough_z = mean(abs(diff(field, 2, 3)), "all", "omitnan") / level;
    rough_x = mean(abs(diff(field, 2, 1)), "all", "omitnan") / level;
    scores = [scores; table(string(v.name), score.edge, results{i_v}.misfit, results{i_v}.misfit_trusted, ...
        rough_z, rough_x, results{i_v}.seconds, results{i_v}.iterations, ...
        VariableNames=["solve" "edge" "misfit" "misfit_30_120km" "roughness_z" "roughness_x" "seconds" "iterations"])]; %#ok<AGROW>
end
clear field frame
disp(scores)
writetable(scores, fullfile(figure_dir, "AR_H_gp_smoothness.csv"))

%% Held-out radars ---------------------------------------------------------------------
% Radars with neighbours on several sides (an edge radar left out is
% extrapolation, not interpolation); the prior is rebuilt without the
% left-out radar's profile (production prior at the production length scale)
held_out = intersect(["KLNX" "KUEX" "KGLD" "KOAX" "KTWX"], string(S.radar_selected.name), "stable");
holdout = table();
for radar_name = held_out(:)'
    for i_v = 1:numel(variants)
        v = variants(i_v);
        score = ar_holdout_score(S, radar_name, @(S_out, keep) ar_solve(S_out, v.name + " without " + radar_name, ...
            prior=ar_profile_prior(S.g, S.radar_scans(keep, :), S.radar_profiles(:, keep, :), S.active_mask_3d, ...
            S.c.inversion.prior_length_scale_m, "asl"), prior_correlation_m=v.correlation_m, ...
            prior_correlation_vertical_m=v.vertical_m, smoothness=v.smoothness));
        holdout = [holdout; table(radar_name, string(v.name), score.relative_error, score.median_abs_log_ratio, ...
            VariableNames=["radar" "solve" "relative_error" "median_abs_log_ratio"])]; %#ok<AGROW>
    end
end
disp(holdout)
disp(groupsummary(holdout, "solve", "mean", ["relative_error" "median_abs_log_ratio"]))
writetable(holdout, fullfile(figure_dir, "AR_H_holdout.csv"))

%% Maps ----------------------------------------------------------------------------
ar_plot_maps(S, results, [1 4 7], fullfile(figure_dir, "AR_H_gp_smoothness.png"), bounds=[-104.5 -97 38.5 45.5]);
