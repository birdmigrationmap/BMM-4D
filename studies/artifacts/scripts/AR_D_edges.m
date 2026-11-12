%% Artefacts D: radar coverage edges in the column maps (dev, fix 3)
% Screenshots 1 and 3 of the us_east video: discs with sharp circular
% edges around KLNX, KGLD, KUDX (High Plains, the -104 deg domain edge)
% and arcs between radars (KEAX, KSGF). The radar-profile prior does not
% honour each radar's own level (its 150 km Gaussian blend mixes in the
% neighbours: at KLNX the prior aloft is 20-40 birds km^-3 where KLNX
% reads 0-5), and the prior term is white: a cell's deviation from the
% prior is not correlated with its neighbours', and the Laplacian
% smoothness beside it is weak (deviations fade within about one 4 km
% cell). So where a radar's rows stop, the field jumps back to the prior.
% Solves of ar_plains (12 radars, 04:30-07:30 UTC):
% - production: white prior;
% - prior 75 km: Gaussian radar weights of 75 km instead of 150 km;
% - smoothed prior 75 km: a second solve whose prior carries the first
%   solve's local deviation past the edge of the data, fading over 75 km
%   (ar_smooth_prior), warm-started;
% - GP 50 / 100 / 150 km: the correlated prior (solve_advection_lean
%   prior_correlation_m, same marginal variance as the white prior), so a
%   deviation set by the data fades over that length.
% Scores: the step of the log column where the data end against the same
% step 40 km inside the data (ar_edge_score), misfit, time; maps at
% 04:30, 06:00, 07:30 (production and GP 100 km).
% Outputs: figures/AR_D_edges.png, AR_D_edges.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_plains", [-104 -94 37.5 47.5], times);
S = ar_load(c);
inv = c.inversion;

%% Solves ------------------------------------------------------------------------
results = {ar_solve(S, "production (white prior)")};
prior_75 = ar_profile_prior(S.g, S.radar_scans, S.radar_profiles, S.active_mask_3d, 75000, "asl");
results{end + 1} = ar_solve(S, "prior 75 km", prior=prior_75);
prior_smooth = ar_smooth_prior(S, S.prior, results{1}.field, 75000);
results{end + 1} = ar_solve(S, "smoothed prior 75 km", prior=prior_smooth, advection_reference=S.prior, ...
    x_initial=results{1}.field);
results{end}.seconds = results{end}.seconds + results{1}.seconds;   % two solves
for length_km = [50 100 150]
    results{end + 1} = ar_solve(S, sprintf("GP %d km", length_km), prior_correlation_m=length_km * 1000); %#ok<SAGROW>
end

%% Scores ------------------------------------------------------------------------
scores = table();
for i_r = 1:numel(results)
    score = ar_edge_score(S, results{i_r}.column);
    scores = [scores; table(results{i_r}.name, score.jump, score.inner, score.edge, results{i_r}.misfit, ...
        results{i_r}.misfit_trusted, results{i_r}.seconds, results{i_r}.iterations, ...
        VariableNames=["solve" "step_at_edge" "step_inside" "edge" "misfit" "misfit_30_120km" "seconds" "iterations"])]; %#ok<AGROW>
end
disp(scores)
writetable(scores, fullfile(figure_dir, "AR_D_edges.csv"))

%% Maps: production and GP 100 km --------------------------------------------------
ar_plot_maps(S, results([1 5]), [1 4 7], fullfile(figure_dir, "AR_D_edges.png"), bounds=[-104.5 -97 38.5 45.5]);
