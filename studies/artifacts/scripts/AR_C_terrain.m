%% Artefacts C: column deficit over the Appalachians (dev, fix 2)
% Screenshot 2 of the us_east video: a dark band along the Blue Ridge and
% the Smokies (column over ridges 0.38 of the valleys) and streaks
% downwind. Two causes:
% - the field is layered in height above sea level: the prior places each
%   radar's profile at the same sea-level heights everywhere, so a ridge
%   cuts off the dense low layer (prior alone: ridge / valley 0.40), and
%   the advection links keep the sea-level height, so birds carried into a
%   ridge are lost and the lee has no inflow;
% - the beam-blockage bug (README): KMRX and KGSP's blocked low beams
%   over the Smokies are kept and read zero density.
% Solves of ar_appalachians (16 radars, 04:30-07:30 UTC), production
% cleaning, then ar_appalachians_bb (blockage fixed, AR_A_cases.m):
% - production: links and prior at height above sea level (ASL);
% - AGL links: links keep the height above the lowest active cell (ar_steps);
% - AGL prior: radar profiles laid from each column's ground up
%   (ar_profile_prior, "agl");
% - AGL links + prior;
% - blockage fixed, ASL / AGL links + prior.
% Scores: ridge / valley column, column ratio per km of ground
% (ar_terrain_score), data misfit; maps at 04:30, 06:00, 07:30 with the DEM.
% Outputs: figures/AR_C_terrain.png, AR_C_terrain.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
box = [-89 -78 33.5 40];

%% Solves ------------------------------------------------------------------------
results = {};
scores = table();
for case_name = ["ar_appalachians" "ar_appalachians_bb"]
    c = ar_test_case(case_name, box, times);
    S = ar_load(c);
    inv = c.inversion;
    % Terrain-following steps and prior
    steps_agl = ar_steps(S.g, S.radar_scans, seconds(diff(c.times)), S.active_mask_3d, inv, "agl");
    prior_agl = ar_profile_prior(S.g, S.radar_scans, S.radar_profiles, S.active_mask_3d, inv.prior_length_scale_m, "agl");
    if case_name == "ar_appalachians"
        S_map = S;
        solves = {ar_solve(S, "production (ASL)")
            ar_solve(S, "AGL links", steps=steps_agl)
            ar_solve(S, "AGL prior", prior=prior_agl)
            ar_solve(S, "AGL links + prior", steps=steps_agl, prior=prior_agl)};
    else
        solves = {ar_solve(S, "blockage fixed (ASL)")
            ar_solve(S, "blockage fixed + AGL", steps=steps_agl, prior=prior_agl)};
    end
    for i_r = 1:numel(solves)
        score = ar_terrain_score(S, solves{i_r}.column);
        scores = [scores; table(solves{i_r}.name, score.ridge_valley, score.per_km, solves{i_r}.misfit, ...
            solves{i_r}.misfit_trusted, solves{i_r}.seconds, solves{i_r}.iterations, ...
            VariableNames=["solve" "ridge_valley" "per_km" "misfit" "misfit_30_120km" "seconds" "iterations"])]; %#ok<AGROW>
    end
    results = [results; solves]; %#ok<AGROW>
end
disp(scores)
writetable(scores, fullfile(figure_dir, "AR_C_terrain.csv"))

%% Maps: production, AGL, blockage fixed + AGL -----------------------------------
ar_plot_maps(S_map, results([1 4 6]), [1 4 7], fullfile(figure_dir, "AR_C_terrain.png"), ...
    bounds=[-88 -79 34 38.5], extra=max(S_map.g.dem, 0), extra_title="Ground (m)", extra_limits=[0 1600]);
