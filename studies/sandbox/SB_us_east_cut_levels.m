%% Sandbox: cut levels on the us_east night (dev)
% The production pipeline for us_east (101 radars east of the High Plains,
% night of 21-22 October 2026, 26 frames) up to the observations, from the
% raw files already listed in prod/us_east/manifest.mat (no download), then
% SB_cut_levels with case_name "us_east". The previous run (23 October 2026,
% before the estimator) is kept in prod/us_east_20261023. The log of the
% whole run goes to prod/us_east/us_east_cut_levels.log.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
case_name = "us_east";
c = case_settings(case_name);
diary(fullfile(c.data_dir, "us_east_cut_levels.log"))

%% Cleaning and observations (production) -------------------------------------
run_tic = tic;
bmm4d.clean_volumes(c);
fprintf("clean_volumes done at %.0f s\n", toc(run_tic))
bmm4d.static_clutter(c);
fprintf("static_clutter done at %.0f s\n", toc(run_tic))
bmm4d.build_observations(c);
fprintf("build_observations done at %.0f s\n", toc(run_tic))

%% Cut levels -------------------------------------------------------------------
run(fullfile(project_root, "studies", "sandbox", "SB_cut_levels.m"))
fprintf("all done at %.0f s\n", toc(run_tic))
diary off
