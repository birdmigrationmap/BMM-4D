%% Regression reference: paper-v1 outputs on a small Key West case
% Three KBYX volumes (23 April 2018, 03:00-03:30 UTC) go through two chains:
% - paper path (regression_case_paper): run on the code at tag paper-v1,
%   checked against today's code, which must give the same result;
% - production path (regression_case_prod): run on today's code (it does
%   not exist at paper-v1).
% The outputs are saved to data/regression/paper_v1_reference.mat (local,
% not tracked) and compared by tests/paperRegressionTest.m. Rerun this
% script only to rebuild the reference, never to accept a change.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"), fullfile(project_root, "tests", "regression"))
reference_dir = fullfile(project_root, "data", "regression");
if ~isfolder(reference_dir)
    mkdir(reference_dir)
end
kernel_dir = fullfile(project_root, "data", "derived", "F_elev", "150km");

%% Case and settings ----------------------------------------------------------
times = datetime("2018-04-23 03:00") + minutes([0 15 30]);
raw_files = load_rd("KBYX", times, closest_accepted=15);
assert(all(isfile(raw_files)), "Missing KBYX raw volumes")

% Paper settings (09_key_west S09_A and S09_C), without the manual removal of
% the north-sector rain artefact; fewer iterations to keep the test short.
paper_settings = struct( ...
    "read_options", {{"removeRain", "rho", "limitElev", 7, "minRange", 3000, "maxRange", 160000, ...
        "maxAlt", 3500, "removeOutliers", 2, "shape", true, "smoothing", 2000, "data_res", [1000 1000], ...
        "waitbar", false, "parallel", false}}, ...
    "grid_bounds", [-83.2, -80.2, 22.7, 25.6], ...
    "heightmax", 3000, ...
    "min_row_mass", 0.75, ...
    "water_override", [-82.05 -81.98 -81.93 -81.97 -82.04; 24.02 24.07 24.02 23.97 23.98], ...
    "inversion", struct("lambda_profile", 10, "lambda_time_profile", 1, "density_prior_idw_neighbors", 3, ...
        "density_prior_idw_power", 2.25, "velocity_idw_neighbors", 1, "velocity_idw_power", 1.6, ...
        "prior", 0.003, "smoothness", 1, "advection", 3, "coeff_horver", -0.5, ...
        "shoreline_crossing_coeff", 0.3, "maxit", 50));

% Production settings as of today, frozen in the reference file
prod_settings = prod_case("keywest");
prod_settings.inversion.maxit = 50;
% Plain PCG: the projected FISTA step (nonnegative) is not reproducible from
% run to run (0.1-0.5% relative difference in the field, single-precision
% parallel sums), while PCG alone is identical.
prod_settings.inversion.nonnegative = false;

%% Paper path on the paper-v1 code ------------------------------------------------
% A detached worktree of paper-v1 replaces functions/ on the path; the
% untracked regtools (get_l) and the data folder come from this checkout.
worktree_dir = fullfile(tempdir, "bmm4d_paper_v1");
if isfolder(worktree_dir)
    system(sprintf('git -C "%s" worktree remove --force "%s"', project_root, worktree_dir));
end
assert(system(sprintf('git -C "%s" worktree add --detach "%s" paper-v1', project_root, worktree_dir)) == 0)
current_path = path;
rmpath(genpath(fullfile(project_root, "functions")), fullfile(project_root, "prod", "functions"))
addpath(genpath(fullfile(worktree_dir, "functions")), genpath(fullfile(project_root, "functions", "regtools")))
try
    assert(startsWith(which("read_rd"), worktree_dir))
    paper_v1 = regression_case_paper(raw_files, kernel_dir, paper_settings);
catch err
    path(current_path)
    rethrow(err)
end
path(current_path)
system(sprintf('git -C "%s" worktree remove --force "%s"', project_root, worktree_dir));

%% Paper path on today's code ---------------------------------------------------
assert(startsWith(which("read_rd"), project_root))
paper_today = regression_case_paper(raw_files, kernel_dir, paper_settings);
difference = regression_compare(paper_v1, paper_today);
disp(difference)

%% Production path on today's code --------------------------------------------
work_dir = tempname;
mkdir(work_dir)
prod_today = regression_case_prod(raw_files, times, kernel_dir, work_dir, prod_settings);
rmdir(work_dir, "s")

%% Save -----------------------------------------------------------------------------
[~, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', project_root));
reference = struct("paper", paper_v1, "prod", prod_today, "raw_files", raw_files, "times", times, ...
    "paper_settings", paper_settings, "prod_settings", prod_settings, ...
    "paper_source", "paper-v1", "prod_source", strtrim(commit), "created", datetime("now"));
save(fullfile(reference_dir, "paper_v1_reference.mat"), "reference", "-v7.3")
