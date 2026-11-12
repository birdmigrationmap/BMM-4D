%% Sandbox: Gaussian-Markov prior on the deviation from the profile prior
% With the profile prior fitted at lambda_profile = 1 (SB_profile_prior),
% the deviation from the prior gets a first-order Matern / SPDE
% Gaussian-Markov prior (solve_advection_lean prior_correlation_m and
% prior_correlation_vertical_m): correlated horizontally over l_h and
% vertically over l_v, instead of independent per cell. Retested on the
% estimator data (no cuts); on 23 October it was tested on the median
% cleaning, where it flattened the imprint but cost 25-35% trusted misfit.
% Outputs: figures/gp_prior.png, gp_prior.csv (sandbox_summary).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
if ~exist("S", "var") || ~isstruct(S)
    S = sandbox_load("sandbox");
end

%% Solve ------------------------------------------------------------------
base = struct("lambda_profile", 1, "video", false);
options = {
    struct("name", "white prior 0.3")
    struct("name", "GP 30 km", "prior_correlation_m", 30000)
    struct("name", "GP 50 km", "prior_correlation_m", 50000)
    struct("name", "GP 100 km", "prior_correlation_m", 100000)
    struct("name", "GP 50 km, vertical 300 m", "prior_correlation_m", 50000, "prior_correlation_vertical_m", 300)
    struct("name", "GP 50 km, prior 0.1", "prior_correlation_m", 50000, "prior", 0.1)};
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    option = options{i_o};
    for name = string(fieldnames(base))'
        option.(name) = base.(name);
    end
    results{i_o} = sandbox_run(S, option);
end

%% Plot results -----------------------------------------------------------
sandbox_summary(results, "gp_prior");

%% Trade-off: white against correlated prior, over prior weights ---------------
tradeoff = {};
for weight = [0.3 0.5 1 2]
    tradeoff{end + 1, 1} = sandbox_run(S, struct("name", sprintf("white %g", weight), "prior", weight, ...
        "lambda_profile", 1, "video", false)); %#ok<SAGROW>
end
for weight = [0.5 1 2 4]
    tradeoff{end + 1, 1} = sandbox_run(S, struct("name", sprintf("GP 50 km, %g", weight), "prior", weight, ...
        "prior_correlation_m", 50000, "lambda_profile", 1, "video", false)); %#ok<SAGROW>
end
sandbox_summary(tradeoff, "gp_tradeoff_normalized");
