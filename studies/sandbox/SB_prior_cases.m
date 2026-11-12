%% Sandbox: prior and reconstruction, side by side, for the prior options
% Two routes for what constrains the field where the data do not:
% 1. A radar-profile prior (mean field) with a white prior around it: the
%    current one (lambda_profile = 10) and the less smoothed profiles
%    (lambda_profile = 1, SB_profile_prior).
% 2. No radar prior: one domain-mean vertical profile per time (no
%    horizontal structure) as the mean, and a Gaussian-Markov covariance on
%    the reconstructed grid (prior_correlation_m, normalized to the white
%    prior's variance) that interpolates between the data: 50 and 150 km.
% (smooth = lambda_profile.) For each case, one frame (02:00 UTC): column density of the prior mean
% and of the reconstruction, the west-east section through KDMX of both, and
% the scores (sandbox_case_figure).
% Outputs: figures/prior_cases.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var") || ~isstruct(S)
    S = sandbox_load("sandbox");
end
cases = {
    struct("name", "1a. Radar prior, smooth 10", "lambda_profile", 10)
    struct("name", "1b. Radar prior, smooth 1", "lambda_profile", 1)
    struct("name", "2a. Flat mean, GP 50 km", "lambda_profile", 1, "prior_mean", "profile", "prior_correlation_m", 50000)
    struct("name", "2b. Flat mean, GP 150 km", "lambda_profile", 1, "prior_mean", "profile", "prior_correlation_m", 150000)};
i_t = 4;

%% Solve ------------------------------------------------------------------
n_case = numel(cases);
results = cell(n_case, 1);
for i_c = 1:n_case
    option = cases{i_c};
    option.video = false;
    results{i_c} = sandbox_run(S, option);
end
disp(table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), cellfun(@(r) r.misfit_trusted, results), ...
    cellfun(@(r) r.detail, results), cellfun(@(r) r.seconds, results), VariableNames=["case" "bias" "trusted_misfit" "detail" "seconds"]))

%% Plot results -----------------------------------------------------------
sandbox_case_figure(S, results, fullfile(figure_dir, "prior_cases.png"), i_t, "Prior options");
