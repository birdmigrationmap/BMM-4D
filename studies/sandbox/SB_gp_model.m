%% Sandbox: Gaussian-covariance model without a radar prior (research)
% The field x (birds km^-3 on the active cells, all frames) minimizes
%   ||W (F x - d)||^2                    data: every observation row
% + w_p^2 (x - m)' Q (x - m)             Gaussian prior around a mean m
% + w_a^2 ||A x||^2                      advection between frames
% (+ w_s^2 ||L x||^2, the smoothness term, switched off here: the
%   covariance takes its role).
% - Mean m: one vertical profile per frame for the whole domain (the mean
%   over the domain of the radar profiles), no horizontal structure, so
%   nothing radar-specific is put in where the data are missing.
% - Covariance: Q = c (I + (l_h/dx)^2 L_h + (l_v/dz)^2 L_v), L_h and L_v the
%   horizontal and vertical graph Laplacians. Q^-1 is a Matern-like
%   covariance (SPDE form) correlated over about l_h horizontally and l_v
%   vertically; c makes each cell's variance that of the white prior, so
%   w_p sets how far the field may move from m and l_h, l_v only how far
%   a deviation spreads. This is a Gaussian process (kriging) on the grid,
%   with the beams as data and advection linking the frames.
% - Assumptions: deviations from m are stationary (same variance and
%   correlation everywhere, in birds km^-3, not relative to the local
%   level); the data error is the same for every row; advection and the
%   take-off/landing relaxation as in production.
% Compared with the radar-prior model (reference): the white prior around
% the interpolated radar profiles (lambda_profile = 1) plus smoothness.
% Outputs: figures/gp_model.png (sandbox_case_figure).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var") || ~isstruct(S)
    S = sandbox_load("sandbox");
end

%% Solve ------------------------------------------------------------------
options = {
    struct("name", "Reference: radar prior (smooth 1), white", "lambda_profile", 1)
    struct("name", "GP 150 km x 300 m, flat mean", "lambda_profile", 1, "prior_mean", "profile", ...
        "prior_correlation_m", 150000, "prior_correlation_vertical_m", 300, "smoothness", 0)
    struct("name", "GP 50 km x 300 m, flat mean", "lambda_profile", 1, "prior_mean", "profile", ...
        "prior_correlation_m", 50000, "prior_correlation_vertical_m", 300, "smoothness", 0)};
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    option = options{i_o};
    option.video = false;
    results{i_o} = sandbox_run(S, option);
end
disp(table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), cellfun(@(r) r.misfit_trusted, results), ...
    cellfun(@(r) r.detail, results), cellfun(@(r) r.negative_share, results), cellfun(@(r) r.seconds, results), ...
    VariableNames=["case" "bias" "trusted_misfit" "detail" "negative_share" "seconds"]))

%% Plot results -----------------------------------------------------------
sandbox_case_figure(S, results, fullfile(figure_dir, "gp_model.png"), 4, "Gaussian covariance without a radar prior");
