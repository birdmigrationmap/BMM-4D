%% Sandbox: plain PCG against PCG + projected FISTA
% The non-negativity step (projected FISTA after PCG) is the slow part of
% the lean solver. This compares the same options with and without it, to
% see whether plain PCG is good enough for fast option tests: speed, radar
% bias, misfit, detail and the share of negative mass. Each run writes a
% video to prod/sandbox/figures/runs; the table and comparison figure go
% to prod/sandbox/figures/fast_solver.*.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"), fullfile(project_root, "prod", "sandbox"))
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end

%% Options ----------------------------------------------------------------
% Data weights from SB_option_sweep: hard cuts (us_east default) and the
% logistic near-range / no-echo transitions
hard_cut = @(r, e, z) double(r >= 25000 & ~(z & r > 75000));
logistic_28 = @(r, e, z) 1 ./ (1 + exp(-(r - 28000) / 3000)) .* (~z + z ./ (1 + exp((r - 75000) / 10000)));
base_options = {
    struct("name", "no cuts")
    struct("name", "hard cuts", "row_weight", hard_cut)
    struct("name", "logistic 28, prior 0.5", "row_weight", logistic_28, "prior", 0.5)};
solver_label = ["PCG" "FISTA"];
options = {};
for i_o = 1:numel(base_options)
    for nonnegative = [false true]
        option = base_options{i_o};
        option.nonnegative = nonnegative;
        option.name = option.name + ", " + solver_label(nonnegative + 1);
        options{end + 1, 1} = option; %#ok<SAGROW>
    end
end

%% Run --------------------------------------------------------------------
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = sandbox_run(S, options{i_o});
    fprintf("[fast] %-34s bias %.3f  trusted %.3f  detail %.2f  negative %.1f%%  %3.0f s (%d PCG it.)\n", ...
        results{i_o}.name, results{i_o}.bias_200km, results{i_o}.misfit_trusted, results{i_o}.detail, ...
        100 * results{i_o}.negative_share, results{i_o}.seconds, results{i_o}.iterations)
end
summary = sandbox_summary(results, "fast_solver");

