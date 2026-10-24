%% Sandbox: radar-imprint option sweep
% Solves the sandbox case (Midwest part of the 21-22 Oct 2026 night,
% 6 frames, 4 km x 100 m, no data cuts in P3) with a list of options and
% scores each one on the radar-distance bias and the data misfit. Each run
% writes a video to studies/sandbox/figures/runs; the table and comparison
% figure (sandbox_summary) go to studies/sandbox/figures. Add options to the list and rerun:
% nothing is decided here.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end

%% Options ----------------------------------------------------------------
% Data weights: hard cuts (current us_east default) and soft tapers.
hard_cut = @(r, e, z) double(r >= 25000 & ~(z & r > 75000));
near_taper = @(r) 1 - exp(-(r / 15000) .^ 2);          % 0 at the radar, 0.67 at 15 km, 0.94 at 25 km
far_zero_taper = @(r, z) ~z + z .* exp(-(r / 75000) .^ 2);   % no-echo rows fade out beyond ~75 km
soft = @(r, e, z) near_taper(r) .* far_zero_taper(r, z);

sweep_name = "option_sweep_2";
near_taper_25 = @(r) 1 - exp(-(r / 25000) .^ 2);       % 0.63 at 25 km, 0.98 at 50 km
soft_25 = @(r, e, z) near_taper_25(r) .* far_zero_taper(r, z);
% Logistic transitions: steep but smooth (near range around 20 km, no-echo
% rows fading around 75 km)
logistic_near = @(r) 1 ./ (1 + exp(-(r - 20000) / 3000));
logistic_zero = @(r, z) ~z + z ./ (1 + exp((r - 75000) / 10000));
logistic = @(r, e, z) logistic_near(r) .* logistic_zero(r, z);
options = { ...
    struct("name", "baseline"), ...
    struct("name", "hard cuts 25-75 km", "row_weight", hard_cut), ...
    struct("name", "soft tapers 15 km", "row_weight", soft), ...
    struct("name", "soft tapers 25 km", "row_weight", soft_25), ...
    struct("name", "soft 25 + shape 1", "row_weight", soft_25, "shape", 1), ...
    struct("name", "GP 100 km", "prior_correlation_m", 100000), ...
    struct("name", "GP 100 km, prior 0.1", "prior_correlation_m", 100000, "prior", 0.1), ...
    struct("name", "soft 25 + GP 100 km, prior 0.1", "row_weight", soft_25, "prior_correlation_m", 100000, "prior", 0.1), ...
    struct("name", "soft 25 + GP 100 km + shape 1", "row_weight", soft_25, "prior_correlation_m", 100000, "shape", 1), ...
    struct("name", "soft + GP 150 km + shape 1", "row_weight", soft, "prior_correlation_m", 150000, "shape", 1), ...
    struct("name", "logistic 20-75 km", "row_weight", logistic), ...
    struct("name", "logistic + shape 1", "row_weight", logistic, "shape", 1), ...
    struct("name", "logistic + GP 50 km", "row_weight", logistic, "prior_correlation_m", 50000), ...
    struct("name", "logistic + GP 50 km + shape 1", "row_weight", logistic, "prior_correlation_m", 50000, "shape", 1)};

% An only_options list in the workspace runs only those options
if exist("only_options", "var")
    options = options(cellfun(@(o) ismember(string(o.name), only_options), options));
end

%% Run --------------------------------------------------------------------
results = cell(numel(options), 1);
for i_o = 1:numel(options)
    results{i_o} = sandbox_run(S, options{i_o});
    fprintf("[sandbox] %-32s bias %.3f  misfit %.3f  trusted %.3f  detail %.2f  %4.0f s\n", results{i_o}.name, ...
        results{i_o}.bias_200km, results{i_o}.misfit, results{i_o}.misfit_trusted, results{i_o}.detail, results{i_o}.seconds)
end
summary = sandbox_summary(results, sweep_name);
