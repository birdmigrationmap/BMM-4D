function [fit, timing] = ar_fit_ranges_check(project_root, n_pairs_per_frame)
%AR_FIT_RANGES_CHECK Range fit on the cached ar_plains_bb case (dev, batch job).
%   [FIT, TIMING] = AR_FIT_RANGES_CHECK(PROJECT_ROOT, N_PAIRS_PER_FRAME)
%   returns the fit and the 10 slowest functions of the profile.
maxNumCompThreads(6);
setup_paths(project_root);
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_plains_bb", [-104 -94 37.5 47.5], times);
S = ar_load(c);
profile on
fit_tic = tic;
fit = ar_fit_ranges(S.observations, S.g, S.prior, S.active_id, n_pairs_per_frame);
fit.seconds = toc(fit_tic);
profile off
info = profile("info");
table_all = struct2table(info.FunctionTable);
timing = sortrows(table_all(:, ["FunctionName" "TotalTime"]), "TotalTime", "descend");
timing = timing(1:min(10, height(timing)), :);
end
