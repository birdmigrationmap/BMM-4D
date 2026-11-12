function scores = ar_matern_check(project_root, out_file)
%AR_MATERN_CHECK Cost and edges of the one-term Matern prior on ar_plains_bb (dev).
%   SCORES = AR_MATERN_CHECK(PROJECT_ROOT, OUT_FILE) solves ar_plains_bb
%   (12 radars, 04:30-07:30 UTC, beam-blockage fix) with:
%   - V0 white prior + smoothness 0.3 (production);
%   - V1 correlated prior 100 km + smoothness 0.3 (fix 3);
%   - B one Matern term, 25 km x 190 m and 50 km x 190 m, no smoothness
%     (PRIOR_OPTIONS.md, option B);
%   all with the production weights, plain PCG (Jacobi), maxit 5000.
%   Saves the scores (iterations, seconds, edge, misfit) and the column
%   maps to OUT_FILE. Made to run as a MATLAB batch job.

maxNumCompThreads(6);
setup_paths(project_root);
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_plains_bb", [-104 -94 37.5 47.5], times);
S = ar_load(c);
S.c.inversion.maxit = 5000;
variants = struct("name", {"V0 white + smoothness", "V1 GP 100 km + smoothness", "B 25 km x 190 m", "B 50 km x 190 m"}, ...
    "correlation_m", {0, 100000, 0, 0}, "matern_m", {[], [], [25000 190], [50000 190]}, "smoothness", {0.3, 0.3, 0, 0});
scores = table();
columns = cell(numel(variants), 1);
for i_v = 1:numel(variants)
    v = variants(i_v);
    result = ar_solve(S, v.name, prior_correlation_m=v.correlation_m, prior_matern_m=v.matern_m, smoothness=v.smoothness);
    edge = ar_edge_score(S, result.column);
    columns{i_v} = result.column;
    scores = [scores; table(string(v.name), result.iterations, result.seconds, edge.edge, result.misfit, result.misfit_trusted, ...
        VariableNames=["solve" "iterations" "seconds" "edge" "misfit" "misfit_30_120km"])]; %#ok<AGROW>
    save(out_file, "scores", "columns")
end
end
