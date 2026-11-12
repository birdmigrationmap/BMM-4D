function ar_terrain_solves(c, variants, held_out, section_time, results_file, n_threads)
%AR_TERRAIN_SOLVES Solves and held-out radars for the vertical rules over terrain (dev).
%   AR_TERRAIN_SOLVES(C, VARIANTS, HELD_OUT, SECTION_TIME, RESULTS_FILE,
%   N_THREADS) loads case C (ar_load), and for each variant (struct array:
%   name, steps, steps_options, prior, prior_options; modes and name-value
%   options of ar_steps and ar_profile_prior) solves the case (ar_solve,
%   white prior), then leaves out each radar of HELD_OUT in turn (rows
%   weight 0, prior rebuilt from the other radars' profiles) and scores
%   the prediction of its rows at 30-120 km (ar_holdout_score), as
%   AR_G_terrain_advection. Meant for a batch job (N_THREADS threads).
%   RESULTS_FILE gets the tables scores and holdout, the column maps
%   (single [nx ny n_time] per variant) and the fields at SECTION_TIME
%   (single [n_active 1] per variant).

maxNumCompThreads(n_threads);
S = ar_load(c);
inv = c.inversion;
dt_seconds = seconds(diff(c.times));
i_section = find(c.times == section_time);

scores = table();
columns = cell(numel(variants), 1);
section_fields = cell(numel(variants), 1);
steps = cell(numel(variants), 1);
for i_v = 1:numel(variants)
    v = variants(i_v);
    steps{i_v} = ar_steps(S.g, S.radar_scans, dt_seconds, S.active_mask_3d, inv, v.steps, v.steps_options{:});
    prior = ar_profile_prior(S.g, S.radar_scans, S.radar_profiles, S.active_mask_3d, inv.prior_length_scale_m, ...
        v.prior, v.prior_options{:});
    result = ar_solve(S, v.name, steps=steps{i_v}, prior=prior);
    score = ar_terrain_score(S, result.column);
    scores = [scores; table(string(v.name), score.ridge_valley, score.per_km, result.misfit, result.misfit_trusted, ...
        result.seconds, result.iterations, VariableNames=["solve" "ridge_valley" "per_km" "misfit" ...
        "misfit_30_120km" "seconds" "iterations"])]; %#ok<AGROW>
    columns{i_v} = single(result.column);
    section_fields{i_v} = single(result.field(:, i_section));
    save(results_file, "scores", "columns", "section_fields", "variants")
end
disp(scores)

holdout = table();
for radar_name = held_out(:)'
    for i_v = 1:numel(variants)
        v = variants(i_v);
        score = ar_holdout_score(S, radar_name, @(S_out, keep) ar_solve(S_out, v.name + " without " + radar_name, ...
            steps=steps{i_v}, prior=ar_profile_prior(S.g, S.radar_scans(keep, :), S.radar_profiles(:, keep, :), ...
            S.active_mask_3d, inv.prior_length_scale_m, v.prior, v.prior_options{:})));
        holdout = [holdout; table(radar_name, string(v.name), score.relative_error, score.median_abs_log_ratio, ...
            VariableNames=["radar" "solve" "relative_error" "median_abs_log_ratio"])]; %#ok<AGROW>
        save(results_file, "scores", "columns", "section_fields", "variants", "holdout")
    end
end
disp(holdout)
end
