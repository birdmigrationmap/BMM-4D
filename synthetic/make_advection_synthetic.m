function truth_sequences = make_advection_synthetic(g, cov_parm, vp_case, trans_finv_case, target_id, n_sim, uv_profile, time_offsets, dt, varargin)
%MAKE_ADVECTION_SYNTHETIC Generate moving truth sequences.
% Observations are rebuilt with forward_observation_sequence(F, truth).

p = inputParser;
addParameter(p, "seed", "shuffle");
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
parse(p, varargin{:});

[generation_grid, crop] = build_advection_generation_grid( ...
    g, uv_profile, time_offsets, seconds(dt));

truth_fields = simulate_gaussian_truth( ...
    generation_grid, cov_parm, trans_finv_case{target_id}, vp_case(:, target_id), ...
    n_sim, seed=p.Results.seed);

steps = bmm4d.advect.build_profile_advection_steps( ...
    generation_grid, uv_profile, seconds(dt), ...
    boundary_mode=p.Results.boundary_mode, ...
    shoreline_crossing_coeff=p.Results.shoreline_crossing_coeff);

truth_sequences = cell(n_sim, 1);
for i_sim = 1:n_sim
    truth_sequence_padded = advect_density_sequence(truth_fields{i_sim}, steps, time_offsets);
    truth_sequences{i_sim} = truth_sequence_padded(crop.x_idx, crop.y_idx, :, :);
end
end
