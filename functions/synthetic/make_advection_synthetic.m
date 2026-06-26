function [truth_sequences, obs_sequences] = make_advection_synthetic(grid, F, cov_parm, vp_case, trans_finv_case, target_id, n_sim, uv_profile, time_offsets, dt, varargin)
%MAKE_ADVECTION_SYNTHETIC Generate moving truth sequences and observations.

p = inputParser;
addParameter(p, "seed", "shuffle");
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
parse(p, varargin{:});

[generation_grid, crop] = build_advection_generation_grid( ...
    grid, uv_profile, time_offsets, seconds(dt));

truth_fields = simulate_gaussian_truth( ...
    generation_grid, cov_parm, trans_finv_case{target_id}, vp_case(:, target_id), ...
    n_sim, seed=p.Results.seed);

steps = build_profile_advection_steps( ...
    generation_grid, uv_profile, seconds(dt), ...
    boundary_mode=p.Results.boundary_mode, ...
    shoreline_crossing_coeff=p.Results.shoreline_crossing_coeff);

truth_sequences = cell(n_sim, 1);
obs_sequences = cell(n_sim, 1);
for i_sim = 1:n_sim
    truth_sequence_padded = advect_density_sequence(truth_fields{i_sim}, steps, time_offsets);
    truth_sequence = truth_sequence_padded(crop.x_idx, crop.y_idx, :, :);
    truth_sequences{i_sim} = truth_sequence;

    obs_rows = cell(numel(time_offsets), 1);
    obs_cols = cell(numel(time_offsets), 1);
    obs_vals = cell(numel(time_offsets), 1);
    for i_t = 1:numel(time_offsets)
        field_t = truth_sequence(:, :, :, i_t);
        obs_t = F * field_t(:);
        obs_rows{i_t} = find(obs_t);
        obs_cols{i_t} = repmat(i_t, numel(obs_rows{i_t}), 1);
        obs_vals{i_t} = obs_t(obs_rows{i_t});
    end
    obs_sequence = sparse(vertcat(obs_rows{:}), vertcat(obs_cols{:}), vertcat(obs_vals{:}), size(F, 1), numel(time_offsets));
    obs_sequences{i_sim} = obs_sequence;
end
end
