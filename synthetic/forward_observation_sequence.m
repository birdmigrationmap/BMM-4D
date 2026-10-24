function obs_sequence = forward_observation_sequence(F, truth_sequence)
%FORWARD_OBSERVATION_SEQUENCE Noise-free sparse observations F*x for each time slice.
% Synthetic observations are not saved in case files; they are rebuilt from
% the saved single-precision truth so generation and reuse stay identical.

n_time = size(truth_sequence, 4);
obs_columns = cell(1, n_time);
for i_t = 1:n_time
    field_t = double(truth_sequence(:, :, :, i_t));
    obs_columns{i_t} = sparse(F * field_t(:));
end
obs_sequence = [obs_columns{:}];
end
