function fields = advect_density_sequence(reference_field, steps, time_offsets)
%ADVECT_DENSITY_SEQUENCE Shift one reference field into a synthetic time sequence.
% steps contains one interval step between consecutive time offsets.

grid_size = size(reference_field);
if numel(grid_size) < 3
    grid_size(3) = 1;
end
n_time = numel(time_offsets);
fields = zeros([grid_size n_time], "like", reference_field);
center_id = find(time_offsets == 0, 1);

if isempty(center_id)
    error("advect_density_sequence:MissingReferenceTime", ...
        "time_offsets must include 0 for the reference field.");
end

if ~iscell(steps) || numel(steps) ~= n_time - 1
    error("advect_density_sequence:InvalidSteps", ...
        "steps must be a cell array with one interval step per consecutive time pair.");
end

fields(:, :, :, center_id) = reference_field;

for i_t = center_id + 1:n_time
    fields(:, :, :, i_t) = shift_field(fields(:, :, :, i_t - 1), steps{i_t - 1}, 1);
end

for i_t = center_id - 1:-1:1
    fields(:, :, :, i_t) = shift_field(fields(:, :, :, i_t + 1), steps{i_t}, -1);
end
end

function field_out = shift_field(field_in, step, direction)

% Apply the same compact mappings used by the advection inversion.
field_out = zeros(size(field_in), "like", field_in);
if direction == 1
    field_out(step.dest_id_valid) = step.coeff_valid .* field_in(step.source_id_valid);
else
    field_out(step.source_id_valid) = step.coeff_valid .* field_in(step.dest_id_valid);
end
end
