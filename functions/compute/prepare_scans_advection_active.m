function system = prepare_scans_advection_active(observations, g, active_mask_3d, steps, varargin)
%PREPARE_SCANS_ADVECTION_ACTIVE Build land/water-subset inversion components.
%   Large per-cell arrays are stored as single precision. F is sliced to
%   active cells, and advection links are remapped from full-grid linear
%   indices to active-domain indices.

p = inputParser;
addParameter(p, "store_class", "single");
parse(p, varargin{:});
store_class = string(p.Results.store_class);

active_id = find(active_mask_3d(:));
n_active = numel(active_id);
n_full = prod(g.sz);
n_time = size(observations.values, 2);
steps = steps(:);
assert(numel(steps) == max(n_time - 1, 0), "prepare_scans_advection_active:InvalidStepCount", ...
    "Expected %d advection steps for %d time slices.", max(n_time - 1, 0), n_time)

F = observations.F(:, active_id);
obs_values = observations.values;
obs_mask = full(obs_values ~= 0);
F_squared = F .^ 2;
row_norm_squared = full(sum(F_squared, 2));
data_scales = single(1 ./ sqrt(max(sum(obs_mask, 1), 1)));

data_diagonal = zeros(n_active, n_time, store_class);
data_rhs = zeros(n_active, n_time, store_class);
data_block_norms = zeros(1, n_time, store_class);
for i_t = 1:n_time
    mask_t = obs_mask(:, i_t);
    data_block_norms(i_t) = single(max(double(data_scales(i_t)) .* sqrt(sum(row_norm_squared(mask_t))), 1));
    data_diagonal(:, i_t) = cast(double(data_scales(i_t)) ^ 2 .* (F_squared' * double(mask_t)), store_class);
    data_rhs(:, i_t) = cast(double(data_scales(i_t)) ^ 2 .* (F' * full(obs_values(:, i_t))), store_class);
end
clear F_squared

steps_active = map_steps_to_active_domain(steps, active_id, n_full);
advection_base_norms = zeros(max(n_time - 1, 0), 1, store_class);
for i_t = 1:n_time - 1
    coeff = steps_active{i_t}.coeff_valid;
    if isscalar(coeff)
        coefficient_norm_squared = numel(steps_active{i_t}.source_id_valid) .* coeff .^ 2;
    else
        coefficient_norm_squared = sum(double(coeff) .^ 2);
    end
    advection_base_norms(i_t) = single(sqrt(2 .* coefficient_norm_squared));
end

system = struct( ...
    "F", F, ...
    "obs_values", obs_values, ...
    "obs_mask", obs_mask, ...
    "data_scales", data_scales, ...
    "data_block_norms", data_block_norms, ...
    "data_diagonal", data_diagonal, ...
    "data_rhs", data_rhs, ...
    "steps", {steps_active}, ...
    "advection_base_norms", advection_base_norms, ...
    "grid_sz", g.sz, ...
    "active_id", uint32(active_id), ...
    "n_full", n_full, ...
    "n_active", n_active, ...
    "n_time", n_time, ...
    "store_class", store_class);
end

function steps_active = map_steps_to_active_domain(steps, active_id, n_full)

active_lookup = zeros(n_full, 1, "uint32");
active_lookup(active_id) = uint32(1:numel(active_id));
steps_active = steps;
for i_step = 1:numel(steps)
    source_active = active_lookup(double(steps{i_step}.source_id_valid));
    dest_active = active_lookup(double(steps{i_step}.dest_id_valid));
    keep = source_active > 0 & dest_active > 0;
    steps_active{i_step}.source_id_valid = source_active(keep);
    steps_active{i_step}.dest_id_valid = dest_active(keep);
    if ~isscalar(steps{i_step}.coeff_valid)
        steps_active{i_step}.coeff_valid = single(steps{i_step}.coeff_valid(keep));
    else
        steps_active{i_step}.coeff_valid = single(steps{i_step}.coeff_valid);
    end
end
end
