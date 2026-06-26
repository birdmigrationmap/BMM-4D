function system = prepare_scans_advection(observations, grid, steps, varargin)
%PREPARE_SCANS_ADVECTION Build reusable matrix-free inversion components.

p = inputParser;
addParameter(p, "water", false);
addParameter(p, "coeff_horver", 0);
addParameter(p, "smoothness", true);
parse(p, varargin{:});

F = observations.F;
obs_values = observations.values;
n_time = size(obs_values, 2);
nxyz = prod(grid.sz);
steps = steps(:);
assert(numel(steps) == max(n_time - 1, 0), "prepare_scans_advection:InvalidStepCount", ...
    "Expected %d advection steps for %d time slices.", max(n_time - 1, 0), n_time)

obs_mask = full(obs_values ~= 0);
F_squared = F .^ 2;
row_norm_squared = full(sum(F_squared, 2));
data_scales = 1 ./ sqrt(max(sum(obs_mask, 1), 1));

data_diagonal = zeros(nxyz, n_time);
data_rhs = zeros(nxyz, n_time);
data_block_norms = zeros(1, n_time);
for i_t = 1:n_time
    mask_t = obs_mask(:, i_t);
    data_block_norms(i_t) = max(data_scales(i_t) .* sqrt(sum(row_norm_squared(mask_t))), 1);
    data_diagonal(:, i_t) = data_scales(i_t) ^ 2 .* (F_squared' * double(mask_t));
    data_rhs(:, i_t) = data_scales(i_t) ^ 2 .* (F' * full(obs_values(:, i_t)));
end
clear F_squared

if p.Results.smoothness
    laplacian_3d = LaplacianMatrix(grid, water=p.Results.water, coeff_horver=p.Results.coeff_horver);
    laplacian_norm = norm(laplacian_3d, "fro");
    laplacian_normal_diagonal = full(sum(laplacian_3d .^ 2, 1))';
else
    laplacian_3d = [];
    laplacian_norm = 0;
    laplacian_normal_diagonal = [];
end

advection_base_norms = zeros(max(n_time - 1, 0), 1);
for i_t = 1:n_time - 1
    coeff = steps{i_t}.coeff_valid;
    if isscalar(coeff)
        coefficient_norm_squared = numel(steps{i_t}.source_id_valid) .* coeff .^ 2;
    else
        coefficient_norm_squared = sum(coeff .^ 2);
    end
    advection_base_norms(i_t) = sqrt(2 .* coefficient_norm_squared);
end

system = struct( ...
    "F", F, ...
    "obs_values", obs_values, ...
    "obs_mask", obs_mask, ...
    "data_scales", data_scales, ...
    "data_block_norms", data_block_norms, ...
    "data_diagonal", data_diagonal, ...
    "data_rhs", data_rhs, ...
    "laplacian", laplacian_3d, ...
    "laplacian_norm", laplacian_norm, ...
    "laplacian_normal_diagonal", laplacian_normal_diagonal, ...
    "steps", {steps}, ...
    "advection_base_norms", advection_base_norms, ...
    "grid_sz", grid.sz, ...
    "nxyz", nxyz, ...
    "n_time", n_time, ...
    "water", p.Results.water, ...
    "coeff_horver", p.Results.coeff_horver);
end
