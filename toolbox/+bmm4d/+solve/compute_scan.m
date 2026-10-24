function [scan, inversion_info] = compute_scan(F, obs_vector, g, varargin)
%COMPUTE_SCAN Reconstruct one 3D density field from one radar observation vector.
%   [SCAN, INVERSION_INFO] = COMPUTE_SCAN(F, OBS_VECTOR, GRID, ...) solves
%   a regularized least-squares inverse problem for a single time step. F
%   is the sparse forward operator, OBS_VECTOR is the sparse observation
%   vector for one time, and GRID is the reconstruction g from
%   CREATE_GRID.
%
%   The unknown is the full 3D density field on the Cartesian g. The
%   optional regularization terms add prior information on:
%   - the overall magnitude of the field (identity)
%   - smoothness in x, y, and z (smoothness)
%   - agreement with a 1D vertical profile in absolute value (vp_absolute)
%   Each regularization block is scaled against the fitted data block so
%   these weights stay comparable across observation counts and operators.

p = inputParser;
addParameter(p,'identity',0);
addParameter(p,'smoothness',0);
addParameter(p,'vp_absolute',0);
addParameter(p,'water',true);
addParameter(p,'tol', 1e-6); % same default as lsqr
addParameter(p,'maxit', 100); % same default as lsqr
addParameter(p, 'coeff_horver', 0); % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. between 0 and .5
parse(p,varargin{:});

num_grid_cells = prod(g.sz);
num_horizontal_cells = prod(g.sz(1:2));

% Keep only non-zero observations. Zero entries do not contribute to the
% least-squares system and are dropped together with the corresponding rows
% of the forward operator.
[obs_values, nonzero_obs_idx] = bmm4d.forward.format_data(obs_vector(:,1));
data_scale = 1 / sqrt(max(numel(obs_values), 1));

% Scale regularization blocks against the fitted data block so lambda
% values stay comparable across observation counts and operators.
A_system = data_scale .* F(nonzero_obs_idx,:);
b_system = data_scale .* obs_values;
data_block_norm = max(norm(A_system, "fro"), 1);

% Penalize large voxel values directly.
if p.Results.identity > 0
    identity_operator = speye(num_grid_cells);
    identity_scale = p.Results.identity .* data_block_norm ./ sqrt(num_grid_cells);
    A_system = [A_system; identity_scale .* identity_operator];
    b_system = [b_system; zeros(num_grid_cells,1)];
end

% Penalize roughness in the 3D field with the project Laplacian operator.
if p.Results.smoothness > 0
    laplacian_3d = bmm4d.grid.LaplacianMatrix(g, water=p.Results.water, coeff_horver=p.Results.coeff_horver);
    smoothness_scale = p.Results.smoothness .* data_block_norm ./ max(norm(laplacian_3d, "fro"), eps);
    A_system = [A_system; smoothness_scale .* laplacian_3d];
    b_system = [b_system; sparse(num_grid_cells,1)];
end

% Penalize deviation from the vertical profile estimated from the same
% observations in every voxel: d(x,y,z) ~= vp(z). The profile is repeated
% across all horizontal cells (operator U in the notes) and set to zero
% over water when the water mask is active.
if p.Results.vp_absolute > 0
    vp = bmm4d.prior.compute_vp(obs_vector(:,1), g, F);
    vp_field = repelem(vp(:), num_horizontal_cells, 1);
    if p.Results.water
        water_mask_3d = repmat(bmm4d.grid.get_water_mask(g), 1, 1, g.sz(3));
        vp_field(water_mask_3d(:) == 1) = 0;
    end

    vp_absolute_scale = p.Results.vp_absolute .* data_block_norm ./ sqrt(num_grid_cells);
    A_system = [A_system; vp_absolute_scale .* speye(num_grid_cells)];
    b_system = [b_system; vp_absolute_scale .* vp_field];
end

% Solve the regularized inverse problem and reshape it back to the canonical
% Cartesian g order [x y z].
[scan_1D, flag, relres, iter] = lsqr(A_system, b_system, p.Results.tol, p.Results.maxit);

% Keep orientation changes out of shared helpers. Plotting code can permute
% arrays locally if a display-oriented layout is needed.
scan = reshape(scan_1D,g.sz);

inversion_info = struct( ...
    "controls", orderfields(struct(p.Results)), ...
    "solver", struct( ...
        "flag", flag, ...
        "relres", relres, ...
        "iter", iter, ...
        "n_data", numel(obs_values), ...
        "n_unknown", num_grid_cells));
end
