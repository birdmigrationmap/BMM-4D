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
%   Every block is a root-mean-square: the data residual is averaged over
%   the observations and each regularization residual over the grid cells
%   (weight / sqrt(n_cells)), with the smoothness operator in physical
%   units. The weights are therefore dimensionless and do not depend on the
%   density magnitude, the number of observations or the grid resolution.
%
%   The least-squares problem is solved through its normal equations with
%   Jacobi-preconditioned conjugate gradients (PCG), started from the
%   vertical-profile field when that prior is active. 'tol' is the PCG
%   relative residual tolerance and 'maxit' the PCG iteration limit.

p = inputParser;
addParameter(p,'identity',0);
addParameter(p,'smoothness',0);
addParameter(p,'vp_absolute',0);
addParameter(p,'water',true);
addParameter(p,'tol', 1e-6);
addParameter(p,'maxit', 1000);
addParameter(p, 'coeff_horver', 0); % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. between 0 and .5
parse(p,varargin{:});

num_grid_cells = prod(g.sz);
num_horizontal_cells = prod(g.sz(1:2));

% Keep only non-zero observations. Zero entries do not contribute to the
% least-squares system and are dropped together with the corresponding rows
% of the forward operator.
[obs_values, nonzero_obs_idx] = bmm4d.forward.format_data(obs_vector(:,1));
data_scale = 1 / sqrt(max(numel(obs_values), 1));

% Data block of the normal equations (mean squared residual).
A_data = data_scale .* F(nonzero_obs_idx,:);
% Regularization blocks are mean squares over the grid cells.
cell_scale = 1 / sqrt(num_grid_cells);
normal_matrix = A_data' * A_data;
normal_rhs = A_data' * (data_scale .* obs_values);
clear A_data

% Penalize large voxel values directly.
if p.Results.identity > 0
    identity_scale = p.Results.identity .* cell_scale;
    normal_matrix = normal_matrix + identity_scale^2 .* speye(num_grid_cells);
end

% Penalize roughness in the 3D field with the physical-unit Laplacian.
if p.Results.smoothness > 0
    laplacian_3d = bmm4d.grid.LaplacianMatrix(g, water=p.Results.water, coeff_horver=p.Results.coeff_horver);
    smoothness_scale = p.Results.smoothness .* cell_scale;
    normal_matrix = normal_matrix + smoothness_scale^2 .* (laplacian_3d' * laplacian_3d);
end

% Penalize deviation from the vertical profile estimated from the same
% observations in every voxel: d(x,y,z) ~= vp(z). The profile is repeated
% across all horizontal cells (operator U in the notes) and set to zero
% over water when the water mask is active.
scan_initial = zeros(num_grid_cells, 1);
if p.Results.vp_absolute > 0
    vp = bmm4d.prior.compute_vp(obs_vector(:,1), g, F);
    vp_field = repelem(vp(:), num_horizontal_cells, 1);
    if p.Results.water
        water_mask_3d = repmat(bmm4d.grid.get_water_mask(g), 1, 1, g.sz(3));
        vp_field(water_mask_3d(:) == 1) = 0;
    end

    vp_absolute_scale = p.Results.vp_absolute .* cell_scale;
    normal_matrix = normal_matrix + vp_absolute_scale^2 .* speye(num_grid_cells);
    normal_rhs = normal_rhs + vp_absolute_scale^2 .* vp_field;
    % Start from the profile field: with small prior weights, voxels far
    % from the radar move away from a zero start only very slowly.
    scan_initial = vp_field;
end

% Jacobi preconditioner; voxels without any constraint keep a unit scale.
normal_diagonal = full(diag(normal_matrix));
normal_diagonal(normal_diagonal <= 0) = 1;
preconditioner = spdiags(normal_diagonal, 0, num_grid_cells, num_grid_cells);

% Solve the regularized inverse problem and reshape it back to the canonical
% Cartesian g order [x y z].
[scan_1D, flag, relres, iter] = pcg(normal_matrix, full(normal_rhs), p.Results.tol, p.Results.maxit, preconditioner, [], scan_initial);

% Keep orientation changes out of shared helpers. Plotting code can permute
% arrays locally if a display-oriented layout is needed.
scan = reshape(scan_1D,g.sz);

inversion_info = struct( ...
    "controls", orderfields(struct(p.Results)), ...
    "solver", struct( ...
        "name", "pcg_normal_jacobi", ...
        "flag", flag, ...
        "relres", relres, ...
        "iter", iter, ...
        "n_data", numel(obs_values), ...
        "n_unknown", num_grid_cells));
end
