function [scan, inversion_info] = compute_scan(F, obs_vector, grid, varargin)
%COMPUTE_SCAN Reconstruct one 3D density field from one radar observation vector.
%   [SCAN, INVERSION_INFO] = COMPUTE_SCAN(F, OBS_VECTOR, GRID, ...) solves
%   a regularized least-squares inverse problem for a single time step. F
%   is the sparse forward operator, OBS_VECTOR is the sparse observation
%   vector for one time, and GRID is the reconstruction grid from
%   CREATE_GRID.
%
%   The unknown is the full 3D density field on the Cartesian grid. The
%   optional regularization terms add prior information on:
%   - the overall magnitude of the field (identity)
%   - smoothness in x, y, and z (smoothness)
%   - agreement with a 1D vertical profile in absolute value (vp_absolute)
%   - agreement with a 1D vertical profile in relative shape (vp_ratio)
%   Each regularization block is scaled against the fitted data block so
%   these weights stay comparable across observation counts and operators.

p = inputParser;
addParameter(p,'identity',0);
addParameter(p,'smoothness',0);
addParameter(p,'vp_absolute',0);
addParameter(p,'vp_ratio',0);
addParameter(p,'water',true);
addParameter(p,'tol', 1e-6); % same default as lsqr
addParameter(p,'maxit', 100); % same default as lsqr
addParameter(p, 'coeff_horver', 0); % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. between 0 and .5
parse(p,varargin{:});

num_grid_cells = prod(grid.sz);
num_horizontal_cells = prod(grid.sz(1:2));
use_vertical_profile = p.Results.vp_absolute > 0 || p.Results.vp_ratio > 0;

if p.Results.water
    if isfield(grid,'mask_water')
        water_mask = grid.mask_water;
    else
        water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
        water_mask = logical(water_mask_data.Fonwater({grid.f_lat(grid.y),grid.f_lon(grid.x)}))';
    end
    water_mask_3d = repmat(water_mask, 1, 1, grid.sz(3));
end

% Keep only non-zero observations. Zero entries do not contribute to the
% least-squares system and are dropped together with the corresponding rows
% of the forward operator.
[obs_values, nonzero_obs_idx] = format_data(obs_vector(:,1));
data_scale = 1 / sqrt(max(numel(obs_values), 1));

% Scale regularization blocks against the fitted data block so lambda
% values stay comparable across observation counts and operators.
A_system = data_scale .* F(nonzero_obs_idx,:);
b_system = data_scale .* obs_values;
data_block_norm = max(norm(A_system, "fro"), 1);

% Penalize large voxel values directly.
if p.Results.identity > 0
    identity_operator = speye(num_grid_cells);
    identity_scale = p.Results.identity .* data_block_norm ./ max(norm(identity_operator, "fro"), eps);
    A_system = [A_system; identity_scale .* identity_operator];
    b_system = [b_system; zeros(num_grid_cells,1)];
end

% Penalize roughness in the 3D field with the project Laplacian operator.
if p.Results.smoothness > 0
    laplacian_3d = LaplacianMatrix(grid, water=p.Results.water, coeff_horver=p.Results.coeff_horver);
    smoothness_scale = p.Results.smoothness .* data_block_norm ./ max(norm(laplacian_3d, "fro"), eps);
    A_system = [A_system; smoothness_scale .* laplacian_3d];
    b_system = [b_system; sparse(num_grid_cells,1)];
end

% Estimate the vertical profile from the same observations whenever a
% VP-based regularization term is active.
if use_vertical_profile
    vp = compute_vp(obs_vector(:,1), grid, F);
end

% Expand a 1D vertical profile to the full 3D grid by repeating each
% altitude value across all horizontal cells. This is the operator U in the
% notes used to map vp(z) to a full field d(x,y,z) = vp(z).
if p.Results.vp_absolute > 0
    vertical_profile_expansion = repelem(speye(grid.sz(3)), num_horizontal_cells, 1);

    if p.Results.water
        vertical_profile_expansion(water_mask_3d(:) == 1, :) = 0;
    end

    % Penalize deviation from the absolute profile value in every voxel:
    % d(x,y,z) ~= vp(z).
    vp_absolute_scale = p.Results.vp_absolute .* data_block_norm ./ max(norm(speye(num_grid_cells), "fro"), eps);
    A_system = [A_system; vp_absolute_scale .* speye(num_grid_cells)];
    b_system = [b_system; vp_absolute_scale .* vertical_profile_expansion * vp];
end

% Penalize deviation from the normalized profile shape while allowing each
% horizontal column to keep its own total magnitude:
% d(:,x,y) ~= s(x,y) * vp_ratio.
if p.Results.vp_ratio > 0
    vp_nonnegative = vp;
    vp_nonnegative(vp_nonnegative < 0) = 0;
    normalized_vp = vp_nonnegative ./ sum(vp_nonnegative);

    % For each voxel, Q*d = 0 enforces equality between the voxel value and
    % the profile-shaped share of the total mass in its vertical column.
    repeated_column_identity = repmat(speye(num_horizontal_cells), grid.sz(3), grid.sz(3));
    profile_ratio_operator = repeated_column_identity .* repelem(normalized_vp', 1, num_horizontal_cells) - speye(num_grid_cells);

    if p.Results.water
        profile_ratio_operator(water_mask_3d(:) == 1, :) = 0;
    end

    vp_ratio_scale = p.Results.vp_ratio .* data_block_norm ./ max(norm(profile_ratio_operator, "fro"), eps);
    A_system = [A_system; vp_ratio_scale .* profile_ratio_operator];
    b_system = [b_system; sparse(num_grid_cells, 1)];
end

% Solve the regularized inverse problem and reshape it back to the canonical
% Cartesian grid order [x y z].
[scan_1D, flag, relres, iter] = lsqr(A_system, b_system, p.Results.tol, p.Results.maxit);

% Keep orientation changes out of shared helpers. Plotting code can permute
% arrays locally if a display-oriented layout is needed.
scan = reshape(scan_1D,grid.sz);

inversion_info = struct( ...
    "controls", orderfields(struct(p.Results)), ...
    "solver", struct( ...
        "flag", flag, ...
        "relres", relres, ...
        "iter", iter, ...
        "n_data", numel(obs_values), ...
        "n_unknown", num_grid_cells));
end
