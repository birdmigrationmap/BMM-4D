function [active_field, info] = solve_separable(observations, g, active_mask_3d, steps, level_reference, shape_patterns, z_agl, varargin)
%SOLVE_SEPARABLE Separable inversion: column amounts times vertical shape patterns.
%   Density on active cells is rho = B_t x with
%       rho(h, z, t) = x_1(h, t) s(h, z, t) + sum_{k=2..K} x_k(h, t) phi_k(z - ground(h)),
%       s = phi_1 + sum_j c_j(h, t) phi_j,
%   phi_k the columns of SHAPE_PATTERNS ([numel(z_agl) n], scaled so that
%   sum(phi_1) dz = 1 km: x_1 is the column amount in birds km^-2 and
%   x_k / x_1 a shape coefficient). s is the regional mean shape: phi_1
%   bent by the coefficients REGIONAL_SHAPE ([nxy J n_time], patterns
%   2..J+1, e.g. from regional_shape.m; empty = the night's mean shape).
%   x_2..x_K are local shape corrections around it. The ground is the DEM
%   low-pass filtered over GROUND_SMOOTHING_M; the full DEM still sets the
%   active cells (ACTIVE_MASK_3D) and the beam blockage (observations.F).
%   Outside the pattern height range phi keeps its end values (constant
%   extrapolation up to the top of the active layer); TOP_TAPER_M > 0
%   instead tapers it linearly to zero above the highest pattern height.
%
%   Objective (all times at once, term weights scaled as in
%   solve_advection_lean so the weights are comparable):
%   - data: ||F B x - d||^2 on observed rows, ROW_WEIGHTS optional;
%   - fine smoothness on x_1: SMOOTHNESS^2 ||W L_h x_1||^2, L_h the
%     horizontal graph Laplacian on active columns (land-water links cut);
%   - shape corrections k >= 2: SHAPE_PRIOR^2 / SHAPE_SD(k-1)^2 *
%     (W x_k)' (I + (SHAPE_LENGTH_M / dx)^2 L_h) (W x_k): small and
%     regionally smooth relative to the level, with a prior sd shrinking
%     with pattern order (penalty growing with order);
%   - optional COLUMN_PRIOR^2 ||W (x_1 - level)||^2 towards the level;
%   - advection on the reconstructed 3D field B x, relaxed where
%     LEVEL_REFERENCE changes (takeoff and landing), as in the lean solver;
%   - non-negativity: x_1 >= 0 by projection and NEGATIVE_PENALTY^2
%     ||min(B x, 0)||^2 on the 3D density, in projected FISTA.
%   W = median(level) / max(level, LEVEL_FLOOR) is the continuous density
%   weighting: structure costs relative to the local level, so dense areas
%   are not over-smoothed and empty areas carry little structure.
%   LEVEL_REFERENCE [n_active n_time] is a smooth 3D field (here the
%   Gaussian radar prior) whose column totals give the level.

p = inputParser;
addParameter(p, "n_pattern", 3);
addParameter(p, "shape_sd", [0.52 0.28 0.14 0.065 0.03]);
addParameter(p, "regional_shape", []);
addParameter(p, "smoothness", 0.3);
addParameter(p, "shape_prior", 0.3);
addParameter(p, "shape_length_m", 100000);
addParameter(p, "column_prior", 0);
addParameter(p, "level_floor", 10);
addParameter(p, "ground_smoothing_m", 30000);
addParameter(p, "top_taper_m", 0);
addParameter(p, "advection", 1);
addParameter(p, "advection_relax_exponent", 2);
addParameter(p, "advection_relax_floor", 5);
addParameter(p, "negative_penalty", 3);
addParameter(p, "row_weights", []);
addParameter(p, "tol", 1e-4);
addParameter(p, "maxit", 500);
addParameter(p, "projected_maxit", 100);
addParameter(p, "projected_tol", 1e-3);
parse(p, varargin{:});
opt = p.Results;
K = opt.n_pattern;
n_time = size(observations.values, 2);

%% Columns and shape basis ----------------------------------------------------
active_id = find(active_mask_3d(:));
n_active = numel(active_id);
[ix, iy, iz] = ind2sub(g.sz, active_id);
[column_xy, ~, column_of_cell] = unique(sub2ind(g.sz(1:2), ix, iy));
n_column = numel(column_xy);
ground_smooth = imgaussfilt(max(g.dem, 0), opt.ground_smoothing_m / mean([g.dx g.dy]));
z = squeeze(g.z);
height = z(iz) - ground_smooth(column_xy(column_of_cell));
if opt.top_taper_m > 0
    z_top = z_agl(end) + opt.top_taper_m;
    phi = interp1([z_agl(:); z_top], [shape_patterns; zeros(1, size(shape_patterns, 2))], min(max(height, z_agl(1)), z_top));
else
    phi = interp1(z_agl(:), shape_patterns, min(max(height, z_agl(1)), z_agl(end)));
end
% x_1 basis per time: the regional mean shape s(h, z, t)
B1 = cell(n_time, 1);
for i_t = 1:n_time
    s = phi(:, 1);
    for j = 1:size(opt.regional_shape, 2)
        s = s + phi(:, j + 1) .* opt.regional_shape(column_xy(column_of_cell), j, i_t);
    end
    B1{i_t} = sparse(1:n_active, column_of_cell, s, n_active, n_column);
end
% Local corrections x_2..x_K (fixed basis)
Br = sparse(repmat((1:n_active)', K - 1, 1), column_of_cell + n_column * (0:K - 2), ...
    reshape(phi(:, 2:K), [], 1), n_active, (K - 1) * n_column);

%% Data term --------------------------------------------------------------------
F = observations.F(:, active_id);
row_keep = any(F, 2);
F = F(row_keep, :);
values = observations.values(row_keep, :);
if ~isempty(opt.row_weights)
    row_scale = sqrt(double(opt.row_weights(row_keep)));
    keep_weighted = row_scale > 0;
    F = row_scale(keep_weighted) .* F(keep_weighted, :);
    values = row_scale(keep_weighted) .* values(keep_weighted, :);
end
obs_mask = full(values ~= 0);
data_scales = 1 ./ sqrt(max(sum(obs_mask, 1), 1));
% Block norms on the 3D operator, as in solve_advection_lean
data_block_norms = max(data_scales .* sqrt(full(sum(F .^ 2, 2))' * obs_mask), 1);
FB1 = cell(n_time, 1);
for i_t = 1:n_time
    FB1{i_t} = F * B1{i_t};
end
FBr = F * Br;
clear F
op = struct("B1", {B1}, "Br", Br, "FB1", {FB1}, "FBr", FBr, "n_column", n_column, "K", K, "n_time", n_time);
data_diagonal = zeros(K * n_column, n_time);
for i_t = 1:n_time
    data_diagonal(:, i_t) = [full((FB1{i_t} .^ 2)' * double(obs_mask(:, i_t))); full((FBr .^ 2)' * double(obs_mask(:, i_t)))];
end
data_diagonal = data_diagonal .* data_scales .^ 2;
rhs = apply_FBt(full(values), op) .* data_scales .^ 2;

%% Level weighting and horizontal operators ----------------------------------------
column_sum_op = sparse(column_of_cell, 1:n_active, g.dz / 1000, n_column, n_active);
level = max(column_sum_op * double(level_reference), opt.level_floor);   % birds km^-2
W = median(level(:)) ./ level;                                           % [n_column n_time]
column_mask = false([g.sz(1:2) 1]);
column_mask(column_xy) = true;
L_h = active_laplacian(column_mask, get_water_mask(g), 0.5) / 0.25;      % unit-weight graph Laplacian
Q = speye(n_column) + (opt.shape_length_m / mean([g.dx g.dy])) ^ 2 * L_h;

smoothness_scales2 = (opt.smoothness .* data_block_norms ./ norm(L_h, "fro")) .^ 2;
shape_scales2 = (opt.shape_prior .* data_block_norms ./ sqrt(n_column)) .^ 2 ./ reshape(opt.shape_sd(1:K - 1), [], 1) .^ 2;   % [K-1 n_time]
column_prior_scales2 = (opt.column_prior .* data_block_norms ./ sqrt(n_column)) .^ 2;
negative_scales2 = (opt.negative_penalty .* data_block_norms ./ sqrt(n_active)) .^ 2;

diagonal = data_diagonal;
first = 1:n_column;
diagonal(first, :) = diagonal(first, :) + full((L_h .^ 2)' * W .^ 2) .* smoothness_scales2 + W .^ 2 .* column_prior_scales2;
rhs(first, :) = rhs(first, :) + W .^ 2 .* level .* column_prior_scales2;
for k = 2:K
    block = (k - 1) * n_column + (1:n_column);
    diagonal(block, :) = diagonal(block, :) + W .^ 2 .* full(diag(Q)) .* shape_scales2(k - 1, :);
end

%% Advection links on the 3D field (as in solve_advection_lean) ------------------------
position = zeros(prod(g.sz), 1, "uint32");
position(active_id) = 1:n_active;
steps = steps(:);
advection_scales2 = zeros(max(n_time - 1, 0), 1);
advection_diagonal = zeros(n_active, n_time);
for i_t = 1:n_time - 1
    source = position(double(steps{i_t}.source_id_valid));
    dest = position(double(steps{i_t}.dest_id_valid));
    keep = source > 0 & dest > 0;
    coeff = steps{i_t}.coeff_valid;
    if ~isscalar(coeff)
        coeff = coeff(keep);
    end
    source = double(source(keep));
    dest = double(dest(keep));
    base_norm = sqrt(2 .* sum(double(coeff) .^ 2 .* ones(nnz(keep), 1)));
    if opt.advection_relax_exponent > 0
        reference_source = double(level_reference(source, i_t));
        reference_dest = double(level_reference(dest, i_t + 1));
        conserved = (min(reference_source, reference_dest) + opt.advection_relax_floor) ./ ...
            (max(reference_source, reference_dest) + opt.advection_relax_floor);
        coeff = double(coeff) .* conserved .^ opt.advection_relax_exponent;
    end
    coeff = double(coeff) .* ones(numel(source), 1);
    advection_scales2(i_t) = (opt.advection .* max(mean(data_block_norms(i_t:i_t + 1)), 1) ./ max(base_norm, eps)) .^ 2;
    w = advection_scales2(i_t) .* coeff .^ 2;
    advection_diagonal(:, i_t) = advection_diagonal(:, i_t) + accumarray(source, w, [n_active 1]);
    advection_diagonal(:, i_t + 1) = advection_diagonal(:, i_t + 1) + accumarray(dest, w, [n_active 1]);
    steps{i_t} = struct("source", source, "dest", dest, "coeff2", coeff .^ 2);
end
% Jacobi approximation of diag(B' A' A B) for the preconditioner
for i_t = 1:n_time
    diagonal(:, i_t) = diagonal(:, i_t) + [(B1{i_t} .^ 2)' * advection_diagonal(:, i_t); (Br .^ 2)' * advection_diagonal(:, i_t)];
end
ridge = 1e-8 * max(diagonal(:));
diagonal = diagonal + ridge;
clear advection_diagonal data_diagonal

op.obs_mask = obs_mask;
op.data_scales2 = data_scales .^ 2;
op.L_h = L_h;
op.Q = Q;
op.W = W;
op.smoothness_scales2 = smoothness_scales2;
op.shape_scales2 = shape_scales2;
op.column_prior_scales2 = column_prior_scales2;
op.advection_scales2 = advection_scales2;
op.steps = steps;
op.ridge = ridge;

%% Preconditioned conjugate gradient --------------------------------------------
x = zeros(K * n_column, n_time);
inverse_diagonal = 1 ./ diagonal;
r = rhs - normal_operator(x, op);
zr = r .* inverse_diagonal;
d = zr;
rz = sum(r .* zr, "all");
rhs_norm = max(norm(rhs(:)), eps);
relres = norm(r(:)) / rhs_norm;
iter = 0;
for iter = 1:opt.maxit
    Ad = normal_operator(d, op);
    alpha = rz / sum(d .* Ad, "all");
    x = x + alpha .* d;
    r = r - alpha .* Ad;
    relres = norm(r(:)) / rhs_norm;
    if relres <= opt.tol
        break
    end
    zr = r .* inverse_diagonal;
    rz_next = sum(r .* zr, "all");
    d = zr + (rz_next / rz) .* d;
    rz = rz_next;
end
solver = struct("name", "separable_pcg", "iterations", iter, "relres", relres, ...
    "negative_share_pcg", mean(apply_B(x, op) < 0, "all"));

%% Projected FISTA: x_1 >= 0 and a penalty on negative 3D densities ------------------
% Largest eigenvalue of D^-1 (H + mu B'B) by power iteration
v = randn(size(x));
v = v ./ norm(v(:));
lambda_max = 0;
for i_power = 1:15
    w = (normal_operator(v, op) + apply_Bt(apply_B(v, op) .* negative_scales2, op)) .* inverse_diagonal;
    lambda_max = norm(w(:));
    v = w ./ lambda_max;
end
step = 1 / (1.05 * lambda_max);
x(first, :) = max(x(first, :), 0);
y = x;
t_momentum = 1;
for i_projected = 1:opt.projected_maxit
    gradient = normal_operator(y, op) - rhs + apply_Bt(min(apply_B(y, op), 0) .* negative_scales2, op);
    x_next = y - step .* gradient .* inverse_diagonal;
    x_next(first, :) = max(x_next(first, :), 0);
    relative_step = norm(x_next(:) - x(:)) / max(norm(x(:)), eps);
    t_next = (1 + sqrt(1 + 4 * t_momentum ^ 2)) / 2;
    y = x_next + ((t_momentum - 1) / t_next) .* (x_next - x);
    x = x_next;
    t_momentum = t_next;
    if relative_step <= opt.projected_tol
        break
    end
end
solver.name = "separable_pcg_projected_fista";
solver.projected_iterations = i_projected;
solver.projected_relative_step = relative_step;

%% Density on active cells ---------------------------------------------------------
rho = apply_B(x, op);
solver.negative_mass_share = -sum(min(rho, 0), "all") / max(sum(max(rho, 0), "all"), eps);
active_field = single(max(rho, 0));
info = struct("controls", rmfield(opt, ["row_weights" "regional_shape"]), "solver", solver, ...
    "pattern_amounts", single(reshape(x, n_column, K, n_time)), "column_xy", column_xy, "level", single(level));
end

function rho = apply_B(X, op)
% 3D density on active cells from the pattern amounts, per time.
n = op.n_column;
rho = op.Br * X(n + 1:end, :);
for i_t = 1:op.n_time
    rho(:, i_t) = rho(:, i_t) + op.B1{i_t} * X(1:n, i_t);
end
end

function Y = apply_Bt(R, op)
Y = zeros(op.K * op.n_column, op.n_time);
for i_t = 1:op.n_time
    Y(1:op.n_column, i_t) = op.B1{i_t}' * R(:, i_t);
end
Y(op.n_column + 1:end, :) = op.Br' * R;
end

function P = apply_FB(X, op)
n = op.n_column;
P = op.FBr * X(n + 1:end, :);
for i_t = 1:op.n_time
    P(:, i_t) = P(:, i_t) + op.FB1{i_t} * X(1:n, i_t);
end
end

function Y = apply_FBt(P, op)
Y = zeros(op.K * op.n_column, op.n_time);
for i_t = 1:op.n_time
    Y(1:op.n_column, i_t) = op.FB1{i_t}' * P(:, i_t);
end
Y(op.n_column + 1:end, :) = op.FBr' * P;
end

function Y = normal_operator(X, op)
% Data, fine smoothness on x_1, shape priors on x_k, advection on B x.
n = op.n_column;
projection = apply_FB(X, op);
projection(~op.obs_mask) = 0;
Y = apply_FBt(projection .* op.data_scales2, op) + op.ridge .* X;
first = 1:n;
X1 = X(first, :);
Y(first, :) = Y(first, :) + op.L_h' * ((op.L_h * X1) .* op.W .^ 2) .* op.smoothness_scales2 ...
    + X1 .* op.W .^ 2 .* op.column_prior_scales2;
for k = 2:op.K
    block = (k - 1) * n + (1:n);
    Y(block, :) = Y(block, :) + op.W .* (op.Q * (op.W .* X(block, :))) .* op.shape_scales2(k - 1, :);
end
if op.n_time > 1
    rho = apply_B(X, op);
    Y_rho = zeros(size(rho));
    for i_t = 1:op.n_time - 1
        step = op.steps{i_t};
        weighted = op.advection_scales2(i_t) .* step.coeff2 .* (rho(step.dest, i_t + 1) - rho(step.source, i_t));
        Y_rho(:, i_t) = Y_rho(:, i_t) - accumarray(step.source, weighted, [size(rho, 1) 1]);
        Y_rho(:, i_t + 1) = Y_rho(:, i_t + 1) + accumarray(step.dest, weighted, [size(rho, 1) 1]);
    end
    Y = Y + apply_Bt(Y_rho, op);
end
end
