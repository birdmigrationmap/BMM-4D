function [active_field, info] = solve_advection_lean(observations, g, active_mask_3d, steps, prior_active, varargin)
%SOLVE_ADVECTION_LEAN Joint advection inversion on active cells, batched.
%   Same objective and term weights as the paper solver (data, prior,
%   3D smoothness with land-water cut, advection), restricted to the active
%   cells (e.g. radar coverage, ground to a flight layer) and optionally in
%   single precision. All times are applied at once (F * X). With every cell
%   active and store_class "double" it follows the paper solver's PCG
%   iterations to rounding (Key West: 3e-12 after 50 iterations); it returns
%   the last iterate, where the paper's MATLAB pcg returned the one with the
%   lowest residual (MIGRATION.md, step 3).
%
%   STEPS are full-grid advection steps (build_radar_advection_steps).
%   PRIOR_ACTIVE is [n_active n_time]. Returns ACTIVE_FIELD [n_active n_time].
%
%   ADVECTION_RELAX_EXPONENT > 0 weakens each advection link where the
%   density prior changes along it (takeoff after sunset, landing before
%   sunrise), where birds are not conserved: the link weight is
%   ((min + floor) / (max + floor))^exponent of the prior at its source and
%   destination, with floor ADVECTION_RELAX_FLOOR (birds km^-3).
%
%   NONNEGATIVE true continues from the CG solution, clipped at zero, with
%   projected FISTA (accelerated projected gradient, diagonal
%   preconditioner) so that every density is >= 0. The step is
%   1 / lambda_max(D^-1 H), estimated by power iteration; it stops after
%   PROJECTED_MAXIT iterations or when the relative step is below
%   PROJECTED_TOL.
%
%   Options for the radar-imprint study (all off by default):
%   - ROW_WEIGHTS [n_rows 1]: weight of each observation row (soft
%     replacement of hard data cuts); rows and values are scaled by sqrt(w).
%   - SHAPE (lambda): separable-model penalty lambda^2 ||x - p * sum_z(x)||^2
%     per column and time, with p the prior's normalized vertical profile
%     of that column: the column total is free, its vertical distribution
%     is held close to the prior's shape (ties unobserved low layers to the
%     observed column).
%   - PRIOR_CORRELATION_M (l): the prior term becomes the Gaussian-Markov
%     form (x - prior)' (I + (l/dx)^2 L_h) (x - prior), L_h the horizontal
%     graph Laplacian on active cells: deviations from the prior are
%     horizontally correlated over about l (first-order Matern / SPDE),
%     so a correction made by the data near a radar spreads instead of
%     stopping at the beam edge.
%   - PRIOR_CORRELATION_VERTICAL_M (l_v): adds (l_v/dz)^2 L_v (vertical graph
%     Laplacian) to that precision: a 3D Gaussian-Markov prior. With
%     smoothness 0 the prior term alone is then a Gaussian-process prior
%     (kriging) around PRIOR_ACTIVE as its mean.
%   - ADVECTION_REFERENCE [n_active n_time]: field used to relax the
%     advection links (default PRIOR_ACTIVE), so a zero or flat GP mean
%     can be used without losing the takeoff/landing relaxation.

p = inputParser;
addParameter(p, "prior", 0.1);
addParameter(p, "smoothness", 0.1);
addParameter(p, "advection", 1);
addParameter(p, "coeff_horver", 0);
addParameter(p, "tol", 1e-6);
addParameter(p, "maxit", 100);
addParameter(p, "store_class", "single");
addParameter(p, "min_row_mass", 0);
addParameter(p, "x_initial", []);
addParameter(p, "advection_relax_exponent", 0);
addParameter(p, "advection_relax_floor", 5);
addParameter(p, "nonnegative", false);
addParameter(p, "row_weights", []);
addParameter(p, "shape", 0);
addParameter(p, "prior_correlation_m", 0);
addParameter(p, "prior_correlation_vertical_m", 0);
addParameter(p, "advection_reference", []);
addParameter(p, "projected_maxit", 100);
addParameter(p, "projected_tol", 1e-3);
parse(p, varargin{:});
opt = p.Results;
cls = char(opt.store_class);
% One advection step per interval between observation times
if numel(steps) ~= max(size(observations.values, 2) - 1, 0)
    error("solve_advection_lean:InvalidStepCount", "Expected %d advection steps, got %d.", ...
        max(size(observations.values, 2) - 1, 0), numel(steps))
end

%% Data term on active columns ---------------------------------------------
active_id = find(active_mask_3d(:));
n_active = numel(active_id);
F = observations.F(:, active_id);
row_keep = full(sum(F, 2)) >= opt.min_row_mass & any(F, 2);
F = F(row_keep, :);
values = observations.values(row_keep, :);
if ~isempty(opt.row_weights)
    % Soft data weights: scale rows and values by sqrt(w); w = 0 removes a row
    row_scale = sqrt(double(opt.row_weights(row_keep)));
    F = spdiags(row_scale, 0, numel(row_scale), numel(row_scale)) * F;
    values = spdiags(row_scale, 0, numel(row_scale), numel(row_scale)) * values;
    keep_weighted = row_scale > 0;
    F = F(keep_weighted, :);
    values = values(keep_weighted, :);
end
n_time = size(values, 2);
obs_mask = full(values ~= 0);
data_scales = 1 ./ sqrt(max(sum(obs_mask, 1), 1));
F_squared = F .^ 2;
row_norm_squared = full(sum(F_squared, 2));
data_block_norms = max(data_scales .* sqrt(row_norm_squared' * obs_mask), 1);
data_diagonal = full(F_squared' * double(obs_mask)) .* data_scales .^ 2;
data_rhs = full(F' * values) .* data_scales .^ 2;
clear F_squared

%% Smoothness and advection operators ---------------------------------------
L = active_laplacian(active_mask_3d, get_water_mask(g), opt.coeff_horver);
laplacian_norm = norm(L, "fro");
laplacian_diagonal = full(sum(L .^ 2, 1))';

% Separable shape term: column membership and the prior's vertical shape
[ix, iy, ~] = ind2sub(g.sz, active_id);
[~, ~, column_of_cell] = unique(sub2ind(g.sz(1:2), ix, iy));
n_column = max(column_of_cell);
column_sum_op = sparse(column_of_cell, 1:n_active, 1, n_column, n_active);
shape_profile = [];
shape_diagonal = 0;
if opt.shape > 0
    prior_column = column_sum_op * double(prior_active);
    shape_profile = double(prior_active) ./ max(prior_column(column_of_cell, :), eps);
    cells_per_column = full(sum(column_sum_op, 2));
    empty = prior_column(column_of_cell, :) <= 0;
    uniform = 1 ./ cells_per_column(column_of_cell);
    shape_profile(empty) = uniform(mod(find(empty) - 1, n_active) + 1);
    % diag(S'S)_i = 1 - 2 p_i + sum_col p^2
    shape_diagonal = 1 - 2 * shape_profile + (column_sum_op' * (column_sum_op * shape_profile .^ 2));
end
% Correlated prior: graph Laplacian of horizontal links (unit weights)
prior_precision = [];
if opt.prior_correlation_m > 0 || opt.prior_correlation_vertical_m > 0
    prior_precision = speye(n_active);
    if opt.prior_correlation_m > 0
        L_horizontal = active_laplacian(active_mask_3d, get_water_mask(g), 0.5) / 0.25;
        prior_precision = prior_precision + (opt.prior_correlation_m / mean([g.dx g.dy])) ^ 2 * L_horizontal;
        clear L_horizontal
    end
    if opt.prior_correlation_vertical_m > 0
        % coeff_horver = -1: horizontal weight 0, vertical weight 1/2
        L_vertical = active_laplacian(active_mask_3d, get_water_mask(g), -1) / 0.5;
        prior_precision = prior_precision + (opt.prior_correlation_vertical_m / g.dz) ^ 2 * L_vertical;
        clear L_vertical
    end
end
advection_reference = prior_active;
if ~isempty(opt.advection_reference)
    advection_reference = opt.advection_reference;
end

position = zeros(prod(g.sz), 1, "uint32");
position(active_id) = 1:n_active;
steps = steps(:);
advection_base_norms = zeros(max(n_time - 1, 0), 1);
for i_t = 1:n_time - 1
    source = position(double(steps{i_t}.source_id_valid));
    dest = position(double(steps{i_t}.dest_id_valid));
    keep = source > 0 & dest > 0;
    coeff = steps{i_t}.coeff_valid;
    if ~isscalar(coeff)
        coeff = coeff(keep);
    end
    % uint32 link indices: half the memory of double on large domains
    source = source(keep);
    dest = dest(keep);
    advection_base_norms(i_t) = sqrt(2 .* sum(double(coeff) .^ 2 .* ones(nnz(keep), 1)));
    if opt.advection_relax_exponent > 0
        prior_source = double(advection_reference(source, i_t));
        prior_dest = double(advection_reference(dest, i_t + 1));
        conserved = (min(prior_source, prior_dest) + opt.advection_relax_floor) ./ ...
            (max(prior_source, prior_dest) + opt.advection_relax_floor);
        coeff = double(coeff) .* conserved .^ opt.advection_relax_exponent;
    end
    steps{i_t} = struct("source", source, "dest", dest, "coeff", cast(coeff, cls));
end

%% Term weights, as in the paper solver ----------------------------------------
prior_scales = opt.prior .* data_block_norms ./ sqrt(n_active);
smoothness_scales = opt.smoothness .* data_block_norms ./ max(laplacian_norm, eps);
advection_scales = zeros(max(n_time - 1, 0), 1);
if opt.advection > 0
    for i_t = 1:n_time - 1
        advection_scales(i_t) = opt.advection .* max(mean(data_block_norms(i_t:i_t + 1)), 1) ./ max(advection_base_norms(i_t), eps);
    end
end

shape_scales = opt.shape .* data_block_norms ./ sqrt(n_active);
if isempty(prior_precision)
    rhs = data_rhs + double(prior_active) .* prior_scales .^ 2;
    diagonal = data_diagonal + prior_scales .^ 2 + laplacian_diagonal .* smoothness_scales .^ 2;
else
    rhs = data_rhs + (prior_precision * double(prior_active)) .* prior_scales .^ 2;
    diagonal = data_diagonal + full(diag(prior_precision)) .* prior_scales .^ 2 + laplacian_diagonal .* smoothness_scales .^ 2;
end
if opt.shape > 0
    diagonal = diagonal + shape_diagonal .* shape_scales .^ 2;
end
for i_t = 1:n_time - 1
    w = advection_scales(i_t) ^ 2 .* double(steps{i_t}.coeff) .^ 2 .* ones(numel(steps{i_t}.source), 1);
    diagonal(:, i_t) = diagonal(:, i_t) + accumarray(steps{i_t}.source, w, [n_active 1]);
    diagonal(:, i_t + 1) = diagonal(:, i_t + 1) + accumarray(steps{i_t}.dest, w, [n_active 1]);
end
diagonal(diagonal == 0) = 1;
clear data_rhs data_diagonal

op = struct("F", cast(F, cls), "L", cast(L, cls), "obs_mask", obs_mask, ...
    "data_scales2", cast(data_scales .^ 2, cls), "prior_scales2", cast(prior_scales .^ 2, cls), ...
    "smoothness_scales2", cast(smoothness_scales .^ 2, cls), "advection_scales2", cast(advection_scales .^ 2, cls), ...
    "steps", {steps}, "n_active", n_active, "n_time", n_time, ...
    "shape_scales2", cast(shape_scales .^ 2, cls), "shape_profile", cast(shape_profile, cls), ...
    "column_sum_op", cast(column_sum_op, cls), "column_of_cell", column_of_cell, ...
    "prior_precision", cast(prior_precision, cls));
clear shape_profile prior_precision
clear F L

%% Preconditioned conjugate gradient ------------------------------------------
x = cast(prior_active, cls);
if ~isempty(opt.x_initial)
    x = cast(opt.x_initial, cls);
end
rhs = cast(rhs, cls);
inverse_diagonal = cast(1 ./ diagonal, cls);
clear diagonal
r = rhs - normal_operator(x, op);
z = r .* inverse_diagonal;
d = z;
rz = sum(double(r) .* double(z), "all");
rhs_norm = max(norm(double(rhs(:))), eps);
relres = norm(double(r(:))) / rhs_norm;
flag = 1;
max_iteration = opt.maxit;
if relres <= opt.tol
    % The initial field already solves the system (e.g. prior only, no data)
    flag = 0;
    max_iteration = 0;
end
for iter = 1:max_iteration
    Ad = normal_operator(d, op);
    alpha = rz / sum(double(d) .* double(Ad), "all");
    x = x + cast(alpha, cls) .* d;
    r = r - cast(alpha, cls) .* Ad;
    relres = norm(double(r(:))) / rhs_norm;
    if relres <= opt.tol
        flag = 0;
        break
    end
    z = r .* inverse_diagonal;
    rz_next = sum(double(r) .* double(z), "all");
    d = z + cast(rz_next / rz, cls) .* d;
    rz = rz_next;
end
if max_iteration == 0
    iter = 0;   % an empty for loop leaves iter empty
end
solver = struct("name", "lean_batched_pcg", "flag", flag, "relres", relres, "iterations", iter, ...
    "n_active", n_active, "n_time", n_time, "n_rows", size(op.F, 1), "store_class", string(cls), ...
    "negative_share_pcg", mean(x(:) < 0));

%% Projected FISTA for non-negative densities ---------------------------------
if opt.nonnegative
    clear r z d Ad
    % Largest eigenvalue of D^-1 H by power iteration (similar to the
    % symmetric D^-1/2 H D^-1/2, so the iteration converges to it)
    v = randn(size(x), cls);
    v = v ./ norm(double(v(:)));
    lambda_max = 0;
    for i_power = 1:15
        w = normal_operator(v, op) .* inverse_diagonal;
        lambda_max = norm(double(w(:)));
        v = w ./ lambda_max;
    end
    clear v w
    step = cast(1 / (1.05 * lambda_max), cls);
    x = max(x, 0);
    y = x;
    t_momentum = 1;
    for i_projected = 1:opt.projected_maxit
        gradient = normal_operator(y, op) - rhs;
        x_next = max(y - step .* gradient .* inverse_diagonal, 0);
        relative_step = norm(double(x_next(:) - x(:))) / max(norm(double(x(:))), eps);
        t_next = (1 + sqrt(1 + 4 * t_momentum ^ 2)) / 2;
        y = x_next + cast((t_momentum - 1) / t_next, cls) .* (x_next - x);
        x = x_next;
        t_momentum = t_next;
        if relative_step <= opt.projected_tol
            break
        end
    end
    clear y gradient x_next
    projected_residual = norm(double(reshape(normal_operator(x, op) - rhs, [], 1))) / rhs_norm;
    solver.name = "lean_batched_pcg_projected_fista";
    solver.projected_iterations = i_projected;
    solver.projected_relative_step = relative_step;
    solver.projected_relres = projected_residual;
    solver.projected_lambda_max = lambda_max;
end
active_field = x;
info = struct("controls", rmfield(opt, ["x_initial" "row_weights" "advection_reference"]), "solver", solver);
end

function Y = normal_operator(X, op)
% Data, prior and smoothness for all times at once, then advection links.
projection = op.F * X;
projection(~op.obs_mask) = 0;
if isempty(op.prior_precision)
    prior_part = X;
else
    prior_part = op.prior_precision * X;
end
Y = (op.F' * projection) .* op.data_scales2 + prior_part .* op.prior_scales2 + (op.L' * (op.L * X)) .* op.smoothness_scales2;
if ~isempty(op.shape_profile)
    % S x = x - p * column_sum(x);  S' v = v - column_sum(p .* v)
    shape_residual = X - op.shape_profile .* (op.column_sum_op' * (op.column_sum_op * X));
    Y = Y + (shape_residual - op.column_sum_op' * (op.column_sum_op * (op.shape_profile .* shape_residual))) .* op.shape_scales2;
end
for i_t = 1:op.n_time - 1
    step = op.steps{i_t};
    weighted = op.advection_scales2(i_t) .* step.coeff .^ 2 .* (X(step.dest, i_t + 1) - X(step.source, i_t));
    Y(:, i_t) = Y(:, i_t) - accumarray(step.source, weighted, [op.n_active 1]);
    Y(:, i_t + 1) = Y(:, i_t + 1) + accumarray(step.dest, weighted, [op.n_active 1]);
end
end
