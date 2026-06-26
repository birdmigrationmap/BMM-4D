function [fields, inversion_info] = compute_scans_advection(system, prior_field, varargin)
%COMPUTE_SCANS_ADVECTION Reconstruct a 4D field with matrix-free PCG.

p = inputParser;
addParameter(p, "advection", 1);
addParameter(p, "prior", 0.003);
addParameter(p, "smoothness", 0.1);
addParameter(p, "solver", "pcg");
addParameter(p, "projected_maxit", 100);
addParameter(p, "projected_tol", 1e-4);
addParameter(p, "projected_step", 0.25);
addParameter(p, "tol", 1e-6);
addParameter(p, "maxit", 100);
parse(p, varargin{:});

prior_4d = expand_prior_field(prior_field, system.grid_sz, system.n_time);
prior_scales = p.Results.prior .* system.data_block_norms ./ sqrt(system.nxyz);

if p.Results.smoothness > 0
    assert(~isempty(system.laplacian), "compute_scans_advection:MissingSmoothnessOperator", ...
        "Prepare the system with smoothness=true before using spatial smoothing.")
    smoothness_scales = p.Results.smoothness .* system.data_block_norms ./ max(system.laplacian_norm, eps);
else
    smoothness_scales = zeros(1, system.n_time);
end

advection_scales = zeros(max(system.n_time - 1, 0), 1);
if p.Results.advection > 0
    for i_t = 1:system.n_time - 1
        reference_norm = max(mean(system.data_block_norms(i_t:i_t + 1)), 1);
        advection_scales(i_t) = p.Results.advection .* reference_norm ./ max(system.advection_base_norms(i_t), eps);
    end
end

rhs = system.data_rhs;
for i_t = 1:system.n_time
    rhs(:, i_t) = rhs(:, i_t) + prior_scales(i_t) ^ 2 .* reshape(prior_4d(:, :, :, i_t), [], 1);
end

normal_diagonal = build_normal_diagonal(system, prior_scales, smoothness_scales, advection_scales);
normal_diagonal(normal_diagonal == 0) = 1;
apply_normal = @(x) apply_normal_operator(x, system, prior_scales, smoothness_scales, advection_scales);
apply_preconditioner = @(x) x ./ normal_diagonal;

x_initial = prior_4d(:);
rhs_vector = rhs(:);
solver_name = string(p.Results.solver);
if solver_name == "pcg"
    [x, flag, relres, iter] = pcg( ...
        apply_normal, rhs_vector, p.Results.tol, p.Results.maxit, ...
        apply_preconditioner, [], x_initial);
    solver = struct( ...
        "name", "matrix_free_pcg_normal", ...
        "flag", flag, ...
        "relres", relres, ...
        "iterations", iter);
elseif solver_name == "projected_fista"
    [x, solver] = projected_fista( ...
        apply_normal, apply_preconditioner, rhs_vector, max(x_initial, 0), ...
        p.Results.projected_maxit, p.Results.projected_tol, p.Results.projected_step, ...
        "matrix_free_projected_fista");
elseif solver_name == "pcg_projected_fista"
    [x_pcg, pcg_flag, pcg_relres, pcg_iter] = pcg( ...
        apply_normal, rhs_vector, p.Results.tol, p.Results.maxit, ...
        apply_preconditioner, [], x_initial);
    [x, projected_solver] = projected_fista( ...
        apply_normal, apply_preconditioner, rhs_vector, max(x_pcg, 0), ...
        p.Results.projected_maxit, p.Results.projected_tol, p.Results.projected_step, ...
        "matrix_free_pcg_projected_fista");
    solver = projected_solver;
    solver.pcg_flag = pcg_flag;
    solver.pcg_relres = pcg_relres;
    solver.pcg_iterations = pcg_iter;
else
    error("compute_scans_advection:UnknownSolver", ...
        "Unknown solver '%s'. Use 'pcg', 'projected_fista', or 'pcg_projected_fista'.", solver_name);
end
fields = reshape(x, [system.grid_sz system.n_time]);

objective = compute_objective(fields, prior_4d, system, prior_scales, smoothness_scales, advection_scales);
inversion_info = struct( ...
    "controls", orderfields(struct(p.Results)), ...
    "solver", orderfields(merge_structs( ...
        solver, struct( ...
        "n_time", system.n_time, ...
        "n_unknown_field", system.nxyz, ...
        "n_unknown_total", numel(x)))), ...
    "objective", objective);
end

function [x, solver] = projected_fista(apply_normal, apply_preconditioner, rhs, x_initial, maxit, tol, step_scale, solver_name)

maxit = max(1, round(maxit));
tol = tol(1);
step_scale = step_scale(1);
x = max(x_initial, 0);
y = x;
t = 1;
rel_step = inf;
projected_gradient_norm = inf;
for iter = 1:maxit
    gradient = apply_normal(y) - rhs;
    x_next = max(y - step_scale .* apply_preconditioner(gradient), 0);
    rel_step = norm(x_next - x) / max(norm(x), 1);
    projected_gradient_norm = norm(y - x_next) / max(step_scale, eps);
    if rel_step <= tol
        x = x_next;
        break
    end

    t_next = 0.5 .* (1 + sqrt(1 + 4 .* t .^ 2));
    y = x_next + ((t - 1) ./ t_next) .* (x_next - x);
    x = x_next;
    t = t_next;
end

residual = apply_normal(x) - rhs;
solver = struct( ...
    "name", solver_name, ...
    "flag", double(rel_step > tol), ...
    "relres", norm(residual) ./ max(norm(rhs), eps), ...
    "iterations", iter, ...
    "projected_rel_step", rel_step, ...
    "projected_gradient_norm", projected_gradient_norm, ...
    "projected_step", step_scale, ...
    "negative_fraction", mean(x < 0));
end

function merged = merge_structs(primary, secondary)

merged = primary;
secondary_fields = fieldnames(secondary);
for i_field = 1:numel(secondary_fields)
    merged.(secondary_fields{i_field}) = secondary.(secondary_fields{i_field});
end
end

function y = apply_normal_operator(x, system, prior_scales, smoothness_scales, advection_scales)

x = reshape(x, system.nxyz, system.n_time);
y = zeros(size(x));
for i_t = 1:system.n_time
    data_projection = system.F * x(:, i_t);
    data_projection(~system.obs_mask(:, i_t)) = 0;
    y(:, i_t) = system.data_scales(i_t) ^ 2 .* (system.F' * data_projection) + ...
        prior_scales(i_t) ^ 2 .* x(:, i_t);
    if smoothness_scales(i_t) > 0
        y(:, i_t) = y(:, i_t) + smoothness_scales(i_t) ^ 2 .* ...
            (system.laplacian' * (system.laplacian * x(:, i_t)));
    end
end

for i_t = 1:system.n_time - 1
    if advection_scales(i_t) == 0
        continue
    end
    step = system.steps{i_t};
    source_id = step.source_id_valid;
    dest_id = step.dest_id_valid;
    coeff = step.coeff_valid;
    residual = coeff .* (x(dest_id, i_t + 1) - x(source_id, i_t));
    weighted = advection_scales(i_t) ^ 2 .* coeff .* residual;
    y(source_id, i_t) = y(source_id, i_t) - weighted;
    y(:, i_t + 1) = y(:, i_t + 1) + accumarray(double(dest_id), weighted, [system.nxyz 1]);
end
y = y(:);
end

function diagonal = build_normal_diagonal(system, prior_scales, smoothness_scales, advection_scales)

diagonal = system.data_diagonal;
for i_t = 1:system.n_time
    diagonal(:, i_t) = diagonal(:, i_t) + prior_scales(i_t) ^ 2;
    if smoothness_scales(i_t) > 0
        diagonal(:, i_t) = diagonal(:, i_t) + ...
            smoothness_scales(i_t) ^ 2 .* system.laplacian_normal_diagonal;
    end
end

for i_t = 1:system.n_time - 1
    if advection_scales(i_t) == 0
        continue
    end
    step = system.steps{i_t};
    coeff_squared = advection_scales(i_t) ^ 2 .* step.coeff_valid .^ 2;
    if isscalar(coeff_squared)
        coeff_squared = repmat(coeff_squared, numel(step.source_id_valid), 1);
    end
    diagonal(step.source_id_valid, i_t) = diagonal(step.source_id_valid, i_t) + coeff_squared;
    diagonal(:, i_t + 1) = diagonal(:, i_t + 1) + ...
        accumarray(double(step.dest_id_valid), coeff_squared, [system.nxyz 1]);
end
diagonal = diagonal(:);
end

function objective = compute_objective(fields, prior_4d, system, prior_scales, smoothness_scales, advection_scales)

fields_2d = reshape(fields, system.nxyz, system.n_time);
prior_2d = reshape(prior_4d, system.nxyz, system.n_time);
data_value = 0;
prior_value = 0;
smoothness_value = 0;
advection_value = 0;
for i_t = 1:system.n_time
    residual = system.F * fields_2d(:, i_t) - full(system.obs_values(:, i_t));
    residual(~system.obs_mask(:, i_t)) = 0;
    data_value = data_value + norm(system.data_scales(i_t) .* residual) ^ 2;
    prior_value = prior_value + norm(prior_scales(i_t) .* (fields_2d(:, i_t) - prior_2d(:, i_t))) ^ 2;
    if smoothness_scales(i_t) > 0
        smoothness_value = smoothness_value + ...
            norm(smoothness_scales(i_t) .* (system.laplacian * fields_2d(:, i_t))) ^ 2;
    end
end
for i_t = 1:system.n_time - 1
    if advection_scales(i_t) > 0
        step = system.steps{i_t};
        residual = step.coeff_valid .* ( ...
            fields_2d(step.dest_id_valid, i_t + 1) - fields_2d(step.source_id_valid, i_t));
        advection_value = advection_value + norm(advection_scales(i_t) .* residual) ^ 2;
    end
end
objective = struct( ...
    "data", data_value, ...
    "prior", prior_value, ...
    "smoothness", smoothness_value, ...
    "advection", advection_value, ...
    "total", data_value + prior_value + smoothness_value + advection_value);
end

function prior_4d = expand_prior_field(prior_field, grid_sz, n_time)

if isscalar(prior_field)
    prior_4d = prior_field .* ones([grid_sz n_time]);
elseif matches_size(prior_field, [grid_sz n_time])
    prior_4d = reshape(prior_field, [grid_sz n_time]);
elseif matches_size(prior_field, grid_sz)
    prior_4d = repmat(reshape(prior_field, grid_sz), 1, 1, 1, n_time);
elseif isvector(prior_field) && numel(prior_field) == grid_sz(3)
    prior_4d = repeat_profile(prior_field(:), grid_sz, n_time);
elseif ismatrix(prior_field) && size(prior_field, 1) == grid_sz(3) && size(prior_field, 2) == n_time
    prior_4d = zeros([grid_sz n_time]);
    for i_t = 1:n_time
        prior_4d(:, :, :, i_t) = repeat_profile(prior_field(:, i_t), grid_sz, 1);
    end
else
    error("compute_scans_advection:InvalidPriorField", ...
        "prior_field must be scalar, nz-by-1, nz-by-n_time, grid.sz, or [grid.sz n_time].");
end
end

function field = repeat_profile(profile, grid_sz, n_time)

nxy = prod(grid_sz(1:2));
profile_field = reshape(repelem(profile(:), nxy), grid_sz);
field = repmat(profile_field, 1, 1, 1, n_time);
end

function tf = matches_size(value, target_size)

value_size = size(value);
if numel(value_size) < numel(target_size)
    value_size(end + 1:numel(target_size)) = 1;
elseif numel(value_size) > numel(target_size)
    target_size(end + 1:numel(value_size)) = 1;
end
tf = isequal(value_size, target_size);
end
