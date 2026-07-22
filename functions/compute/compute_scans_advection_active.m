function [active_fields, inversion_info] = compute_scans_advection_active(system, prior_field, varargin)
%COMPUTE_SCANS_ADVECTION_ACTIVE Solve active-domain 4D inversion with PCG.
%   This variant keeps large vectors in single precision and supports the
%   data, prior, and advection terms used by the full-US synthetic workflow.

p = inputParser;
addParameter(p, "advection", 1);
addParameter(p, "prior", 0.003);
addParameter(p, "tol", 1e-4);
addParameter(p, "maxit", 50);
addParameter(p, "nonnegative", true);
parse(p, varargin{:});

prior_2d = expand_active_prior(prior_field, system);
prior_scales = single(p.Results.prior .* double(system.data_block_norms) ./ sqrt(system.n_active));

advection_scales = zeros(max(system.n_time - 1, 0), 1, "single");
if p.Results.advection > 0
    for i_t = 1:system.n_time - 1
        reference_norm = max(mean(double(system.data_block_norms(i_t:i_t + 1))), 1);
        advection_scales(i_t) = single(p.Results.advection .* reference_norm ./ max(double(system.advection_base_norms(i_t)), eps));
    end
end

rhs = system.data_rhs;
system.data_rhs = [];
for i_t = 1:system.n_time
    rhs(:, i_t) = rhs(:, i_t) + prior_scales(i_t) ^ 2 .* prior_2d(:, i_t);
end
normal_diagonal = build_normal_diagonal(system, prior_scales, advection_scales);
system.data_diagonal = [];
normal_diagonal(normal_diagonal == 0) = 1;

apply_normal = @(x) apply_normal_operator(x, system, prior_scales, advection_scales);
apply_preconditioner = @(x) x ./ normal_diagonal;
x_initial = prior_2d(:);
rhs_vector = rhs(:);

[x, solver] = pcg_single(apply_normal, rhs_vector, p.Results.tol, p.Results.maxit, apply_preconditioner, x_initial);
if p.Results.nonnegative
    x = max(x, 0);
end
active_fields = reshape(x, system.n_active, system.n_time);

controls = struct(p.Results);
inversion_info = struct( ...
    "controls", orderfields(controls), ...
    "solver", orderfields(merge_structs(solver, struct( ...
        "n_time", system.n_time, ...
        "n_unknown_field", system.n_active, ...
        "n_unknown_total", numel(x)))), ...
    "scales", struct( ...
        "prior", prior_scales, ...
        "advection", advection_scales));
end

function [x, solver] = pcg_single(apply_normal, rhs, tol, maxit, apply_preconditioner, x_initial)

tol = tol(1);
maxit = max(1, round(maxit(1)));
x = single(x_initial);
rhs = single(rhs);
r = rhs - apply_normal(x);
z = apply_preconditioner(r);
p = z;
rz = sum(double(r) .* double(z), "all");
rhs_norm = max(norm(double(rhs)), eps);
relres = norm(double(r)) ./ rhs_norm;
flag = 1;
for iter = 1:maxit
    Ap = apply_normal(p);
    denom = sum(double(p) .* double(Ap), "all");
    if denom <= eps
        break
    end
    alpha = single(rz ./ denom);
    x = x + alpha .* p;
    r = r - alpha .* Ap;
    relres = norm(double(r)) ./ rhs_norm;
    if double(relres) <= double(tol)
        flag = 0;
        break
    end
    z = apply_preconditioner(r);
    rz_next = sum(double(r) .* double(z), "all");
    beta = single(rz_next ./ max(rz, eps));
    p = z + beta .* p;
    rz = rz_next;
end
solver = struct( ...
    "name", "active_single_pcg", ...
    "flag", flag, ...
    "relres", relres, ...
    "iterations", iter);
end

function y = apply_normal_operator(x, system, prior_scales, advection_scales)

x = reshape(single(x), system.n_active, system.n_time);
y = zeros(size(x), "single");
for i_t = 1:system.n_time
    data_projection = system.F * double(x(:, i_t));
    data_projection(~system.obs_mask(:, i_t)) = 0;
    y(:, i_t) = single(double(system.data_scales(i_t)) ^ 2 .* (system.F' * data_projection)) + ...
        prior_scales(i_t) ^ 2 .* x(:, i_t);
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
    y(:, i_t + 1) = y(:, i_t + 1) + single(accumarray(double(dest_id), double(weighted), [system.n_active 1]));
end
y = y(:);
end

function diagonal = build_normal_diagonal(system, prior_scales, advection_scales)

diagonal = system.data_diagonal;
for i_t = 1:system.n_time
    diagonal(:, i_t) = diagonal(:, i_t) + prior_scales(i_t) ^ 2;
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
        single(accumarray(double(step.dest_id_valid), double(coeff_squared), [system.n_active 1]));
end
diagonal = diagonal(:);
end

function prior_2d = expand_active_prior(prior_field, system)

if isscalar(prior_field)
    prior_2d = single(prior_field) .* ones(system.n_active, system.n_time, "single");
elseif isequal(size(prior_field), [system.n_active system.n_time])
    prior_2d = single(prior_field);
elseif numel(prior_field) == system.n_active * system.n_time
    prior_2d = reshape(single(prior_field), system.n_active, system.n_time);
elseif isequal(size(prior_field), [system.grid_sz system.n_time])
    prior_field_2d = reshape(prior_field, [], system.n_time);
    prior_2d = single(prior_field_2d(system.active_id, :));
else
    error("compute_scans_advection_active:InvalidPriorField", ...
        "prior_field must be scalar, n_active-by-n_time, or [g.sz n_time].")
end
end

function merged = merge_structs(primary, secondary)

merged = primary;
secondary_fields = fieldnames(secondary);
for i_field = 1:numel(secondary_fields)
    merged.(secondary_fields{i_field}) = secondary.(secondary_fields{i_field});
end
end
