function [fields, info] = solve_advection_fast(system, prior_field, varargin)
%SOLVE_ADVECTION_FAST Joint advection inversion with batched operators.
%   Prod copy of the PCG path of functions/compute/compute_scans_advection.m
%   (same scales, preconditioner and start), with the normal operator
%   applied to all times at once: F, F' and the Laplacian are the same for
%   every time, so F * X with one column per time reads each sparse operator
%   once per iteration instead of once per time. SYSTEM comes from
%   prepare_scans_advection; PRIOR_FIELD is [g.sz n_time].

p = inputParser;
addParameter(p, "advection", 1);
addParameter(p, "prior", 0.003);
addParameter(p, "smoothness", 0.1);
addParameter(p, "tol", 1e-6);
addParameter(p, "maxit", 100);
parse(p, varargin{:});

n_time = system.n_time;
nxyz = system.nxyz;
prior_2d = reshape(prior_field, nxyz, n_time);

% Term weights, as in compute_scans_advection.
prior_scales = p.Results.prior .* system.data_block_norms ./ sqrt(nxyz);
smoothness_scales = zeros(1, n_time);
if p.Results.smoothness > 0
    smoothness_scales = p.Results.smoothness .* system.data_block_norms ./ max(system.laplacian_norm, eps);
end
advection_scales = zeros(max(n_time - 1, 0), 1);
if p.Results.advection > 0
    for i_t = 1:n_time - 1
        reference_norm = max(mean(system.data_block_norms(i_t:i_t + 1)), 1);
        advection_scales(i_t) = p.Results.advection .* reference_norm ./ max(system.advection_base_norms(i_t), eps);
    end
end

rhs = system.data_rhs + prior_2d .* prior_scales .^ 2;
diagonal = system.data_diagonal + prior_scales .^ 2;
if any(smoothness_scales > 0)
    diagonal = diagonal + system.laplacian_normal_diagonal .* smoothness_scales .^ 2;
end
for i_t = 1:n_time - 1
    if advection_scales(i_t) == 0
        continue
    end
    step = system.steps{i_t};
    coeff_squared = advection_scales(i_t) ^ 2 .* step.coeff_valid .^ 2;
    if isscalar(coeff_squared)
        coeff_squared = repmat(coeff_squared, numel(step.source_id_valid), 1);
    end
    diagonal(step.source_id_valid, i_t) = diagonal(step.source_id_valid, i_t) + coeff_squared;
    diagonal(:, i_t + 1) = diagonal(:, i_t + 1) + accumarray(double(step.dest_id_valid), coeff_squared, [nxyz 1]);
end
diagonal(diagonal == 0) = 1;

apply_normal = @(x) normal_operator(x, system, prior_scales, smoothness_scales, advection_scales);
[x, flag, relres, iter] = pcg(apply_normal, rhs(:), p.Results.tol, p.Results.maxit, @(x) x ./ diagonal(:), [], prior_2d(:));

fields = reshape(x, [system.grid_sz n_time]);
info = struct("controls", p.Results, ...
    "solver", struct("name", "batched_pcg_normal", "flag", flag, "relres", relres, "iterations", iter, ...
    "n_time", n_time, "n_unknown_total", numel(x)));
end

function y = normal_operator(x, system, prior_scales, smoothness_scales, advection_scales)
X = reshape(x, system.nxyz, system.n_time);
projection = system.F * X;
projection(~system.obs_mask) = 0;
Y = (system.F' * projection) .* system.data_scales .^ 2 + X .* prior_scales .^ 2;
if any(smoothness_scales > 0)
    Y = Y + (system.laplacian' * (system.laplacian * X)) .* smoothness_scales .^ 2;
end
for i_t = 1:system.n_time - 1
    if advection_scales(i_t) == 0
        continue
    end
    step = system.steps{i_t};
    residual = step.coeff_valid .* (X(step.dest_id_valid, i_t + 1) - X(step.source_id_valid, i_t));
    weighted = advection_scales(i_t) ^ 2 .* step.coeff_valid .* residual;
    Y(step.source_id_valid, i_t) = Y(step.source_id_valid, i_t) - weighted;
    Y(:, i_t + 1) = Y(:, i_t + 1) + accumarray(double(step.dest_id_valid), weighted, [system.nxyz 1]);
end
y = Y(:);
end
