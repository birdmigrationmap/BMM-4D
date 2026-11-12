function result = ar_solve(S, name, varargin)
%AR_SOLVE Solve a test case with the production inversion settings (dev).
%   RESULT = AR_SOLVE(S, NAME, ...) solves S (ar_load) with
%   bmm4d.solve.solve_advection_lean and the settings of c.inversion
%   (production weights, advection relaxation), as bmm4d.invert, except:
%   - "prior": the prior mean [n_active n_time] (default S.prior);
%   - "steps": the advection steps (default S.steps; steps with link
%     ratios, ar_steps "sigma", are solved by ar_solve_advection_lean);
%   - "nonnegative": default false (plain PCG, faster; production true);
%   - "x_initial": warm start;
%   - "prior_correlation_m": correlated prior (default 0, white);
%   - "prior_correlation_vertical_m": its vertical correlation (default 0);
%   - "prior_matern_m": [l_h l_v] (m), one Matern prior term (option B of
%     PRIOR_OPTIONS.md; dev solver; pass smoothness 0);
%   - "smoothness": Laplacian smoothness weight (default c.inversion.smoothness);
%   - "advection_reference": field whose ratios relax the advection links
%     (takeoff, landing; default the prior mean used).
%   RESULT has the name, the field, the column maps [nx ny n_time] (lowest
%   cell by its air share, NaN outside active columns), the relative data
%   misfit on all rows and on rows 30-120 km from their radar, the solve
%   time and iterations.

p = inputParser;
addParameter(p, "prior", S.prior);
addParameter(p, "steps", S.steps);
addParameter(p, "nonnegative", false);
addParameter(p, "x_initial", []);
addParameter(p, "prior_correlation_m", 0);
addParameter(p, "prior_correlation_vertical_m", 0);
addParameter(p, "prior_matern_m", []);
addParameter(p, "smoothness", S.c.inversion.smoothness);
addParameter(p, "advection_reference", []);
parse(p, varargin{:});
opt = p.Results;
if isempty(opt.advection_reference)
    opt.advection_reference = opt.prior;
end
inv = S.c.inversion;

% Links with density ratios (ar_steps "sigma") need the dev solver
% and so does the one-term Matern prior
solver = @bmm4d.solve.solve_advection_lean;
extra = {};
if isfield(opt.steps{1}, "ratio_valid") || ~isempty(opt.prior_matern_m)
    solver = @ar_solve_advection_lean;
    extra = {"prior_matern_m", opt.prior_matern_m};
end
solve_tic = tic;
[field, info] = solver(S.observations, S.g, S.active_mask_3d, opt.steps, opt.prior, ...
    "prior", inv.prior, "smoothness", opt.smoothness, "advection", inv.advection, "coeff_horver", inv.coeff_horver, ...
    "tol", inv.tol, "maxit", inv.maxit, "store_class", "single", "advection_relax_exponent", inv.advection_relax_exponent, ...
    "advection_reference", opt.advection_reference, "nonnegative", opt.nonnegative, "projected_maxit", inv.projected_maxit, ...
    "row_weights", S.row_weights, "x_initial", opt.x_initial, "prior_correlation_m", opt.prior_correlation_m, ...
    "prior_correlation_vertical_m", opt.prior_correlation_vertical_m, extra{:});
seconds_used = toc(solve_tic);

% Misfit on all observed rows and on rows 30-120 km from their radar
predicted = S.observations.F(:, S.active_id) * double(field);
observed = full(S.observations.values);
mask = observed ~= 0;
misfit = norm(predicted(mask) - observed(mask)) / norm(observed(mask));
trusted = mask & S.observations.row_range_m >= 30000 & S.observations.row_range_m <= 120000;
misfit_trusted = norm(predicted(trusted) - observed(trusted)) / norm(observed(trusted));

result = struct("name", string(name), "field", field, "column", ar_columns(S, field), ...
    "misfit", misfit, "misfit_trusted", misfit_trusted, "seconds", seconds_used, ...
    "iterations", info.solver.iterations);
fprintf("%s: %.0f s, %d iterations, misfit %.3f (30-120 km %.3f)\n", name, seconds_used, ...
    info.solver.iterations, misfit, misfit_trusted)
end
