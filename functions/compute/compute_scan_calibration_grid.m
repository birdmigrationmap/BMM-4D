function rmse_grid = compute_scan_calibration_grid(F, obs_vector, grid, target_scan, varargin)
%COMPUTE_SCAN_CALIBRATION_GRID Evaluate COMPUTE_SCAN on a parameter grid.
%   RMSE_GRID = COMPUTE_SCAN_CALIBRATION_GRID(F, OBS_VECTOR, GRID, TARGET_SCAN, ...)
%   returns the RMSE of each parameter combination, shaped like the
%   NDGRID of vp_absolute_values, smoothness_values, and coeff_horver_values.
%   The candidate weights are passed to COMPUTE_SCAN, where each residual
%   block is normalized by its row count before applying the weight.
%
%   Optional name-value arguments:
%     'vp_absolute_values' : vector of RMS-normalized candidate weights
%     'smoothness_values'  : vector of RMS-normalized candidate weights
%     'coeff_horver_values': vector of candidate anisotropy values
%     'water'              : passed to COMPUTE_SCAN, default false
%     'tol'                : LSQR tolerance, default 1e-6
%     'maxit'              : LSQR iterations, default 200
%     'r_lim'              : max horizontal radius in km for scoring, default 50
%     'verbose'            : print progress, default true

p = inputParser;
addParameter(p, 'vp_absolute_values', 1);
addParameter(p, 'smoothness_values', 10);
addParameter(p, 'coeff_horver_values', 0);
addParameter(p, 'water', false);
addParameter(p, 'tol', 1e-6);
addParameter(p, 'maxit', 200);
addParameter(p, 'r_lim', 50);
addParameter(p, 'verbose', true);
parse(p, varargin{:});

[vp_absolute_grid, smoothness_grid, coeff_horver_grid] = ndgrid( ...
    p.Results.vp_absolute_values, ...
    p.Results.smoothness_values, ...
    p.Results.coeff_horver_values ...
);

x_grid = grid.x(:) - grid.x(end) / 2;
y_grid = grid.y(:) - grid.y(end) / 2;
[X_grid, Y_grid] = ndgrid(x_grid, y_grid);
score_mask = repmat(hypot(X_grid, Y_grid) <= p.Results.r_lim * 1000, 1, 1, grid.sz(3));

rmse_grid = nan(size(vp_absolute_grid));

for i_combination = 1:numel(vp_absolute_grid)
    if p.Results.verbose
        fprintf("Calibration %d / %d: vp_absolute = %.4g, smoothness = %.4g, coeff_horver = %.4g\n", ...
            i_combination, ...
            numel(vp_absolute_grid), ...
            vp_absolute_grid(i_combination), ...
            smoothness_grid(i_combination), ...
            coeff_horver_grid(i_combination));
    end

    % COMPUTE_SCAN owns VP estimation and Laplacian construction so every
    % calibration point follows the same production inversion workflow.
    scan = compute_scan(F, obs_vector, grid, ...
        "vp_absolute", vp_absolute_grid(i_combination), ...
        "smoothness", smoothness_grid(i_combination), ...
        "coeff_horver", coeff_horver_grid(i_combination), ...
        "water", p.Results.water, ...
        "tol", p.Results.tol, ...
        "maxit", p.Results.maxit);

    err = scan - target_scan;
    err = err(score_mask);
    rmse_grid(i_combination) = sqrt(mean(err.^2));

end
end
