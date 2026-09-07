function [rmse_grid, metrics] = compute_scan_calibration_grid(F, obs_vector, g, target_scan, varargin)
%COMPUTE_SCAN_CALIBRATION_GRID Evaluate COMPUTE_SCAN-style calibration weights.
%   RMSE_GRID = COMPUTE_SCAN_CALIBRATION_GRID(F, OBS_VECTOR, GRID, TARGET_SCAN, ...)
%   returns the RMSE of each parameter combination, shaped like the
%   NDGRID of vp_absolute_values, smoothness_values, and coeff_horver_values.
%   The candidate weights are passed to COMPUTE_SCAN, where each residual
%   block is normalized by its row count before applying the weight.
%
%   [RMSE_GRID, METRICS] also returns normalized RMSE, bias, and pattern
%   correlation grids using the same score mask.
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

x_grid = g.x(:) - g.x(end) / 2;
y_grid = g.y(:) - g.y(end) / 2;
[X_grid, Y_grid] = ndgrid(x_grid, y_grid);
score_mask = repmat(hypot(X_grid, Y_grid) <= p.Results.r_lim * 1000, 1, 1, g.sz(3));

rmse_grid = nan(size(vp_absolute_grid));
nrmse_grid = nan(size(vp_absolute_grid));
bias_grid = nan(size(vp_absolute_grid));
correlation_grid = nan(size(vp_absolute_grid));

num_grid_cells = prod(g.sz);
num_horizontal_cells = prod(g.sz(1:2));
[obs_values, nonzero_obs_idx] = format_data(obs_vector(:, 1));
data_scale = 1 / sqrt(max(numel(obs_values), 1));
A_data = data_scale .* F(nonzero_obs_idx, :);
b_data = data_scale .* obs_values;
data_block_norm = max(norm(A_data, "fro"), 1);

identity_operator = speye(num_grid_cells);
identity_norm = sqrt(num_grid_cells);

laplacian_operators = cell(numel(p.Results.coeff_horver_values), 1);
laplacian_norms = zeros(numel(p.Results.coeff_horver_values), 1);
if any(p.Results.smoothness_values > 0)
    for i_coeff = 1:numel(p.Results.coeff_horver_values)
        laplacian_operators{i_coeff} = LaplacianMatrix(g, water=p.Results.water, coeff_horver=p.Results.coeff_horver_values(i_coeff));
        laplacian_norms(i_coeff) = max(norm(laplacian_operators{i_coeff}, "fro"), eps);
    end
end

use_vertical_profile = any(p.Results.vp_absolute_values > 0);
if use_vertical_profile
    vp = compute_vp(obs_vector(:, 1), g, F);
    vertical_profile_expansion = repelem(speye(g.sz(3)), num_horizontal_cells, 1);
    if p.Results.water
        if isfield(g, "mask_water")
            water_mask = g.mask_water;
        else
            water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
            water_mask = logical(water_mask_data.Fonwater({g.f_lat(g.y), g.f_lon(g.x)}))';
        end
        water_mask_3d = repmat(water_mask, 1, 1, g.sz(3));
        vertical_profile_expansion(water_mask_3d(:) == 1, :) = 0;
    end
    vertical_profile_field = vertical_profile_expansion * vp;
else
    vertical_profile_field = sparse(num_grid_cells, 1);
end

target_scan = double(target_scan);
target_values = target_scan(score_mask);
target_scale = sqrt(mean(target_values .^ 2));
target_mean = mean(target_values);
scan_initial = zeros(num_grid_cells, 1);

for i_combination = 1:numel(vp_absolute_grid)
    if p.Results.verbose
        fprintf("Calibration %d / %d: vp_absolute = %.4g, smoothness = %.4g, coeff_horver = %.4g\n", ...
            i_combination, ...
            numel(vp_absolute_grid), ...
            vp_absolute_grid(i_combination), ...
            smoothness_grid(i_combination), ...
            coeff_horver_grid(i_combination));
    end

    A_parts = cell(3, 1);
    b_parts = cell(3, 1);
    n_part = 1;
    A_parts{n_part} = A_data;
    b_parts{n_part} = b_data;

    if smoothness_grid(i_combination) > 0
        i_coeff = find(p.Results.coeff_horver_values == coeff_horver_grid(i_combination), 1);
        smoothness_scale = smoothness_grid(i_combination) .* data_block_norm ./ laplacian_norms(i_coeff);
        n_part = n_part + 1;
        A_parts{n_part} = smoothness_scale .* laplacian_operators{i_coeff};
        b_parts{n_part} = sparse(num_grid_cells, 1);
    end

    if vp_absolute_grid(i_combination) > 0
        vp_absolute_scale = vp_absolute_grid(i_combination) .* data_block_norm ./ identity_norm;
        n_part = n_part + 1;
        A_parts{n_part} = vp_absolute_scale .* identity_operator;
        b_parts{n_part} = vp_absolute_scale .* vertical_profile_field;
    end

    A_system = vertcat(A_parts{1:n_part});
    b_system = vertcat(b_parts{1:n_part});
    scan_1d = lsqr(A_system, b_system, p.Results.tol, p.Results.maxit, [], [], scan_initial);
    scan_initial = scan_1d;
    scan = reshape(scan_1d, g.sz);

    err = scan - target_scan;
    err_values = err(score_mask);
    scan_values = scan(score_mask);
    rmse_grid(i_combination) = sqrt(mean(err_values .^ 2));
    nrmse_grid(i_combination) = rmse_grid(i_combination) ./ target_scale;
    bias_grid(i_combination) = mean(err_values) ./ max(abs(target_mean), eps);
    correlation_grid(i_combination) = corr(scan_values(:), target_values(:), Rows="complete");

end
correlation_grid(~isfinite(correlation_grid)) = nan;

metrics = struct( ...
    "rmse", rmse_grid, ...
    "nrmse", nrmse_grid, ...
    "bias", bias_grid, ...
    "correlation", correlation_grid, ...
    "vp_absolute_values", p.Results.vp_absolute_values, ...
    "smoothness_values", p.Results.smoothness_values, ...
    "coeff_horver_values", p.Results.coeff_horver_values, ...
    "r_lim", p.Results.r_lim);
end
