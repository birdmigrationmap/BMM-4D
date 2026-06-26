function [vp, rho, eta, lambda, std_vp] = compute_vp_thiko(obs_vector, grid, F, varargin)
%COMPUTE_VP_THIKO Estimate a vertical profile with standard-form Tikhonov.
%
%   [vp, rho, eta, lambda] = compute_vp_thiko(obs_vector, grid, F)
%   estimates the vertical profile vp using standard-form Tikhonov
%   regularization.
%
%   [vp, rho, eta, lambda, std_vp] = compute_vp_thiko(...)
%   also estimates the standard deviation of each vertical bin from the
%   regularized inverse covariance approximation.
%
%   Optional name-value arguments:
%     'lambda' : regularization parameter, default NaN
%                If NaN, lambda is selected by the L-curve.
%     'vp_0'   : prior/background vertical profile, default zeros.

p = inputParser;
addParameter(p, 'lambda', nan);
addParameter(p, 'vp_0', zeros(grid.sz(3), 1));
parse(p, varargin{:});

% Keep only non-zero observations and matching rows of F.
[obs_values, nonzero_obs_idx] = format_data(obs_vector(:, 1));

% Map each 3D grid cell to its vertical bin.
vertical_aggregation = repelem( ...
    speye(grid.sz(3)), ...
    prod(grid.sz(1:2)), ...
    1 ...
);

% Observation model:
%   b = A * vp
% where A maps a vertical profile to the observed values.
A_system = F(nonzero_obs_idx, :) * vertical_aggregation;
b_system = obs_values;

% Compact SVD used by RegTools standard-form Tikhonov routines.
[U, S, V] = csvd(A_system);

% Select regularization strength.
if isnan(p.Results.lambda)
    % L-curve.
    % l_curve(U, S, b_system); % to plot
    lambda = l_curve_wo_fig(U, S, b_system, 'Tikh');

    % Alternative criteria:
    % lambda = gcv_wo_fig(U, S, b_system);
    % delta = sigma * sqrt(length(b_system));
    % lambda = discrep(U, S, b_system, delta);
    % lambda = quasiopt(U, S, b_system);
else
    lambda = p.Results.lambda;
end

% Regularized vertical profile.
% rho = ||A * vp - b||_2
% eta = ||vp - vp_0||_2
[vp, rho, eta] = tikhonov(U, S, V, b_system, lambda, p.Results.vp_0);

% Optional uncertainty estimate.
% Only computed when the caller requests the fifth output.
if nargout >= 5
    singular_values = S(:);

    residuals = b_system - A_system * vp;
    residuals = residuals(:);

    % Effective degrees of freedom for standard-form Tikhonov.
    filter_factors = singular_values.^2 ./ ...
        (singular_values.^2 + lambda^2);

    dof = sum(filter_factors, "all");

    % Residual-based observation noise variance estimate.
    sigma2 = sum(residuals.^2, "all") / max(numel(b_system) - dof, eps);

    % Diagonal weights of the regularized inverse covariance.
    w = singular_values.^2 ./ ...
        (singular_values.^2 + lambda^2).^2;

    % csvd returns singular values as a vector and V with matching columns.
    % Cov(vp) ≈ sigma2 * V * diag(w) * V'
    k = length(w);
    cov_vp = sigma2 * V(:, 1:k) * diag(w) * V(:, 1:k)';

    std_vp = sqrt(diag(cov_vp));
end

end
