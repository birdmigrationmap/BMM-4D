function [vp, rho, eta, lambda, std_vp] = compute_vp(obs_vector, g, F, varargin)
%COMPUTE_VP Estimate a vertical profile with Tikhonov regularization.
%
%   [vp, rho, eta, lambda] = compute_vp(obs_vector, g, F)
%   estimates the vertical profile vp with derivative-penalized Tikhonov
%   regularization.
%
%   [vp, rho, eta, lambda, std_vp] = compute_vp(...)
%   also estimates the standard deviation of each vertical bin from the
%   regularized inverse covariance approximation.
%
%   Optional name-value arguments:
%     'order_derivative_smoothness' : derivative order used in L, default 2.
%     'lambda'                      : regularization parameter, default NaN.
%                                     If NaN, lambda is selected by the
%                                     L-curve.
%     'vp_0'                        : prior/background profile, default zeros.

p = inputParser;
addParameter(p, 'order_derivative_smoothness', 1);
addParameter(p, 'lambda', nan);
addParameter(p, 'vp_0', zeros(g.sz(3), 1));
parse(p, varargin{:});

% Keep only non-zero observations and matching rows of F.
[obs_values, nonzero_obs_idx] = bmm4d.forward.format_data(obs_vector(:, 1));

% Map each 3D g cell to its vertical bin.
vertical_aggregation = repelem( ...
    speye(g.sz(3)), ...
    prod(g.sz(1:2)), ...
    1 ...
);

% Observation model:
%   b = A * vp
% where A maps a vertical profile to the observed values.
A_system = F(nonzero_obs_idx, :) * vertical_aggregation;
b_system = obs_values;

% Smoothness penalty matrix L.
% The solution minimizes:
%   ||A * vp - b||^2 + lambda^2 * ||L * (vp - vp_0)||^2
L = bmm4d.util.difference_operator(g.sz(3), p.Results.order_derivative_smoothness);

% Generalized SVD used by RegTools Tikhonov routines.
[U, sm, X, ~, ~] = cgsvd(A_system, L);

% Select regularization strength.
if isnan(p.Results.lambda)
    % L-curve.
    % l_curve(U, sm, b_system); % to plot
    regparam = bmm4d.prior.l_curve_wo_fig(U, sm, b_system);
    lambda = regparam(1);

else
    lambda = p.Results.lambda;
end

% Regularized vertical profile.
% rho = ||A * vp - b||_2
% eta = ||L * (vp - vp_0)||_2
[vp, rho, eta] = tikhonov(U, sm, X, b_system, lambda, p.Results.vp_0);

% Optional uncertainty estimate.
% Only computed when the caller requests the fifth output.
if nargout >= 5
    sigma = sm(:, 1);
    mu = sm(:, 2);
    sigma = sigma(:);
    mu = mu(:);

    residuals = b_system - A_system * vp;
    residuals = residuals(:);

    % Effective degrees of freedom for Tikhonov regularization.
    filter_factors = sigma.^2 ./ (sigma.^2 + lambda^2 * mu.^2);
    dof = sum(filter_factors, "all");

    % Residual-based observation noise variance estimate.
    sigma2 = sum(residuals.^2, "all") / max(numel(b_system) - dof, eps);

    % Diagonal weights of the regularized inverse covariance.
    w = sigma.^2 ./ (sigma.^2 + lambda^2 * mu.^2).^2;

    % Cov(vp) ≈ sigma2 * X * diag(w) * X'
    k = length(w);
    cov_vp = sigma2 * X(:, 1:k) * diag(w) * X(:, 1:k)';

    std_vp = sqrt(diag(cov_vp));
end

end
