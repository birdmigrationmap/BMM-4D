function T = ar_prior_pcg_toy()
%AR_PRIOR_PCG_TOY PCG iterations of the normal equations for each prior form (dev toy).
% Used by PRIOR_OPTIONS.md (section 5). Seconds; no data.
% 2D lattice 160 x 160 at 4 km, four radars 200 km apart observing 50% of
% the cells within 120 km through a 3 x 3 beam-like average; weights scaled
% as in solve_advection_lean (data by observation count, prior by sqrt(n),
% smoothness by the Laplacian's Frobenius norm). Jacobi PCG to 1e-4 (the
% production tolerance), or PCG with the prior-aware preconditioner
% M = prior part + diag(data), factorized exactly (best case).
rng(1)
nx = 160; n = nx ^ 2; dx_km = 4;
e = ones(nx, 1);
D1 = spdiags([-e 2 * e -e], -1:1, nx, nx);
D1(1, 1) = 1; D1(end, end) = 1;                 % graph Laplacian (Neumann)
Lg = kron(speye(nx), D1) + kron(D1, speye(nx)); % unit 2D graph Laplacian
Ls = Lg / 4;                                    % smoothness operator, unit diagonal

% Observations: 3 x 3 average around observed cells within 120 km of a radar
[ix, iy] = ndgrid(1:nx);
radar = [55 55; 105 55; 55 105; 105 105];
covered = false(nx);
for r = 1:4
    covered = covered | hypot(ix - radar(r, 1), iy - radar(r, 2)) <= 30;
end
obs_cell = find(covered & rand(nx) < 0.5);
[ox, oy] = ind2sub([nx nx], obs_cell);
[sx, sy] = ndgrid(-1:1);
rows = repelem((1:numel(obs_cell))', 9);
cols = sub2ind([nx nx], min(max(repelem(ox, 9) + repmat(sx(:), numel(ox), 1), 1), nx), ...
    min(max(repelem(oy, 9) + repmat(sy(:), numel(oy), 1), 1), nx));
F = sparse(rows, cols, 1 / 9, numel(obs_cell), n);
truth = smooth_field(nx, 8) + 0.5 * smooth_field(nx, 40);
y = F * truth(:);
data_scale2 = 1 / numel(obs_cell);
data_part = data_scale2 * (F' * F);
rhs_data = data_scale2 * (F' * y);
lp = 0.3; ls = 0.3;
prior_scale2 = lp ^ 2 / n;
smooth_scale2 = ls ^ 2 / norm(Ls, "fro") ^ 2;

% Marginal-variance normalization on a periodic lattice, as in the solver
w = 2 * pi * (0:255) / 256;
[kx, ky] = ndgrid(w, w);
lam = 4 - 2 * cos(kx) - 2 * cos(ky);
cases = {
    "white + smoothness 0.3 (v2 production)", speye(n), true
    "GP alpha=1 100 km + smoothness 0.3 (dev)", gp1(Lg, 100 / dx_km, lam), true
    "SPDE alpha=2 20 km", gp2(Lg, 20 / dx_km, lam), false
    "SPDE alpha=2 50 km", gp2(Lg, 50 / dx_km, lam), false
    };
out = cell(0, 1);
for i = 1:size(cases, 1)
    [name, Qp, use_smooth] = cases{i, :};
    A = data_part + prior_scale2 * Qp + use_smooth * smooth_scale2 * (Ls' * Ls);
    b = rhs_data;                       % prior mean 0
    jac = full(diag(A));
    tic
    [~, ~, ~, it_j] = pcg(A, b, 1e-4, 5000, @(v) v ./ jac);
    t_j = toc;
    % Prior-aware preconditioner: the prior (and smoothness) part kept whole
    M = prior_scale2 * Qp + use_smooth * smooth_scale2 * (Ls' * Ls) + spdiags(full(diag(data_part)), 0, n, n);
    R = chol(M);
    tic
    [~, ~, ~, it_p] = pcg(A, b, 1e-4, 5000, R', R);
    t_p = toc;
    out{end + 1, 1} = table(name, it_j, it_p, t_j, t_p, ...
        VariableNames=["prior" "iterations_jacobi" "iterations_prior_aware" "s_jacobi" "s_prior_aware"]); %#ok<AGROW>
end
T = vertcat(out{:});
end

function Q = gp1(Lg, l, lam)
Q = speye(size(Lg, 1)) + l ^ 2 * Lg;
Q = mean(1 ./ (1 + l ^ 2 * lam), "all") * Q;
end

function Q = gp2(Lg, l, lam)
K = speye(size(Lg, 1)) + l ^ 2 * Lg;
Q = mean(1 ./ (1 + l ^ 2 * lam) .^ 2, "all") * (K' * K);
end

function f = smooth_field(nx, l)
% Unit-variance smooth random field (Gaussian filter of white noise, periodic)
[kx, ky] = ndgrid([0:nx / 2 - 1, -nx / 2:-1] * 2 * pi / nx);
f = real(ifft2(fft2(randn(nx)) .* exp(-(kx .^ 2 + ky .^ 2) * l ^ 2 / 4)));
f = f / std(f(:));
end
