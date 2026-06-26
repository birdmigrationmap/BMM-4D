function Dsim = simulate_gaussian_truth(grid, cov_parm, trans_finv, vp, n_sim, varargin)
%SIMULATE_GAUSSIAN_TRUTH Generate synthetic density fields in [x y z] order.
% Gaussian fields are transformed back to density units using the fitted
% altitude-dependent inverse transform and vertical profile, then clipped
% at zero. No permutation is applied: outputs follow the project [x y z]
% grid convention and can be forward-projected with F * Dsim{i}(:).

p = inputParser;
addParameter(p, "seed", "shuffle");
parse(p, varargin{:});

if exist("FGS", "file") ~= 2
    error("simulate_gaussian_truth:FGSNotFound", "FastGaussianSimulation is not on the MATLAB path. Add it in setup_paths.m or external/FastGaussianSimulation.")
end

cov_parm = cov_parm(:);

sim.s = grid.sz;
sim.n = n_sim;
sim.seed = p.Results.seed;
sim.DisplayProgression = false;

covar(1).model = 'k-bessel';
covar(1).azimuth = [0 0];
covar(1).var = cov_parm(1);
covar(1).alpha = cov_parm(4);
covar(1).range = cov_parm([3 3 2]) ./ [grid.dx; grid.dy; grid.dz];

z_grid = repmat(reshape(squeeze(grid.z), 1, 1, []), grid.sz(1), grid.sz(2), 1);
vp_grid = repmat(reshape(vp(:), 1, 1, []), grid.sz(1), grid.sz(2), 1);

Zsim = FGS(sim, covar);
Dsim = cell(n_sim, 1);
for i_sim = 1:n_sim
    Dsim{i_sim} = trans_finv(Zsim{i_sim}, z_grid) + vp_grid;
    Dsim{i_sim}(Dsim{i_sim} < 0) = 0;
end

if false
    % Inspect orthogonal slices through one simulated density field.
    x_plot = grid.x - grid.x(end)/2;
    y_plot = grid.y - grid.y(end)/2;
    z_plot = squeeze(grid.z);
    figure;
    h = slice(x_plot/1000, y_plot/1000, z_plot, permute(Dsim{1}, [2 1 3]), 0, 0, 200);
    set(h, EdgeColor="none")
    xlabel("x (km)"); ylabel("y (km)"); zlabel("Altitude (m)"); axis tight
end
end
