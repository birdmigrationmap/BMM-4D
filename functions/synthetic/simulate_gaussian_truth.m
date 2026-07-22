function Dsim = simulate_gaussian_truth(g, cov_parm, trans_finv, vp, n_sim, varargin)
%SIMULATE_GAUSSIAN_TRUTH Generate synthetic density fields in [x y z] order.
% Gaussian fields are transformed back to density units using the fitted
% altitude-dependent inverse transform and vertical profile, then clipped
% at zero. No permutation is applied: outputs follow the project [x y z]
% g convention and can be forward-projected with F * Dsim{i}(:).

p = inputParser;
addParameter(p, "seed", "shuffle");
parse(p, varargin{:});

if exist("FGS", "file") ~= 2
    error("simulate_gaussian_truth:FGSNotFound", "FastGaussianSimulation is not on the MATLAB path. Add it in setup_paths.m or external/FastGaussianSimulation.")
end

cov_parm = cov_parm(:);

sim.s = g.sz;
sim.n = n_sim;
sim.seed = p.Results.seed;
sim.DisplayProgression = false;

covar(1).model = 'k-bessel';
covar(1).azimuth = [0 0];
covar(1).var = cov_parm(1);
covar(1).alpha = cov_parm(4);
covar(1).range = cov_parm([3 3 2]) ./ [g.dx; g.dy; g.dz];

z_grid = repmat(reshape(squeeze(g.z), 1, 1, []), g.sz(1), g.sz(2), 1);
vp_grid = repmat(reshape(vp(:), 1, 1, []), g.sz(1), g.sz(2), 1);

Zsim = FGS(sim, covar);
Dsim = cell(n_sim, 1);
for i_sim = 1:n_sim
    Dsim{i_sim} = trans_finv(Zsim{i_sim}, z_grid) + vp_grid;
    Dsim{i_sim}(Dsim{i_sim} < 0) = 0;
end
end
