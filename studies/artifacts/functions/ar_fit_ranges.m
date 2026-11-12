function fit = ar_fit_ranges(observations, g, prior_active, active_id, n_pairs_per_bin)
%AR_FIT_RANGES Fit the prior's correlation ranges on the data, beam-aware (dev).
%   FIT = AR_FIT_RANGES(OBSERVATIONS, G, PRIOR_ACTIVE, ACTIVE_ID,
%   N_PAIRS_PER_BIN) estimates the horizontal and vertical ranges of the
%   anomaly from the prior mean (the profile field PRIOR_ACTIVE, [n_active
%   n_time]) for two prior forms used by the inversion:
%   - "gp": precision I + (l_h/dx)^2 L_h + (l_v/dz)^2 L_v (the correlated
%     prior of fix 3, with a vertical part);
%   - "matern": its square (one Matern term, PRIOR_OPTIONS.md option B).
%   Method: the residual r_i = y_i - (F mu)_i of observations at 20-120 km
%   from their radar; for pairs (i, j) of the same frame, sampled with at
%   most N_PAIRS_PER_BIN (default 150) per frame and distance bin, the
%   product r_i r_j is compared with the model covariance of the two
%   beam-averaged observations, s2 * F_i R F_j' (R the model's correlation
%   on the solve lattice, computed by FFT of 1 / Q), so the beam width is
%   accounted for. Pairs are binned by horizontal and vertical distance of
%   their beam centres; s2 and the ranges are fitted by weighted least
%   squares on the bin means (pairs i = j excluded: the nugget, i.e. the
%   observation error, is left free).
%   FIT has, per form, the ranges l_h (m), l_v (m), s2 and the residual,
%   and the binned empirical covariance for plotting.

if nargin < 5
    n_pairs_per_bin = 150;
end
rng(7)
n_time = size(observations.values, 2);
F = observations.F(:, active_id);
[ix, iy, iz] = ind2sub(g.sz, active_id);
cell_xyz = [double(ix) double(iy) double(iz)];
row_mass = full(sum(F, 2));
centre = (F * cell_xyz) ./ max(row_mass, eps);             % beam centres (cells)
trusted = observations.row_range_m >= 20000 & observations.row_range_m <= 120000 & row_mass > 0;

% Bins by beam-centre distance: horizontal (km) and vertical (m)
h_edges = [0 2 4 6 8 12 16 20 25 30 40 50 60 80 100 125 150 200 250 300];
v_edges = [0 100 200 300 500 800 1200 2000 3000];
n_bin = (numel(h_edges) - 1) * (numel(v_edges) - 1);

% Sample pairs per frame: many random candidates (long lags) and partners
% among the nearest rows (short lags), then at most N_PAIRS_PER_BIN per
% distance bin, so that every lag range is sampled alike
pair_i = []; pair_j = []; pair_product = []; bin = [];
for i_t = 1:n_time
    rows = find(trusted & observations.values(:, i_t) ~= 0);
    if numel(rows) < 1000
        continue
    end
    residual = full(observations.values(rows, i_t)) - F(rows, :) * double(prior_active(:, i_t));
    residual = residual - mean(residual);
    xyz = centre(rows, :);
    a = randi(numel(rows), 2e6, 1);
    b = randi(numel(rows), 2e6, 1);
    anchor = randi(numel(rows), 1e4, 1);
    near = knnsearch(xyz(:, 1:2), xyz(anchor, 1:2), K=min(300, numel(rows)));
    a = [a; repmat(anchor, size(near, 2), 1)]; %#ok<AGROW>
    b = [b; near(:)]; %#ok<AGROW>
    h_km = hypot(xyz(a, 1) - xyz(b, 1), xyz(a, 2) - xyz(b, 2)) * g.dx / 1000;
    v_m = abs(xyz(a, 3) - xyz(b, 3)) * g.dz;
    [~, ~, h_bin] = histcounts(h_km, h_edges);
    [~, ~, v_bin] = histcounts(v_m, v_edges);
    valid = find(h_bin > 0 & v_bin > 0 & a ~= b);
    pair_bin = sub2ind([numel(h_edges) - 1, numel(v_edges) - 1], h_bin(valid), v_bin(valid));
    % Random order, then the first N_PAIRS_PER_BIN of each bin
    order = randperm(numel(valid))';
    valid = valid(order);
    pair_bin = pair_bin(order);
    [pair_bin, sorted] = sort(pair_bin);
    valid = valid(sorted);
    first = [true; diff(pair_bin) > 0];
    group_start = cummax((1:numel(pair_bin))' .* first);
    rank_in_bin = (1:numel(pair_bin))' - group_start + 1;
    take = rank_in_bin <= n_pairs_per_bin;
    pick = valid(take);
    pair_i = [pair_i; rows(a(pick))]; %#ok<AGROW>
    pair_j = [pair_j; rows(b(pick))]; %#ok<AGROW>
    pair_product = [pair_product; residual(a(pick)) .* residual(b(pick))]; %#ok<AGROW>
    bin = [bin; pair_bin(take)]; %#ok<AGROW>
end
n_pair = numel(pair_i);
bin_count = accumarray(bin, 1, [n_bin 1]);
bin_cov = accumarray(bin, pair_product, [n_bin 1]) ./ max(bin_count, 1);

% Lag list of every pair: products of the two rows' cell weights
Ft = F';
lag = cell(numel(pair_i), 1);
weight = cell(numel(pair_i), 1);
for k = 1:numel(pair_i)
    [ci, ~, wi] = find(Ft(:, pair_i(k)));
    [cj, ~, wj] = find(Ft(:, pair_j(k)));
    wi = wi / sum(wi); wj = wj / sum(wj);
    [I, J] = ndgrid(1:numel(ci), 1:numel(cj));
    lag{k} = int16(abs(cell_xyz(ci(I(:)), :) - cell_xyz(cj(J(:)), :)));
    weight{k} = single(wi(I(:)) .* wj(J(:)));
end
pair_of = repelem((1:numel(pair_i))', cellfun(@numel, weight));
lag = vertcat(lag{:});
weight = vertcat(weight{:});
pair_bin = bin;

% Model correlation on the solve lattice (periodic, large enough)
n_h = 2 ^ nextpow2(2 * max(double(lag(:, 1:2)), [], "all") + 64);
n_v = 2 ^ nextpow2(2 * g.sz(3));
wh = 2 * pi * (0:n_h - 1) / n_h;
wv = 2 * pi * (0:n_v - 1) / n_v;
[kx, ky, kz] = ndgrid(wh, wh, wv);
lambda_h = 4 - 2 * cos(kx) - 2 * cos(ky);
lambda_v = 2 - 2 * cos(kz);
clear kx ky kz
lag_index = sub2ind([n_h n_h n_v], double(lag(:, 1)) + 1, double(lag(:, 2)) + 1, double(lag(:, 3)) + 1);
clear lag

forms = ["gp" "matern"];
powers = [1 2];
fit = struct("n_pairs", n_pair, "h_edges", h_edges, "v_edges", v_edges, ...
    "bin_count", reshape(bin_count, numel(h_edges) - 1, []), "bin_cov", reshape(bin_cov, numel(h_edges) - 1, []));
use = bin_count >= 30;
for i_f = 1:2
    objective = @(log_l) model_misfit(exp(log_l), powers(i_f));
    log_l = fminsearch(objective, log([50000 300]), optimset("TolX", 1e-3, "MaxFunEvals", 200));
    [misfit, s2, model] = model_misfit(exp(log_l), powers(i_f));
    fit.(forms(i_f)) = struct("l_h_m", exp(log_l(1)), "l_v_m", exp(log_l(2)), "s2", s2, "misfit", misfit, ...
        "bin_model", reshape(model, numel(h_edges) - 1, []));
    fprintf("%s: l_h %.0f km, l_v %.0f m, s2 %.3g, weighted misfit %.3g\n", forms(i_f), ...
        exp(log_l(1)) / 1000, exp(log_l(2)), s2, misfit)
end

    function [misfit, s2, model] = model_misfit(l, power)
        % Correlation of the beam-averaged pairs under ranges l = [l_h l_v]
        q = (1 + (l(1) / g.dx) ^ 2 * lambda_h + (l(2) / g.dz) ^ 2 * lambda_v) .^ power;
        correlation = real(ifftn(1 ./ q));
        correlation = correlation / correlation(1);
        pair_corr = accumarray(pair_of, double(weight) .* correlation(lag_index), [numel(pair_bin) 1]);
        model = accumarray(pair_bin, pair_corr, [n_bin 1]) ./ max(bin_count, 1);
        % Best variance for these ranges, then the weighted squared misfit
        s2 = sum(bin_count(use) .* model(use) .* bin_cov(use)) / sum(bin_count(use) .* model(use) .^ 2);
        misfit = sum(bin_count(use) .* (bin_cov(use) - s2 * model(use)) .^ 2) / sum(bin_count(use) .* bin_cov(use) .^ 2);
        model = s2 * model;
    end
end
