function [active_field, info] = solve_two_level(observations, g, active_mask_3d, steps, advection_reference, shape_patterns, z_agl, varargin)
%SOLVE_TWO_LEVEL Two-level separable inversion: coarse level and shape, fine deviation.
%   Density on active cells (MODEL_REVIEW section 5, prototype step 3):
%       rho(h, z, t) = x_1(h, t) s(h, z, t),   x_1 = P v + e,
%       s = phi_1 + sum_{k=2..K} (P c_k)(h, t) phi_k(z - smoothed ground(h)),
%   - v: the large-scale level V_L (birds km^-2) on a coarse grid of nodes
%     NODE_SPACING_M apart, P the bilinear interpolation to the 4 km columns;
%   - e: the fine-scale deviation on the columns;
%   - c_k: the shape coefficients on the same coarse grid (c = 0 is the
%     night's mean shape phi_1), phi_k the columns of SHAPE_PATTERNS scaled
%     so that sum(phi_1) dz = 1 km.
%   The ground is the DEM low-pass filtered over GROUND_SMOOTHING_M; the
%   full DEM still sets the active cells and the beam blockage.
%
%   rho is bilinear in (x_1, c), so the solver alternates N_OUTER rounds:
%   1. amount step, c fixed, linear in (v, e):
%      - data ||F rho - d||^2 (ROW_WEIGHTS optional);
%      - level: LEVEL_SMOOTHNESS^2 (LEVEL_LENGTH_M / spacing)^2 v' L_c v, an
%        intrinsic GMRF on the coarse grid (no mean, no radar prior);
%      - fine deviation: FINE_PRIOR^2 (W e)' (I + (FINE_LENGTH_M / dx)^2 L_h) (W e),
%        relative to a reference amount: W = median(ref) / max(ref, LEVEL_FLOOR)
%        from the previous round (uniform in the first round), so the
%        deviation is multiplicative (heteroscedastic) and empty areas
%        carry little structure. WEIGHT_REFERENCE "level": ref = P v;
%        "amount": ref = x_1 smoothed over WEIGHT_SMOOTHING_M (Gaussian);
%        "uniform": W = 1. FINE_AMPLITUDE "uniform" drops W from the
%        amplitude part only, FINE_PRIOR^2 (||e||^2 + (W e)' L (W e)): the
%        relative amplitude makes V_L a mean weighted by 1/level^2, which
%        sits below dense patches (+12% on observed columns) and leaves the
%        unobserved far columns low. Default "uniform" (24 October 2026);
%      - advection on the 3D field, relaxed where ADVECTION_REFERENCE
%        changes (takeoff and landing), as in solve_advection_lean.
%        ADVECTION_ON "deviation" advects only the deviation B1 e: the
%        level's own change through the night is then the source term
%        (takeoff, landing), and data-free areas follow the level instead
%        of keeping earlier, lower densities. Default "deviation";
%      - NEGATIVE_PENALTY^2 ||min(rho, 0)||^2 by FISTA after the PCG.
%   2. shape step, x_1 fixed, linear in c, one direct solve per time:
%      data plus SHAPE_PRIOR^2 / SHAPE_SD(k-1)^2 (c_k - c0_k)' (I + (SHAPE_LENGTH_M /
%      spacing)^2 L_c) (c_k - c0_k) (penalty growing with pattern order). Advection
%      is ignored in this step.
%   SHAPE_INITIAL ([nxy K-1 n_time], e.g. regional_shape.m) starts c and
%   is the mean of its prior; empty uses the night's mean shape (c = 0).
%
%   Observation model (step 4), each switch independent, updated between
%   rounds from the previous round's prediction:
%   - NEAR_RANGE "fixed" / "estimate": predicted = a(r) F rho with the
%     recovery curve a(r) = 1 - exp(-(r / r0)^k) (start NEAR_RANGE_R0_M,
%     NEAR_RANGE_K, from the measured near-range suppression); "estimate"
%     leaves the rows within NEAR_RANGE_FIT_MAX_M out of the first round,
%     then fits r0 and k (grid search) on their echo values against that
%     held-out prediction, and uses all rows with a(r) from the second
%     round on. One curve per elevation group of NEAR_RANGE_ELEV_EDGES
%     (degrees; [0 90] = one global curve). "off": a = 1. Default:
%     "estimate" per elevation group (< 1.2, 1.2-2, 2-3, >= 3 deg) within
%     50 km, decided 23 October 2026 (MODEL_REVIEW section 6, decision 4).
%   - NO_ECHO "one_sided": a no-echo gate only says predicted <= L(r), the
%     detection limit: rows predicted below L drop out, rows above it are
%     pulled to L (active set, first round all active). L(r) is
%     DETECTION_LIMIT "empirical" (the 1% quantile of the echo values by
%     range, interpolated) or "clear_air" (-19.5 dBZ at 50 km plus
%     20 log10(r / 50 km), with an 11 cm^2 cross-section), times
%     DETECTION_SCALE. "zero": no-echo gates are zeros (as before).
%   - ERROR_CV > 0: multiplicative error, sigma = ERROR_CV * max(predicted,
%     0) + ERROR_FLOOR; data weights 1 / sigma^2 normalized to median 1 on
%     the echo rows (iterative reweighting from the second round).

p = inputParser;
addParameter(p, "n_pattern", 3);
addParameter(p, "shape_sd", [0.52 0.28 0.14 0.065 0.03]);
addParameter(p, "shape_initial", []);
addParameter(p, "node_spacing_m", 40000);
addParameter(p, "level_smoothness", 1);
addParameter(p, "level_length_m", 600000);
addParameter(p, "fine_prior", 0.3);
addParameter(p, "fine_length_m", 30000);
addParameter(p, "shape_prior", 1);
addParameter(p, "shape_length_m", 300000);
addParameter(p, "level_floor", 10);
addParameter(p, "weight_reference", "level");
addParameter(p, "weight_smoothing_m", 100000);
addParameter(p, "fine_amplitude", "uniform");
addParameter(p, "advection_on", "deviation");
addParameter(p, "ground_smoothing_m", 30000);
addParameter(p, "advection", 1);
addParameter(p, "advection_relax_exponent", 2);
addParameter(p, "advection_relax_floor", 5);
addParameter(p, "negative_penalty", 3);
addParameter(p, "row_weights", []);
addParameter(p, "n_outer", 3);
addParameter(p, "near_range", "estimate");
addParameter(p, "near_range_r0_m", 10000);
addParameter(p, "near_range_k", 3.8);
addParameter(p, "near_range_fit_max_m", 50000);
addParameter(p, "near_range_elev_edges", [0 1.2 2 3 90]);
addParameter(p, "no_echo", "zero");
addParameter(p, "detection_limit", "empirical");
addParameter(p, "detection_scale", 1);
addParameter(p, "error_cv", 0);
addParameter(p, "error_floor", 1);
addParameter(p, "tol", 1e-4);
addParameter(p, "maxit", 500);
addParameter(p, "projected_maxit", 100);
addParameter(p, "projected_tol", 1e-3);
parse(p, varargin{:});
opt = p.Results;
K = opt.n_pattern;
n_time = size(observations.values, 2);

%% Columns and shape patterns ----------------------------------------------------
active_id = find(active_mask_3d(:));
n_active = numel(active_id);
[ix, iy, iz] = ind2sub(g.sz, active_id);
[column_xy, ~, column_of_cell] = unique(sub2ind(g.sz(1:2), ix, iy));
n_column = numel(column_xy);
ground_smooth = imgaussfilt(max(g.dem, 0), opt.ground_smoothing_m / mean([g.dx g.dy]));
z = squeeze(g.z);
height = z(iz) - ground_smooth(column_xy(column_of_cell));
phi = interp1(z_agl(:), shape_patterns(:, 1:K), min(max(height, z_agl(1)), z_agl(end)));

%% Coarse grid: bilinear interpolation P and graph Laplacian L_c ------------------
step_cells = round(opt.node_spacing_m / mean([g.dx g.dy]));
[column_ix, column_iy] = ind2sub(g.sz(1:2), column_xy);
node_sz = ceil((g.sz(1:2) - 1) / step_cells) + 1;
fx = (column_ix - 1) / step_cells;
fy = (column_iy - 1) / step_cells;
x0 = min(floor(fx), node_sz(1) - 2);
y0 = min(floor(fy), node_sz(2) - 2);
wx = fx - x0;
wy = fy - y0;
corner_id = [sub2ind(node_sz, x0 + 1, y0 + 1), sub2ind(node_sz, x0 + 2, y0 + 1), ...
    sub2ind(node_sz, x0 + 1, y0 + 2), sub2ind(node_sz, x0 + 2, y0 + 2)];
corner_w = [(1 - wx) .* (1 - wy), wx .* (1 - wy), (1 - wx) .* wy, wx .* wy];
[node_id, ~, node_of_corner] = unique(corner_id(:));
n_node = numel(node_id);
P = sparse(repmat((1:n_column)', 4, 1), node_of_corner, corner_w(:), n_column, n_node);
% 4-neighbour links between used nodes
node_position = zeros(prod(node_sz), 1);
node_position(node_id) = 1:n_node;
[node_ix, node_iy] = ind2sub(node_sz, node_id);
links = zeros(0, 2);
for offset = [1 0; 0 1]'
    nx = node_ix + offset(1);
    ny = node_iy + offset(2);
    inside = nx <= node_sz(1) & ny <= node_sz(2);
    neighbour = zeros(n_node, 1);
    neighbour(inside) = node_position(sub2ind(node_sz, nx(inside), ny(inside)));
    keep = neighbour > 0;
    links = [links; find(keep) neighbour(keep)]; %#ok<AGROW>
end
A_c = sparse(links(:, 1), links(:, 2), 1, n_node, n_node);
A_c = A_c + A_c';
L_c = spdiags(full(sum(A_c, 2)), 0, n_node, n_node) - A_c;
% Cells to coarse nodes, for the shape step
cell_to_node = P(column_of_cell, :);
node_lon = g.lon(min((node_ix - 1) * step_cells + 1, g.sz(1)));
node_lat = g.lat(min((node_iy - 1) * step_cells + 1, g.sz(2)));

%% Data term --------------------------------------------------------------------
% Row x time arrays are single (about 0.5 GB each on the eastern US night);
% they are cast to double one time column at a time where they meet the
% sparse operators.
F = observations.F(:, active_id);
row_keep = any(F, 2);
F = F(row_keep, :);
values = observations.values(row_keep, :);
row_range_m = double(observations.row_range_m(row_keep));
row_elev_deg = double(observations.row_elev_deg(row_keep));
if ~isempty(opt.row_weights)
    row_scale = sqrt(double(opt.row_weights(row_keep)));
    keep_weighted = row_scale > 0;
    F = row_scale(keep_weighted) .* F(keep_weighted, :);
    values = row_scale(keep_weighted) .* values(keep_weighted, :);
    row_range_m = row_range_m(keep_weighted);
    row_elev_deg = row_elev_deg(keep_weighted);
end
values = single(full(values));
obs_mask = values ~= 0;
is_no_echo = obs_mask & values < 1e-6;      % no-echo gates are stored as eps
is_echo = obs_mask & ~is_no_echo;
values(is_no_echo) = 0;
data_scales2 = 1 ./ max(sum(obs_mask, 1), 1);
data_block_norms = max(sqrt(data_scales2) .* sqrt(full(sum(F .^ 2, 2))' * obs_mask), 1);

%% Observation model: detection limit and near-range recovery ----------------------
if opt.detection_limit == "clear_air"
    detection = z_to_refl(idb(-19.5)) / 11 .* (row_range_m / 50000) .^ 2;   % birds km^-3
else
    % 1% quantile of the echo values in 10 km range bins
    edges = 0:10000:max(row_range_m) + 10000;
    row_bin = discretize(row_range_m, edges);
    bin_of_echo = repmat(row_bin, 1, n_time);
    low = accumarray(bin_of_echo(is_echo), double(values(is_echo)), [numel(edges) - 1 1], @(v) prctile(v, 1), nan);
    clear bin_of_echo
    centers = edges(1:end - 1)' + 5000;
    ok = ~isnan(low);
    detection = interp1(centers(ok), low(ok), row_range_m, "linear", "extrap");
end
detection = opt.detection_scale .* max(detection, 0);
% Recovery curve per elevation group
elev_group = discretize(row_elev_deg, opt.near_range_elev_edges);
elev_group(isnan(elev_group)) = numel(opt.near_range_elev_edges) - 1;
n_group = numel(opt.near_range_elev_edges) - 1;
near_r0 = repmat(opt.near_range_r0_m, n_group, 1);
near_k = repmat(opt.near_range_k, n_group, 1);
recovery = @(r0, k) 1 - exp(-(row_range_m ./ r0(elev_group)) .^ k(elev_group));
near_rows = find(row_range_m < opt.near_range_fit_max_m);

%% Advection links on the 3D field (as in solve_advection_lean) ------------------------
% uint32 links and single weights: half the memory of double
position = zeros(prod(g.sz), 1, "uint32");
position(active_id) = 1:n_active;
steps = steps(:);
advection_scales2 = zeros(max(n_time - 1, 0), 1);
advection_diagonal = zeros(n_active, n_time, "single");
for i_t = 1:n_time - 1
    source = position(double(steps{i_t}.source_id_valid));
    dest = position(double(steps{i_t}.dest_id_valid));
    keep = source > 0 & dest > 0;
    coeff = steps{i_t}.coeff_valid;
    if ~isscalar(coeff)
        coeff = coeff(keep);
    end
    source = source(keep);
    dest = dest(keep);
    base_norm = sqrt(2 .* sum(double(coeff) .^ 2 .* ones(nnz(keep), 1)));
    if opt.advection_relax_exponent > 0
        reference_source = double(advection_reference(source, i_t));
        reference_dest = double(advection_reference(dest, i_t + 1));
        conserved = (min(reference_source, reference_dest) + opt.advection_relax_floor) ./ ...
            (max(reference_source, reference_dest) + opt.advection_relax_floor);
        coeff = double(coeff) .* conserved .^ opt.advection_relax_exponent;
        clear reference_source reference_dest conserved
    end
    coeff2 = single(double(coeff) .^ 2 .* ones(numel(source), 1));
    advection_scales2(i_t) = (opt.advection .* max(mean(data_block_norms(i_t:i_t + 1)), 1) ./ max(base_norm, eps)) .^ 2;
    w = single(advection_scales2(i_t)) .* coeff2;
    advection_diagonal(:, i_t) = advection_diagonal(:, i_t) + accumarray(source, w, [n_active 1]);
    advection_diagonal(:, i_t + 1) = advection_diagonal(:, i_t + 1) + accumarray(dest, w, [n_active 1]);
    steps{i_t} = struct("source", source, "dest", dest, "coeff2", coeff2);
end
clear position coeff coeff2 w

%% Prior operators and term scales --------------------------------------------------
column_mask = false([g.sz(1:2) 1]);
column_mask(column_xy) = true;
L_h = bmm4d.grid.active_laplacian(column_mask, bmm4d.grid.get_water_mask(g), 0.5) / 0.25;      % unit-weight graph Laplacian
L_fine = (opt.fine_length_m / mean([g.dx g.dy])) ^ 2 * L_h;   % Q_fine = I + L_fine
R_level = (opt.level_length_m / opt.node_spacing_m) ^ 2 * L_c;
fine_scales2 = (opt.fine_prior .* data_block_norms ./ sqrt(n_column)) .^ 2;
level_scales2 = (opt.level_smoothness .* data_block_norms ./ sqrt(n_node)) .^ 2;
negative_scales2 = (opt.negative_penalty .* data_block_norms ./ sqrt(n_active)) .^ 2;
Q_shape = speye(n_node) + (opt.shape_length_m / opt.node_spacing_m) ^ 2 * L_c;

%% Initial shape coefficients on the nodes -----------------------------------------
shape = zeros(n_node, K - 1, n_time);
if ~isempty(opt.shape_initial)
    node_xy = sub2ind(g.sz(1:2), min((node_ix - 1) * step_cells + 1, g.sz(1)), min((node_iy - 1) * step_cells + 1, g.sz(2)));
    shape = opt.shape_initial(node_xy, 1:K - 1, :);
end
shape_mean = shape;

%% Alternate amount and shape steps -------------------------------------------------
X = zeros(n_node + n_column, n_time);          % [v; e] per time
predicted = [];
history = table('Size', [opt.n_outer 7], 'VariableTypes', repmat("double", 1, 7), ...
    'VariableNames', ["round" "pcg_iterations" "relres" "projected_iterations" "shape_rms" "no_echo_active_share" "seconds"]);
for i_outer = 1:opt.n_outer
    round_tic = tic;
    % Observation model from the previous round's prediction
    a = ones(numel(row_range_m), 1);
    if opt.near_range ~= "off"
        a = recovery(near_r0, near_k);
    end
    data_weight = single(obs_mask);
    target = values;
    if opt.no_echo == "one_sided"
        % Active set: no-echo gates predicted above the detection limit are pulled to it
        no_echo_limit = repmat(single(detection), 1, n_time);
        target(is_no_echo) = no_echo_limit(is_no_echo);
        if i_outer > 1
            data_weight(is_no_echo & predicted <= no_echo_limit) = 0;
        end
        clear no_echo_limit
    end
    if opt.near_range == "estimate" && i_outer == 1
        data_weight(near_rows, :) = 0;   % held out to estimate a(r)
    end
    if opt.error_cv > 0 && i_outer > 1
        error_weight = 1 ./ (opt.error_cv .* max(predicted, 0) + opt.error_floor) .^ 2;
        error_weight = error_weight ./ median(error_weight(is_echo));
        data_weight = data_weight .* error_weight;
        clear error_weight
    end

    % Basis of x_1 per time (the regional shape s on each active cell) and the data operator
    basis = zeros(n_active, n_time, "single");
    FB1 = cell(n_time, 1);
    for i_t = 1:n_time
        s = phi(:, 1);
        column_shape = P * shape(:, :, i_t);
        for k = 2:K
            s = s + phi(:, k) .* column_shape(column_of_cell, k - 1);
        end
        basis(:, i_t) = single(s);
        FB1{i_t} = a .* (F * sparse(1:n_active, column_of_cell, s, n_active, n_column));
    end
    clear s
    % Heteroscedastic weights relative to the level of the previous round
    if i_outer == 1 || opt.weight_reference == "uniform"
        W = ones(n_column, n_time);
    else
        W = median(weight_reference(:)) ./ weight_reference;
    end
    % Amplitude of the deviation: relative (as W) or unweighted, which keeps
    % V_L the plain local mean of x_1 (no pull towards the low columns)
    W_amplitude = W;
    if opt.fine_amplitude == "uniform"
        W_amplitude = ones(n_column, n_time);
    end
    op = struct("P", P, "basis", basis, "column_of_cell", column_of_cell, "n_column", n_column, ...
        "FB1", {FB1}, "data_weight", data_weight, "data_scales2", data_scales2, ...
        "L_fine", L_fine, "W", W, "W_amplitude", W_amplitude, "advect_deviation", opt.advection_on == "deviation", ...
        "fine_scales2", fine_scales2, "R_level", R_level, "level_scales2", level_scales2, ...
        "advection_scales2", advection_scales2, "steps", {steps}, "n_node", n_node, "n_time", n_time);
    clear basis

    % Right-hand side and Jacobi preconditioner
    g_data = zeros(n_column, n_time);
    diag_column = zeros(n_column, n_time);
    for i_t = 1:n_time
        weight_t = double(data_weight(:, i_t));
        g_data(:, i_t) = FB1{i_t}' * (weight_t .* double(target(:, i_t))) .* data_scales2(i_t);
        diag_column(:, i_t) = full((FB1{i_t} .^ 2)' * weight_t) .* data_scales2(i_t) ...
            + accumarray(column_of_cell, double(op.basis(:, i_t) .^ 2 .* advection_diagonal(:, i_t)), [n_column 1]);
    end
    clear target weight_t
    rhs = [P' * g_data; g_data];
    diagonal = [(P .^ 2)' * diag_column + full(diag(R_level)) .* level_scales2; ...
        diag_column + (W_amplitude .^ 2 + W .^ 2 .* full(diag(L_fine))) .* fine_scales2];
    ridge = 1e-8 * max(diagonal(:));
    op.ridge = ridge;
    inverse_diagonal = 1 ./ (diagonal + ridge);

    % Preconditioned conjugate gradient from the previous round
    r = rhs - amount_operator(X, op);
    zr = r .* inverse_diagonal;
    d = zr;
    rz = sum(r .* zr, "all");
    rhs_norm = max(norm(rhs(:)), eps);
    relres = norm(r(:)) / rhs_norm;
    for iter = 1:opt.maxit
        Ad = amount_operator(d, op);
        alpha = rz / sum(d .* Ad, "all");
        X = X + alpha .* d;
        r = r - alpha .* Ad;
        relres = norm(r(:)) / rhs_norm;
        if relres <= opt.tol
            break
        end
        zr = r .* inverse_diagonal;
        rz_next = sum(r .* zr, "all");
        d = zr + (rz_next / rz) .* d;
        rz = rz_next;
    end

    % FISTA with the penalty on negative densities
    v = randn(size(X));
    v = v ./ norm(v(:));
    for i_power = 1:15
        w = (amount_operator(v, op) + density_t(density(v, op) .* single(negative_scales2), op)) .* inverse_diagonal;
        lambda_max = norm(w(:));
        v = w ./ lambda_max;
    end
    step = 1 / (1.05 * lambda_max);
    Y = X;
    t_momentum = 1;
    for i_projected = 1:opt.projected_maxit
        gradient = amount_operator(Y, op) - rhs + density_t(min(density(Y, op), 0) .* single(negative_scales2), op);
        X_next = Y - step .* gradient .* inverse_diagonal;
        relative_step = norm(X_next(:) - X(:)) / max(norm(X(:)), eps);
        t_next = (1 + sqrt(1 + 4 * t_momentum ^ 2)) / 2;
        Y = X_next + ((t_momentum - 1) / t_next) .* (X_next - X);
        X = X_next;
        t_momentum = t_next;
        if relative_step <= opt.projected_tol
            break
        end
    end
    clear Y gradient X_next v w r zr d Ad
    x1 = P * X(1:n_node, :) + X(n_node + 1:end, :);
    level = max(P * X(1:n_node, :), opt.level_floor);
    weight_reference = level;
    if opt.weight_reference == "amount"
        % Previous amount smoothed by normalized convolution over the active columns
        sigma_cells = opt.weight_smoothing_m / mean([g.dx g.dy]);
        is_column = zeros(g.sz(1:2));
        is_column(column_xy) = 1;
        support = imgaussfilt(is_column, sigma_cells);
        for i_t = 1:n_time
            amount_grid = zeros(g.sz(1:2));
            amount_grid(column_xy) = max(x1(:, i_t), 0);
            smoothed = imgaussfilt(amount_grid, sigma_cells) ./ max(support, eps);
            weight_reference(:, i_t) = max(smoothed(column_xy), opt.level_floor);
        end
    end
    rho_round = density(X, op);
    if opt.near_range == "estimate" && i_outer == 1
        % Recovery curve per elevation group: grid search of r0 and k on the held-out near echo rows
        F_near = F(near_rows, :);
        unscaled = zeros(numel(near_rows), n_time);
        for i_t = 1:n_time
            unscaled(:, i_t) = F_near * double(rho_round(:, i_t));   % prediction without the near-range factor
        end
        clear F_near
        r_all = repmat(row_range_m(near_rows), 1, n_time);
        group_all = repmat(elev_group(near_rows), 1, n_time);
        echo_near = is_echo(near_rows, :);
        values_near = double(values(near_rows, :));
        for i_g = 1:n_group
            near = echo_near & group_all == i_g & unscaled > 0;
            u = unscaled(near);
            d_near = values_near(near);
            r_near = r_all(near);
            best = inf;
            for r0 = 2000:500:40000
                for k = 1:0.25:6
                    cost = sum((d_near - (1 - exp(-(r_near / r0) .^ k)) .* u) .^ 2);
                    if cost < best
                        best = cost;
                        near_r0(i_g) = r0;
                        near_k(i_g) = k;
                    end
                end
            end
        end
        clear unscaled r_all group_all echo_near values_near near u d_near r_near
    end
    % Prediction with this round's a(r), needed by the one-sided no-echo and error models
    if opt.no_echo == "one_sided" || opt.error_cv > 0
        predicted = zeros(numel(row_range_m), n_time, "single");
        for i_t = 1:n_time
            predicted(:, i_t) = single(a .* (F * double(rho_round(:, i_t))));
        end
    end
    clear rho_round

    % Shape step: coefficients on the nodes, one direct solve per time (all
    % rows, unobserved ones have zero weight)
    if K > 1 && i_outer < opt.n_outer
        for i_t = 1:n_time
            amount_cell = max(x1(column_of_cell, i_t), 0);
            G = cell(1, K - 1);
            for k = 2:K
                G{k - 1} = (amount_cell .* phi(:, k)) .* cell_to_node;
            end
            row_sqrt_weight = sqrt(double(data_weight(:, i_t)));
            FG = (row_sqrt_weight .* a) .* (F * [G{:}]);
            clear G
            shape_target = row_sqrt_weight .* (double(values(:, i_t)) - a .* (F * (amount_cell .* phi(:, 1))));
            if opt.no_echo == "one_sided"
                limit_t = detection .* is_no_echo(:, i_t);
                shape_target = shape_target + row_sqrt_weight .* limit_t;
            end
            normal = FG' * FG;
            information = full(mean(diag(normal(1:n_node, 1:n_node))));
            R_shape = kron(spdiags(opt.shape_sd(1) ^ 2 ./ reshape(opt.shape_sd(1:K - 1), [], 1) .^ 2, 0, K - 1, K - 1), Q_shape);
            % Prior centred on the initial shape (the radar profiles' regional shape, or the mean shape)
            prior_mean = reshape(shape_mean(:, :, i_t), [], 1);
            shape(:, :, i_t) = reshape((normal + opt.shape_prior ^ 2 * information * R_shape) \ ...
                (FG' * shape_target + opt.shape_prior ^ 2 * information * (R_shape * prior_mean)), n_node, K - 1);
            clear FG normal
        end
    end
    history(i_outer, :) = {i_outer, iter, relres, i_projected, sqrt(mean(shape .^ 2, "all")), ...
        nnz(data_weight(is_no_echo) > 0) / max(nnz(is_no_echo), 1), toc(round_tic)};
    clear FB1
end

%% Density on active cells ---------------------------------------------------------
rho = density(X, op);
solver = struct("name", "two_level_alternating", "iterations", sum(history.pcg_iterations), "history", history, ...
    "negative_mass_share", -sum(min(rho, 0), "all") / max(sum(max(rho, 0), "all"), eps));
active_field = max(rho, 0);
info = struct("controls", rmfield(opt, ["row_weights" "shape_initial"]), "solver", solver, ...
    "level_nodes", single(X(1:n_node, :)), "shape_nodes", single(shape), "node_lon", node_lon, "node_lat", node_lat, ...
    "column_amount", single(x1), "column_level", single(P * X(1:n_node, :)), "column_xy", column_xy, ...
    "near_curve", table(opt.near_range_elev_edges(1:end - 1)', opt.near_range_elev_edges(2:end)', near_r0, near_k, ...
    VariableNames=["elev_from" "elev_to" "r0_m" "k"]), "detection_bins", detection_summary(row_range_m, detection));
end

function rho = density(X, op)
% 3D density (single) from [v; e]: rho_t = s_t .* x_1(column of the cell).
x1 = single(op.P * X(1:op.n_node, :) + X(op.n_node + 1:end, :));
rho = op.basis .* x1(op.column_of_cell, :);
end

function Y = density_t(R, op)
% Adjoint of density: sum s_t .* R over the cells of each column, then [P'; I].
g1 = zeros(op.n_column, op.n_time);
for i_t = 1:op.n_time
    g1(:, i_t) = accumarray(op.column_of_cell, double(op.basis(:, i_t) .* R(:, i_t)), [op.n_column 1]);
end
Y = [op.P' * g1; g1];
end

function Y = amount_operator(X, op)
% Data, level GMRF, fine deviation GMRF and advection for [v; e].
V = X(1:op.n_node, :);
E = X(op.n_node + 1:end, :);
x1 = op.P * V + E;
g1 = zeros(size(x1));
for i_t = 1:op.n_time
    projection = (op.FB1{i_t} * x1(:, i_t)) .* double(op.data_weight(:, i_t));
    g1(:, i_t) = op.FB1{i_t}' * projection .* op.data_scales2(i_t);
end
g_advection = zeros(size(x1));
if op.n_time > 1
    % Advected 3D field: the whole density, or only its deviation from the level
    advected = x1;
    if op.advect_deviation
        advected = E;
    end
    advected = single(advected);
    rho = op.basis .* advected(op.column_of_cell, :);
    Y_rho = zeros(size(rho), "single");
    for i_t = 1:op.n_time - 1
        step = op.steps{i_t};
        weighted = single(op.advection_scales2(i_t)) .* step.coeff2 .* (rho(step.dest, i_t + 1) - rho(step.source, i_t));
        Y_rho(:, i_t) = Y_rho(:, i_t) - accumarray(step.source, weighted, [size(rho, 1) 1]);
        Y_rho(:, i_t + 1) = Y_rho(:, i_t + 1) + accumarray(step.dest, weighted, [size(rho, 1) 1]);
    end
    clear rho
    for i_t = 1:op.n_time
        g_advection(:, i_t) = accumarray(op.column_of_cell, double(op.basis(:, i_t) .* Y_rho(:, i_t)), [op.n_column 1]);
    end
end
if op.advect_deviation
    Y = [op.P' * g1; g1 + g_advection];
else
    Y = [op.P' * (g1 + g_advection); g1 + g_advection];
end
Y = Y + [(op.R_level * V) .* op.level_scales2; ...
    (op.W_amplitude .^ 2 .* E + op.W .* (op.L_fine * (op.W .* E))) .* op.fine_scales2] + op.ridge .* X;
end

function T = detection_summary(row_range_m, detection)
% Detection limit at a few ranges (nearest row), for reporting.
range_km = [10 25 50 100 150]';
limit = nan(size(range_km));
for i_r = 1:numel(range_km)
    [~, nearest] = min(abs(row_range_m - 1000 * range_km(i_r)));
    limit(i_r) = detection(nearest);
end
T = table(range_km, limit);
end
