function cov_fit = fit_kbessel_covariance(x, y, z, D_norm, varargin)
%FIT_KBESSEL_COVARIANCE Fit an anisotropic k-Bessel covariance model.
% Pair statistics are accumulated in chunks to avoid storing several full
% pdist vectors. This reproduces the original lag-bin means using bounded
% memory, while still evaluating every eligible sample pair.

p = inputParser;
addParameter(p, "time", ones(size(D_norm)));
addParameter(p, "max_time_lag", inf);
addParameter(p, "v_bins", [0 10 20 30 40:20:100 200 300 400 700]);
addParameter(p, "h_bins", [0 500 1000:1000:5000 7000 10000 15000 20000:10000:50000]);
addParameter(p, "min_pairs", 5);
addParameter(p, "chunk_size", 1000);
parse(p, varargin{:});

cov_fit.v = p.Results.v_bins(:)';
cov_fit.h = p.Results.h_bins(:)';
[cov_fit.H, cov_fit.V] = meshgrid(cov_fit.h(1:end-1) + diff(cov_fit.h)/2, cov_fit.v(1:end-1) + diff(cov_fit.v)/2);

time = p.Results.time;
id = isfinite(x(:)) & isfinite(y(:)) & isfinite(z(:)) & isfinite(D_norm(:)) & isfinite(time(:));
x = x(id);
y = y(id);
z = z(id);
time = time(id);
D_norm = D_norm(id);

cov_fit.var = var(D_norm, "omitnan");
D_norm = D_norm - mean(D_norm, "omitnan");

% Accumulate covariance and mean lag distance within each horizontal and
% vertical bin. Time pairs outside the selected window are ignored.
n_bins = (numel(cov_fit.v) - 1) * (numel(cov_fit.h) - 1);
sum_cov = zeros(n_bins, 1);
sum_h = zeros(n_bins, 1);
sum_v = zeros(n_bins, 1);
n = zeros(n_bins, 1);

for i_start = 1:p.Results.chunk_size:numel(D_norm)
    i = i_start:min(i_start + p.Results.chunk_size - 1, numel(D_norm));
    for j_start = i_start:p.Results.chunk_size:numel(D_norm)
        j = j_start:min(j_start + p.Results.chunk_size - 1, numel(D_norm));
        time_gap = max([min(time(i))-max(time(j)), min(time(j))-max(time(i)), 0]);
        if time_gap >= p.Results.max_time_lag
            continue
        end

        dist_h = hypot(x(i) - x(j)', y(i) - y(j)');
        dist_v = abs(z(i) - z(j)');
        dist_t = abs(time(i) - time(j)');
        cov_pair = D_norm(i) .* D_norm(j)';

        if i_start == j_start
            keep = triu(true(numel(i), numel(j)), 1);
        else
            keep = true(numel(i), numel(j));
        end
        keep = keep & dist_h < cov_fit.h(end) & dist_v < cov_fit.v(end) & dist_t < p.Results.max_time_lag;

        dist_h_pair = dist_h(keep);
        dist_v_pair = dist_v(keep);
        cov_pair = cov_pair(keep);
        h_bin = discretize(dist_h_pair, cov_fit.h);
        v_bin = discretize(dist_v_pair, cov_fit.v);
        bin = v_bin + (h_bin - 1) * (numel(cov_fit.v) - 1);
        valid = ~isnan(bin) & isfinite(cov_pair);

        sum_cov = sum_cov + accumarray(bin(valid), cov_pair(valid), [n_bins 1], @sum, 0);
        sum_h = sum_h + accumarray(bin(valid), dist_h_pair(valid), [n_bins 1], @sum, 0);
        sum_v = sum_v + accumarray(bin(valid), dist_v_pair(valid), [n_bins 1], @sum, 0);
        n = n + accumarray(bin(valid), 1, [n_bins 1], @sum, 0);
    end
end

cov_fit.emp_grid = reshape(sum_cov ./ n, size(cov_fit.H));
cov_fit.H_emp = reshape(sum_h ./ n, size(cov_fit.H));
cov_fit.V_emp = reshape(sum_v ./ n, size(cov_fit.H));
cov_fit.n = reshape(n, size(cov_fit.H));
cov_fit.n(cov_fit.n == 0) = nan;
cov_fit.emp_grid(cov_fit.n < p.Results.min_pairs) = nan;

% Fit the same anisotropic k-Bessel model used in the original script.
cov_fit.fx_dist = @(v, range_v, h, range_h) sqrt((v ./ range_v).^2 + (h ./ range_h).^2);
cov_fit.fx_cov = @(d, alpha) 1 / (2^(alpha - 1) * gamma(alpha)) .* max(d, eps).^alpha .* besselk(alpha, max(d, eps));
cov_fit.fx = @(v, h, parm) parm(1) .* cov_fit.fx_cov(cov_fit.fx_dist(v, parm(2), h, parm(3)), parm(4));

parm0 = [cov_fit.emp_grid(1,1) 100 5000 1];
if isnan(parm0(1))
    parm0(1) = cov_fit.var;
end
parm_min = [0 0 0 0]';
parm_max = [cov_fit.var 1000 20000 10]';

% Give more weight to larger covariance values and better-sampled lag bins.
err_w = @(c, n_pair) max(c, 0.1) .* log(n_pair + 10);
err1 = @(parm) cov_fit.fx(cov_fit.V_emp, cov_fit.H_emp, parm) - max(0, cov_fit.emp_grid);
err2 = @(parm) sum(err_w(cov_fit.emp_grid, cov_fit.n) .* err1(parm).^2, "all", "omitnan");
cov_fit.parm = fmincon(err2, parm0, [], [], [], [], parm_min, parm_max, [], optimoptions("fmincon", "Display", "off"));

if false
    % Check empirical lag-bin covariance against the fitted model.
    id = ~isnan(cov_fit.emp_grid);
    w = err_w(cov_fit.emp_grid(id), cov_fit.n(id));
    w = w ./ max(w) .* 500;
    figure; hold on; grid on;
    scatter3(cov_fit.H_emp(id)/1000, cov_fit.V_emp(id), cov_fit.emp_grid(id), w, cov_fit.emp_grid(id), "filled")
    surf(cov_fit.H/1000, cov_fit.V, cov_fit.fx(cov_fit.V, cov_fit.H, cov_fit.parm))
    scatter3(0, 0, cov_fit.var, max(w), "k", "filled", MarkerFaceAlpha=.5)
    xlabel("Horizontal distance (km)"); ylabel("Vertical distance (m)"); zlabel("Covariance")
    view(3)
end
end
