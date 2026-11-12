function curve = as_aspect_curve(D, rows, n_bin)
%AS_ASPECT_CURVE Non-parametric aspect function from radar-pair log ratios (dev).
%   CURVE = AS_ASPECT_CURVE(D, ROWS, N_BIN) estimates log T(psi) on N_BIN
%   bins of psi = az - track (default 24 of 15 deg) from the pairs ROWS of
%   as_overlap_design: regressors I(psi_a in k) - I(psi_b in k), with the
%   pair-frame-layer fixed effects and the beam slopes of as_overlap_fit.
%   The bins are identified up to a constant: log T is centred to a mean of
%   zero over the bins (T has an azimuthal mean of about 1).
%   CURVE: psi (bin centres, deg), log_t, se (clustered by radar pair).

if nargin < 2 || isempty(rows)
    rows = true(size(D.y));
end
if nargin < 3
    n_bin = 24;
end
width = 360 / n_bin;
bin_a = min(n_bin, floor(D.psi_a(rows) / width) + 1);
bin_b = min(n_bin, floor(D.psi_b(rows) / width) + 1);
n = nnz(rows);
indicator = sparse(1:n, bin_a, 1, n, n_bin) - sparse(1:n, bin_b, 1, n, n_bin);
% Bin 1 is the reference (identified up to a constant)
E = D;
E.y = D.y(rows);
E.X = [full(indicator(:, 2:end)) D.X(rows, D.beam_columns)];
E.names = ["bin" + (2:n_bin) D.names(D.beam_columns)];
E.group = D.group(rows);
E.pair = D.pair(rows);
E.var_est = D.var_est(rows);
fit = as_overlap_fit(E, E.names);
[~, id] = ismember("bin" + (2:n_bin), fit.names);
log_t = [0 fit.beta(id)];
se = [0 fit.se(id)];
offset = mean(log_t);
curve = struct("psi", (width / 2:width:360)', "log_t", (log_t - offset)', "se", se', "fit", rmfield(fit, ["resid" "y_tilde" "weights"]));
end
