function T = ar_prior_spectrum()
%AR_PRIOR_SPECTRUM Correlation implied by each prior and smoothness combination (dev).
% Used by PRIOR_OPTIONS.md (section 3). Seconds; no data.
% Periodic 256 x 256 x 32 lattice, symbols of the solver operators with the
% solver's Frobenius scaling (the data-block factor is common and cancels).
nx = 256; nz = 32;
w = 2 * pi * (0:nx - 1) / nx;
wz = 2 * pi * (0:nz - 1) / nz;
[kx, ky, kz] = ndgrid(w, w, wz);
lh = 4 - 2 * cos(kx) - 2 * cos(ky);   % unit horizontal graph Laplacian
lv = 2 - 2 * cos(kz);                 % unit vertical graph Laplacian
dx_km = 4; dz_m = 100;

cfg = {
    % name, lambda_p, lambda_s, coeff_horver, gp_length_km, alpha2 [lr_km lz_m]
    "paper static single radar (0.01/0.03)", 0.01, 0.03, 0, 0, []
    "paper synthetic mosaic (0.1/0.1)", 0.1, 0.1, 0, 0, []
    "paper Midwest (0.1/0.3)", 0.1, 0.3, 0, 0, []
    "paper Key West (0.1/0.1, horver -0.5)", 0.1, 0.1, -0.5, 0, []
    "v2 production (0.3/0.3)", 0.3, 0.3, 0, 0, []
    "dev GP 100 km + smoothness 0.3", 0.3, 0.3, 0, 100, []
    "dev GP 100 km, no smoothness", 0.3, 0, 0, 100, []
    "SPDE alpha=2, 20 km x 190 m", 0.3, 0, 0, 0, [20 190]
    "SPDE alpha=2, 50 km x 190 m", 0.3, 0, 0, 0, [50 190]
    };
rows = cell(size(cfg, 1), 1);
for i = 1:size(cfg, 1)
    [name, lp, ls, c, gp_km, a2] = cfg{i, :};
    Lsym = ((1 + c) * lh + (1 - 2 * c) * lv) / 6;
    rowsq = ((1 + c) * 4 / 6 + (1 - 2 * c) * 2 / 6) ^ 2 + 4 * ((1 + c) / 6) ^ 2 + 2 * ((1 - 2 * c) / 6) ^ 2;
    if ~isempty(a2)
        K = 1 + (a2(1) / dx_km) ^ 2 * lh + (a2(2) / dz_m) ^ 2 * lv;
        Qp = K .^ 2;
        Qp = Qp * mean(1 ./ Qp, "all");     % same marginal variance as the white prior
        Qp = lp ^ 2 * Qp;
    elseif gp_km > 0
        a = (gp_km / dx_km) ^ 2;
        Qp = 1 + a * lh;
        Qp = Qp * mean(1 ./ Qp, "all");
        Qp = lp ^ 2 * Qp;
    else
        Qp = lp ^ 2 * ones(size(lh));
    end
    Q = Qp + ls ^ 2 * Lsym .^ 2 / rowsq;
    C = real(ifftn(1 ./ Q));
    c0 = C(1, 1, 1);
    rx = squeeze(C(1:nx / 2, 1, 1)) / c0;
    rz = squeeze(C(1, 1, 1:nz / 2)) / c0;
    lag = @(km) rx(round(km / dx_km) + 1);
    rows{i} = table(name, lag(4), lag(20), lag(52), lag(100), efold(rz) * dz_m, min(rx), sqrt(c0) * lp, ...
        max(Q, [], "all") / min(Q, [], "all"), ...
        VariableNames=["prior" "r_4km" "r_20km" "r_52km" "r_100km" "efold_v_m" "min_corr_h" "sd_vs_white" "spectral_cond"]);
end
T = vertcat(rows{:});
end

function d = efold(r)
% Lag (cells) where the correlation first falls below 1/e, linear interpolation
j = find(r < exp(-1), 1);
if isempty(j)
    d = inf;
    return
end
d = (j - 2) + (r(j - 1) - exp(-1)) / (r(j - 1) - r(j));
end
