function radar_fit = as_fit_radar_offsets(D, rows, radar_height, aspect)
%AS_FIT_RADAR_OFFSETS Per-radar side-on amplitude around a shared model (dev).
%   RADAR_FIT = AS_FIT_RADAR_OFFSETS(D, ROWS, RADAR_HEIGHT, ASPECT) fits, for
%   each radar on its overlap pairs ROWS of as_overlap_design, an offset of
%   the side-on amplitude b2 to the shared model ASPECT (as_fit_aspect_model,
%   used as an offset of the log ratios; the other radar of each pair at the
%   shared model), with errors clustered by pair and frame. With the track
%   nearly the same everywhere in one night, the offset of one radar also
%   takes any fixed azimuthal pattern of that radar (partial blockage,
%   calibration by sector) on the side-on axis, so it is not the aspect
%   alone. RADAR_FIT: estimate, se, and tau (sd of the offsets beyond
%   their standard errors).

height_km = (D.height_asl_m - (radar_height(D.radar_a) + radar_height(D.radar_b)) / 2) / 1000;
b2_model = max(0, aspect.b2 + aspect.b2_per_km * (height_km - aspect.h0_m / 1000) + ...
    aspect.b2_per_ms * (D.speed - aspect.v0_ms));
offset_y = b2_model .* D.X(:, 1) + aspect.b1 * D.X(:, 3);
[~, ~, pair_frame] = unique([D.pair D.frame], "rows");
n_radar = max([D.radar_a; D.radar_b]);
[estimate, se] = deal(nan(n_radar, 1));
F = D;
F.y = D.y - offset_y;
F.pair = pair_frame;
F.names = ["own" D.names(D.beam_columns)];
for r = 1:n_radar
    rows_r = rows & (D.radar_a == r | D.radar_b == r);
    if nnz(rows_r) < 1000
        continue
    end
    % Own side-on regressor: +(-cos 2 psi) as radar a, -(-cos 2 psi) as radar b
    own = zeros(size(D.y));
    own(D.radar_a == r) = -cosd(2 * D.psi_a(D.radar_a == r));
    own(D.radar_b == r) = cosd(2 * D.psi_b(D.radar_b == r));
    F.X = [own, D.X(:, D.beam_columns)];
    fit_r = as_overlap_fit(F, F.names, rows_r);
    estimate(r) = fit_r.beta(1);
    se(r) = fit_r.se(1);
end
ok = isfinite(estimate);
radar_fit = struct("estimate", estimate, "se", se, "tau", sqrt(max(0, var(estimate(ok)) - mean(se(ok) .^ 2))));
end
