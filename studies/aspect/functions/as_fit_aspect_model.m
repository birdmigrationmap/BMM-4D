function aspect = as_fit_aspect_model(D, rows, radar_height, form)
%AS_FIT_ASPECT_MODEL Shared aspect function from the radar overlap (dev).
%   ASPECT = AS_FIT_ASPECT_MODEL(D, ROWS, RADAR_HEIGHT, FORM) fits on the
%   pairs ROWS of as_overlap_design the parameters of as_aspect_template,
%   with the pair-frame-layer fixed effects and the beam slopes of
%   as_overlap_fit. FORM:
%   - "symmetric": one side-on amplitude b2, no front/back term (b1 = 0);
%   - "constant": one b2 and one front/back term b1;
%   - "height": b2 proportional to the height above the antennas
%     (b2 = b2_per_km h, zero at the ground), one b1;
%   - "height_speed": b2 at h0 = 2 km and v0 = 9 m/s with slopes per km
%     and per m/s of the birds' speed, one b1.
%   Height: cell height minus the mean antenna height of the two radars
%   (RADAR_HEIGHT, m, by radar id); speed and track: mean of the two
%   radars' profiles at the cell height. ASPECT.se: standard errors of the
%   fitted parameters (clustered by radar pair).

height_km = (D.height_asl_m - (radar_height(D.radar_a) + radar_height(D.radar_b)) / 2) / 1000;
c2 = D.X(:, 1);
c1 = D.X(:, 3);
beam = D.X(:, D.beam_columns);
aspect = struct("b2", 0, "b2_per_km", 0, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, "b1", 0, "form", string(form));
switch form
    case "symmetric"
        names = "b2";
        X = c2;
    case "constant"
        names = ["b2" "b1"];
        X = [c2 c1];
    case "height"
        names = ["b2_per_km" "b1"];
        X = [c2 .* height_km c1];
    case "height_speed"
        aspect.h0_m = 2000;
        aspect.v0_ms = 9;
        names = ["b2" "b2_per_km" "b2_per_ms" "b1"];
        X = [c2, c2 .* (height_km - 2), c2 .* (D.speed - 9), c1];
end
E = D;
E.X = [X beam];
E.names = [names D.names(D.beam_columns)];
fit = as_overlap_fit(E, E.names, rows);
for i = 1:numel(names)
    aspect.(names(i)) = fit.beta(i);
end
aspect.se = fit.se(1:numel(names));
end
