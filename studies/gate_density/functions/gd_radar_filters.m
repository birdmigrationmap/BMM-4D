function [eta_out, cfp_code] = gd_radar_filters(eta, range_m, mode, varargin)
%GD_RADAR_FILTERS Simple stand-ins for the radar's clutter processing.
%   The WSR-88D rules are not public in detail, so these are assumptions
%   meant to test the estimator, not to reproduce the radar exactly.
%   MODE "point": a gate more than POINT_DB above both range neighbours and
%   both azimuth neighbours is a point target; it is replaced by the mean
%   of these four neighbours and flagged with CFP code 1.
%   MODE "clutter": within FILTER_RANGE_M, every gate is clutter filtered
%   (CFP code >= 8). Ground clutter (lognormal, median CLUTTER_DB above
%   1 bird km^-3 at 2 km, falling 20 dB per decade of range) is removed
%   with a random share of the bird power, uniform in [0, LOSS_MAX], where
%   LOSS_MAX falls linearly from LOSS_NEAR at 2 km to 0 at FILTER_RANGE_M.
%   CFP = 10 log10((birds + clutter) / output), as in the Level II files.
%   ETA_OUT is the output reflectivity (cm^2 km^-3); CFP_CODE is 0
%   elsewhere.

p = inputParser;
addParameter(p, "point_db", 8);
addParameter(p, "filter_range_m", 25000);
addParameter(p, "clutter_db", 30);
addParameter(p, "loss_near", 0.8);
parse(p, varargin{:});
opt = p.Results;
range_m = range_m(:);
eta_out = eta;
cfp_code = zeros(size(eta), "uint8");
switch mode
    case "point"
        level = 10 * log10(max(eta, 1e-6));
        neighbours = cat(3, circshift(level, 1, 1), circshift(level, -1, 1), circshift(level, 1, 2), circshift(level, -1, 2));
        is_point = all(level - neighbours > opt.point_db, 3);
        is_point([1 end], :) = false;
        neighbour_mean = mean(cat(3, circshift(eta, 1, 1), circshift(eta, -1, 1), circshift(eta, 1, 2), circshift(eta, -1, 2)), 3);
        eta_out(is_point) = neighbour_mean(is_point);
        cfp_code(is_point) = 1;
    case "clutter"
        rows = range_m < opt.filter_range_m;
        % 1 bird km^-3 is 11 cm^2 km^-3
        clutter_median = 11 * 10 .^ ((opt.clutter_db - 20 * log10(range_m(rows) / 2000)) / 10);
        clutter = clutter_median .* exp(1.5 * randn(nnz(rows), size(eta, 2)));
        loss_max = opt.loss_near * max(0, 1 - (range_m(rows) - 2000) / (opt.filter_range_m - 2000));
        kept = eta(rows, :) .* (1 - loss_max .* rand(nnz(rows), size(eta, 2)));
        removed_db = 10 * log10((eta(rows, :) + clutter) ./ max(kept, 1e-6));
        eta_out(rows, :) = kept;
        cfp_code(rows, :) = uint8(min(max(round(removed_db), 0), 247) + 8);
end
end
