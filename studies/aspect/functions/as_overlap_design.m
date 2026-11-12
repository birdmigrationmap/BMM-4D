function D = as_overlap_design(pairs, method, axis)
%AS_OVERLAP_DESIGN Log ratios and aspect regressors of radar pairs (dev).
%   D = AS_OVERLAP_DESIGN(PAIRS, METHOD) for the pairs of as_overlap_pairs
%   and the density version METHOD ("raw", "ell", "fix", or a two-column
%   [a b] matrix of densities) returns the pairs with both densities > 0:
%   - y: log density of radar a minus radar b (the same cell, same time);
%   - group: pair and frame (fixed effect: calibration, per-volume scale);
%   - pair: pair of radars (night-constant calibration);
%   - psi_a, psi_b: azimuth from each radar minus the birds' track (deg);
%     track from the mean of the two radars' velocity at the cell height;
%   - X: regressors, differences a - b of
%       c2 = -cos 2 psi (side-on positive: a > 0 means side-on brighter),
%       s2 = sin 2 psi (axis rotated from track +- 90 deg),
%       c1 = cos psi (tail-on brighter than head-on: psi = 0 is tail-on,
%            the bird flies away from the radar),
%       s1 = sin psi (left / right),
%       c4 = cos 4 psi,
%       beam = log(range_a / range_b) (beam broadening);
%   - var_est: variance of y from the estimator sd of each cell mean.

if isnumeric(method)
    d_a = double(method(:, 1));
    d_b = double(method(:, 2));
else
    d_a = double(pairs.(method + "_a"));
    d_b = double(pairs.(method + "_b"));
end
u = double(pairs.u_a + pairs.u_b) / 2;
v = double(pairs.v_a + pairs.v_b) / 2;
keep = d_a > 0 & d_b > 0 & isfinite(d_a) & isfinite(d_b) & isfinite(u) & isfinite(v) & hypot(u, v) > 0;
track = atan2d(u(keep), v(keep));
% AXIS (optional, one value per pair, deg): flight direction used instead of the track
if nargin > 2 && ~isempty(axis)
    keep = keep & isfinite(axis);
    track = axis(keep);
end
D = struct();
D.y = log(d_a(keep)) - log(d_b(keep));
D.keep = keep;
D.frame = double(pairs.frame(keep));
D.radar_a = double(pairs.radar_a(keep));
D.radar_b = double(pairs.radar_b(keep));
[~, ~, D.pair] = unique([D.radar_a D.radar_b], "rows");
% Fixed effect per pair, frame and 500 m layer: calibration, per-volume
% scale and the different beam smoothing of the vertical profile
D.layer = floor(double(pairs.height_asl_m(keep)) / 500);
[~, ~, D.group] = unique([D.pair D.frame D.layer], "rows");
D.psi_a = mod(double(pairs.az_a(keep)) - track, 360);
D.psi_b = mod(double(pairs.az_b(keep)) - track, 360);
D.speed = hypot(u(keep), v(keep));
D.track = mod(track, 360);
D.height_asl_m = double(pairs.height_asl_m(keep));
D.range_a = double(pairs.range_m_a(keep));
D.range_b = double(pairs.range_m_b(keep));
D.elev_a = double(pairs.elev_a(keep));
D.elev_b = double(pairs.elev_b(keep));
D.density = sqrt(d_a(keep) .* d_b(keep));
if ismember("sd2_a", pairs.Properties.VariableNames)
    D.var_est = double(pairs.sd2_a(keep)) ./ double(pairs.n_a(keep)) ./ d_a(keep) .^ 2 + ...
        double(pairs.sd2_b(keep)) ./ double(pairs.n_b(keep)) ./ d_b(keep) .^ 2;
else
    % No estimator sd (raw reflectivity): uniform weights
    D.var_est = zeros(nnz(keep), 1);
end
% Beam broadening: log range ratio, its slope free in each 500 m layer
beam = log(D.range_a ./ D.range_b);
layers = 0:7;
beam_layers = (D.layer == layers) .* beam;
D.names = ["c2" "s2" "c1" "s1" "c4" "beam_" + layers];
D.X = [-(cosd(2 * D.psi_a) - cosd(2 * D.psi_b)), sind(2 * D.psi_a) - sind(2 * D.psi_b), ...
    cosd(D.psi_a) - cosd(D.psi_b), sind(D.psi_a) - sind(D.psi_b), ...
    cosd(4 * D.psi_a) - cosd(4 * D.psi_b), beam_layers];
D.beam_columns = 6:size(D.X, 2);
end
