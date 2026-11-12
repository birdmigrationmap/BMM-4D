function fit = gd_fit_range_bias(density, weight, range_m, elevation_deg, scan_id, antenna_agl_m, beam, varargin)
%GD_FIT_RANGE_BIAS Vertical profile per scan and one bias per range bin.
%   Each observation i is the mean DENSITY over the gates of one range bin
%   of one sweep (RANGE_M, ELEVATION_DEG) of scan SCAN_ID, with WEIGHT its
%   number of gates. It is modelled as
%       density_i = b(k_i) * sum_z W_i(z) rho_s(z),
%   with rho_s the vertical profile of scan s in layers of LAYER_EDGES_M
%   above ground (zero below ground and above the top), W_i the share of the
%   two-way beam (f^4, Gaussian in elevation with sd width / (4 sqrt(ln2)))
%   in each layer at the bin's range, beam height with the 4/3 earth radius,
%   and b(k) a bias per bin of BIN_EDGES_KM, the same for all scans and cuts.
%   Different cuts see the same layer at different ranges, which separates
%   b from rho. Fitted by multiplicative (Richardson-Lucy) updates; each scan
%   is weighted by the inverse of its mean density so that dense nights do
%   not dominate. b is scaled to 1 over REFERENCE_KM. FIT has bin_km, bias,
%   layer_mid_m, profile [n_layer n_scan], and per observation the beam
%   height_m, predicted density and cut_ratio (observed / predicted without b).
%   GROUP (optional, one positive integer per observation) gives each group
%   of cuts its own bias b(k, g), sharing the profiles; bias is then
%   [n_bin n_group], all scaled by group 1's mean over REFERENCE_KM, so the
%   other groups keep their level relative to group 1. With
%   NORMALIZE_GROUPS true, every group's bias is held at mean 1 over
%   REFERENCE_KM during the fit (all cuts agree far from the radar).
%   BIAS_BOUNDS [n_bin n_group 2] (optional) keeps each bias within
%   [lower upper] at every update (e.g. the share of the echo the clutter
%   filter left, GD_CFP_POWER; equal bounds fix it). With bounds, give
%   group 1 a fixed reference bin so that the scale is set.

p = inputParser;
addParameter(p, "layer_edges_m", 0:100:4000);
addParameter(p, "bin_edges_km", [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80]);
addParameter(p, "reference_km", [30 60]);
addParameter(p, "n_iteration", 300);
addParameter(p, "group", []);
addParameter(p, "normalize_groups", false);
addParameter(p, "bias_bounds", []);
parse(p, varargin{:});
opt = p.Results;

group = opt.group(:);
if isempty(group)
    group = ones(numel(density), 1);
end
density = density(:); weight = weight(:); range_m = range_m(:);
elevation_deg = elevation_deg(:); scan_id = scan_id(:); antenna_agl_m = antenna_agl_m(:);
keep = weight > 0 & ~isnan(density) & group > 0;
[density, weight, range_m, elevation_deg, scan_id, antenna_agl_m, group] = deal(density(keep), weight(keep), ...
    range_m(keep), elevation_deg(keep), scan_id(keep), antenna_agl_m(keep), group(keep));
[scans, ~, i_scan] = unique(scan_id);
n_scan = numel(scans);
edges = opt.layer_edges_m(:)';
n_layer = numel(edges) - 1;
bin_edges = opt.bin_edges_km(:)';
n_bin = numel(bin_edges) - 1;
i_bin = discretize(range_m / 1000, bin_edges);
keep = ~isnan(i_bin);
[density, weight, range_m, elevation_deg, i_scan, antenna_agl_m, i_bin, group] = deal(density(keep), weight(keep), ...
    range_m(keep), elevation_deg(keep), i_scan(keep), antenna_agl_m(keep), i_bin(keep), group(keep));
% One bias per range bin and group
n_group = max(group);
i_bias = sub2ind([n_bin n_group], i_bin, group);
n_obs = numel(density);

% Beam share per layer
earth_radius_m = 4 / 3 * 6371000;
height_m = antenna_agl_m + range_m .* sind(elevation_deg) + range_m .^ 2 / (2 * earth_radius_m);
height_sd_m = range_m * beam.width_rad / (4 * sqrt(log(2)));
W = normcdf((edges(2:end) - height_m) ./ height_sd_m) - normcdf((edges(1:end - 1) - height_m) ./ height_sd_m);

% Weights: gates, divided by the scan's mean density
scan_mean = accumarray(i_scan, density .* weight) ./ accumarray(i_scan, weight);
w = weight ./ max(scan_mean(i_scan), 1);

% Multiplicative updates
bin_mid_km = (bin_edges(1:end - 1) + bin_edges(2:end)) / 2;
reference_bin = bin_mid_km >= opt.reference_km(1) & bin_mid_km <= opt.reference_km(2);
profile = repmat(scan_mean', n_layer, 1);
bias = ones(n_bin * n_group, 1);
for it = 1:opt.n_iteration
    profile_obs = profile(:, i_scan)';                  % [n_obs n_layer]
    clean = sum(W .* profile_obs, 2);                    % prediction without bias
    predicted = bias(i_bias) .* clean;
    ratio = density ./ max(predicted, 1e-9);
    % Profile: per scan and layer
    numerator = accumarray([repmat(i_scan, n_layer, 1) kron((1:n_layer)', ones(n_obs, 1))], ...
        reshape(w .* bias(i_bias) .* W .* ratio, [], 1), [n_scan n_layer])';
    denominator = accumarray([repmat(i_scan, n_layer, 1) kron((1:n_layer)', ones(n_obs, 1))], ...
        reshape(w .* bias(i_bias) .* W, [], 1), [n_scan n_layer])';
    profile = profile .* numerator ./ max(denominator, 1e-12);
    % Bias: per bin
    clean = sum(W .* profile(:, i_scan)', 2);
    bias = accumarray(i_bias, w .* density, [n_bin * n_group 1]) ./ accumarray(i_bias, w .* clean, [n_bin * n_group 1]);
    bias(isnan(bias)) = 1;
    if ~isempty(opt.bias_bounds)
        bias = min(max(bias, reshape(opt.bias_bounds(:, 1:n_group, 1), [], 1)), reshape(opt.bias_bounds(:, 1:n_group, 2), [], 1));
    end
    if opt.normalize_groups
        bias = reshape(bias, n_bin, n_group);
        bias = bias ./ mean(bias(reference_bin, :), 1, "omitnan");
        bias = bias(:);
    end
end
% Bins without observations have no bias estimate
bias(accumarray(i_bias, w, [n_bin * n_group 1]) == 0) = nan;
bias = reshape(bias, n_bin, n_group);
bin_km = (bin_edges(1:end - 1) + bin_edges(2:end)) / 2;
in_reference = bin_km >= opt.reference_km(1) & bin_km <= opt.reference_km(2);
scale = mean(bias(in_reference, 1), "omitnan");
bias = bias / scale;
profile = profile * scale;
predicted = bias(i_bias) .* sum(W .* profile(:, i_scan)', 2);
fit = struct("bin_km", bin_km(:), "bias", bias, "layer_mid_m", (edges(1:end - 1) + edges(2:end))' / 2, ...
    "profile", profile, "scans", scans, "predicted", predicted, "density", density, "weight", weight, ...
    "range_m", range_m, "elevation_deg", elevation_deg, "height_m", height_m, "i_scan", i_scan, "i_bin", i_bin, ...
    "group", group, "cut_ratio", density ./ (predicted ./ bias(i_bias)));
end
