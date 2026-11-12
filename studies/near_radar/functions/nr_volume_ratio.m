function [rows, profile, gate_list] = nr_volume_ratio(volume, cuts, options)
%NR_VOLUME_RATIO Each sweep against the radar's split-cut profile, no solve (dev).
%   [ROWS, PROFILE] = NR_VOLUME_RATIO(VOLUME, CUTS) compares, in one
%   cleaned volume (bmm4d.radar.clean_volume: volume.dens.sweeps with
%   data, range, elev, az), the density of every sweep and range ring with
%   what a vertical profile fitted to the split cuts alone (cuts whose
%   reflectivity comes from a contiguous surveillance scan) at 20-60 km
%   predicts there. The prediction weights the profile by the beam (two-way
%   Gaussian, 0.95 deg one-way width, 4/3 earth, no birds below the
%   antenna), so wide far beams and narrow near beams are compared on equal
%   terms. A ratio below 1 on the steep cuts near the radar is the
%   near-radar deficit; the split cuts at the same range are the control.
%   CUTS is the cut table of nr_read_vcp (elev, waveform_name, split); a
%   sweep takes the type of the cut closest in elevation (within 0.1 deg).
%   Name-value options:
%     reference_range_km - ring of the profile fit (default [20 60]);
%     range_edges_km     - range rings of ROWS (default [3 5 7.5 10 12.5
%                          15 20 25 30 40 50 60 80 100]);
%     dz_m               - height bin of the profile (default 100);
%     beam_width_deg     - one-way 3 dB beam width (default 0.95);
%     smoothing          - weight of the profile's curvature penalty,
%                          relative to the mean gates per fitted ring
%                          (default 1).
%   ROWS has one row per sweep and ring: station, datetime, elev, waveform
%   ("split", "batch", "doppler", "staggered", ...), split, range_lo_km,
%   range_hi_km, n (gates), height_m (mean beam centre above the antenna),
%   density and predicted (sums over the gates). Only gates whose beam lies
%   within the heights the fit sees are counted. density ./ predicted is
%   the ratio; summing both over volumes pools them.
%   PROFILE has the fitted profile: height_m, density, covered.
%   GATE_LIST (only when asked for) has one row per counted gate: sweep
%   (index into volume.dens.sweeps), range_m, az_deg, density, predicted,
%   for per-gate tests (e.g. against the clutter filter's removed power).

arguments
    volume struct
    cuts table
    options.reference_range_km (1, 2) double = [20 60]
    options.range_edges_km (1, :) double = [3 5 7.5 10 12.5 15 20 25 30 40 50 60 80 100]
    options.dz_m (1, 1) double = 100
    options.beam_width_deg (1, 1) double = 0.95
    options.smoothing (1, 1) double = 1
end

rows = table();
profile = table();
gate_list = table();
sweeps = volume.dens.sweeps;
n_sweep = numel(sweeps);

% Cut type of each sweep from the scan strategy
waveform = strings(n_sweep, 1);
split = false(n_sweep, 1);
for i_s = 1:n_sweep
    [gap, i_cut] = min(abs(cuts.elev - sweeps(i_s).elev));
    if isempty(gap) || gap > 0.1
        waveform(i_s) = "unknown";
        continue
    end
    split(i_s) = cuts.split(i_cut);
    waveform(i_s) = cuts.waveform_name(i_cut);
    if split(i_s)
        waveform(i_s) = "split";
    end
end

% Finite gates of every sweep: range (m) and density
gates = cell(n_sweep, 1);
for i_s = 1:n_sweep
    [i_range, i_az] = find(isfinite(sweeps(i_s).data));
    gates{i_s} = [sweeps(i_s).range(i_range) sweeps(i_s).data(sub2ind(size(sweeps(i_s).data), i_range, i_az)) ...
        sweeps(i_s).az(i_az)];
end

% Beam: two-way Gaussian pattern sampled across +-1.5 one-way widths
offset_deg = linspace(-1.5, 1.5, 21) * options.beam_width_deg;
beam_weight = exp(-8 * log(2) * offset_deg .^ 2 / options.beam_width_deg ^ 2);
beam_weight = beam_weight / sum(beam_weight);
reference_m = 1000 * options.reference_range_km;
max_height = 0;
for i_s = find(split)'
    [~, top] = slant2ground(reference_m(2), sweeps(i_s).elev + offset_deg(end));
    max_height = max(max_height, top);
end
if max_height == 0
    return
end
n_bin = ceil(max_height / options.dz_m);
bin_centre = ((1:n_bin)' - 0.5) * options.dz_m;
% Beam-weighted linear-interpolation basis of the profile for one gate
% geometry: [n_gate n_bin]; parts of the beam below the antenna see no birds
basis = @(range_m, elev) beam_basis(range_m, elev, offset_deg, beam_weight, bin_centre);

% Profile fit: mean density of the split cuts per 2.5 km ring of the reference range
fit_edges = reference_m(1):2500:reference_m(2);
A = zeros(0, n_bin);
y = zeros(0, 1);
n_fit = zeros(0, 1);
for i_s = find(split)'
    g = gates{i_s};
    ring = discretize(g(:, 1), fit_edges);
    ok = ~isnan(ring);
    if ~any(ok), continue, end
    n = accumarray(ring(ok), 1, [numel(fit_edges) - 1 1]);
    mean_range = accumarray(ring(ok), g(ok, 1), size(n)) ./ n;
    mean_density = accumarray(ring(ok), g(ok, 2), size(n)) ./ n;
    has = n > 0;
    A = [A; basis(mean_range(has), sweeps(i_s).elev)]; %#ok<AGROW>
    y = [y; mean_density(has)]; %#ok<AGROW>
    n_fit = [n_fit; n(has)]; %#ok<AGROW>
end
if numel(y) < 4
    return
end
% Non-negative least squares, rings weighted by their gate count, with a
% curvature penalty that spreads the profile over bins the beams average
curvature = diff(eye(n_bin), 2);
lambda = options.smoothing * mean(n_fit);
fitted = lsqnonneg([sqrt(n_fit) .* A; sqrt(lambda) * curvature], [sqrt(n_fit) .* y; zeros(n_bin - 2, 1)]);
% Heights the fit sees: bins with at least 1% of the largest beam weight
sensitivity = (n_fit' * A)';
covered = sensitivity >= 0.01 * max(sensitivity);
covered(1:find(covered, 1) - 1) = true;   % from the antenna up to the first seen bin
covered(find(covered, 1, "last") + 1:end) = false;
profile = table(bin_centre, fitted, covered, VariableNames=["height_m" "density" "covered"]);

% Each sweep and ring against the beam-weighted prediction
edges_m = 1000 * options.range_edges_km;
n_ring = numel(edges_m) - 1;
parts = cell(n_sweep, 1);
gate_parts = cell(n_sweep, 1);
for i_s = 1:n_sweep
    g = gates{i_s};
    ring = discretize(g(:, 1), edges_m);
    g = g(~isnan(ring), :);
    ring = ring(~isnan(ring));
    if isempty(g), continue, end
    B = basis(g(:, 1), sweeps(i_s).elev);
    % Keep gates whose beam weight lies (95%) on covered heights
    inside = B(:, covered) * ones(nnz(covered), 1) >= 0.95 * sum(B, 2) & sum(B, 2) > 0;
    if ~any(inside), continue, end
    predicted = B(inside, :) * fitted;
    [~, centre_height] = slant2ground(g(inside, 1), sweeps(i_s).elev);
    r = ring(inside);
    n = accumarray(r, 1, [n_ring 1]);
    k = find(n > 0);
    m = numel(k);
    sum_of = @(values) accumarray(r, values, [n_ring 1]);
    height_sum = sum_of(centre_height);
    density = sum_of(g(inside, 2));
    predicted_sum = sum_of(predicted);
    parts{i_s} = table(repmat(string(volume.station), m, 1), repmat(volume.datetime, m, 1), ...
        repmat(round(sweeps(i_s).elev, 2), m, 1), repmat(waveform(i_s), m, 1), repmat(split(i_s), m, 1), ...
        options.range_edges_km(k)', options.range_edges_km(k + 1)', n(k), height_sum(k) ./ n(k), ...
        density(k), predicted_sum(k), ...
        VariableNames=["station" "datetime" "elev" "waveform" "split" "range_lo_km" "range_hi_km" "n" ...
        "height_m" "density" "predicted"]);
    if nargout > 2
        gate_parts{i_s} = table(repmat(i_s, nnz(inside), 1), g(inside, 1), g(inside, 3), g(inside, 2), predicted, ...
            VariableNames=["sweep" "range_m" "az_deg" "density" "predicted"]);
    end
end
rows = vertcat(parts{:});
if nargout > 2
    gate_list = vertcat(gate_parts{:});
end
end

function B = beam_basis(range_m, elev, offset_deg, beam_weight, bin_centre)
% Beam-weighted interpolation weights of the profile bins for each gate
B = zeros(numel(range_m), numel(bin_centre));
for k = 1:numel(offset_deg)
    [~, height_m] = slant2ground(range_m(:), elev + offset_deg(k));
    % Linear interpolation between bin centres; below the antenna: no birds
    position = height_m / (bin_centre(2) - bin_centre(1)) + 0.5;
    lower = floor(position);
    frac = position - lower;
    for shift = 0:1
        j = lower + shift;
        w = beam_weight(k) * (shift * frac + (1 - shift) * (1 - frac));
        ok = j >= 1 & j <= numel(bin_centre) & height_m >= 0;
        B = B + sparse(find(ok), j(ok), w(ok), numel(range_m), numel(bin_centre));
    end
end
B = full(B);
end
