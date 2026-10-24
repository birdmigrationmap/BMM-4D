function rings = gd_estimate_rings(log10_lo, log10_hi, range_m, beam, table_s, range_edges_m, density_grid)
%GD_ESTIMATE_RINGS Homogeneous-density estimate per range ring of a sweep.
%   Pools the gates of each ring [RANGE_EDGES_M(k), RANGE_EDGES_M(k+1))
%   over all azimuths and estimates its density with GD_ESTIMATE_DENSITY.
%   RINGS is a table with the ring centre (km), the estimate and its 95%
%   interval (birds km^-3), the number of gates and the maximum log
%   likelihood (summed over rings for calibration).

range_m = range_m(:);
n_ring = numel(range_edges_m) - 1;
center_km = (range_edges_m(1:end - 1) + range_edges_m(2:end))' / 2000;
[density, lower, upper, n_gate, log_likelihood] = deal(nan(n_ring, 1));
for k = 1:n_ring
    rows = range_m >= range_edges_m(k) & range_m < range_edges_m(k + 1);
    lo = log10_lo(rows, :);
    hi = log10_hi(rows, :);
    known = ~(isinf(lo) & isinf(hi));        % gates without any information are left out
    if nnz(known) < 50
        continue
    end
    volume = mean(beam.volume_sample_km3(range_m(rows)));
    [density(k), interval, curve] = gd_estimate_density(lo(known), hi(known), volume, table_s, density_grid);
    log_likelihood(k) = max(curve);
    lower(k) = interval(1);
    upper(k) = interval(2);
    n_gate(k) = nnz(known);
end
rings = table(center_km, density, lower, upper, n_gate, log_likelihood);
end
