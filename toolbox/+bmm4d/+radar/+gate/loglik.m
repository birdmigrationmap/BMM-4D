function [log_p, known] = loglik(log10_lo, log10_hi, range_m, beam, table_s, density_grid)
%LOGLIK Log likelihood of every gate of a sweep at every density.
%   Each gate is an interval [LOG10_LO, LOG10_HI] of its summed bird
%   contribution s (BMM4D.RADAR.GATE.BOUNDS); -Inf/Inf for both ends means
%   no information. At a density, the number of birds in a gate's counting
%   volume (V_sample of BMM4D.RADAR.GATE.BEAM) is Poisson with mean
%   lambda = density * V_sample(r); the distribution of s at that lambda
%   is interpolated in log lambda in the table of
%   BMM4D.RADAR.GATE.BIRD_SUM_TABLE. The interval ends are interpolated
%   linearly in the tabulated cdf, so the +-0.25 dB of a stored value is
%   kept exactly.
%   LOG_P [n_range n_az n_density] is the log probability of each gate's
%   interval, 0 where the gate has no information (KNOWN false).

range_m = range_m(:);
[n_range, n_az] = size(log10_lo);
grid_d = density_grid(:)';
n_d = numel(grid_d);
n_s = numel(table_s.log10_s);
n_lambda = numel(table_s.log10_lambda);
known = ~(isinf(log10_lo) & isinf(log10_hi));
% Fractional position of each interval end on the s grid, and of each lambda on the lambda grid
step_s = table_s.log10_s(2) - table_s.log10_s(1);
position_of = @(x) min(max((x - table_s.log10_s(1)) / step_s + 1, 1), n_s);
lo_position = position_of(log10_lo);
hi_position = position_of(log10_hi);
step_lambda = table_s.log10_lambda(2) - table_s.log10_lambda(1);
lambda_position = min(max((log10(grid_d .* beam.volume_sample_km3(range_m)) - table_s.log10_lambda(1)) / step_lambda + 1, 1), n_lambda);
cdf_by_lambda = table_s.cdf';
log_p = zeros(n_range, n_az, n_d);
for i_r = 1:n_range
    if ~any(known(i_r, :))
        continue
    end
    cdf = interpolate_rows(cdf_by_lambda, lambda_position(i_r, :)')';
    cdf_lo = interpolate_rows(cdf, lo_position(i_r, :)');
    cdf_hi = interpolate_rows(cdf, hi_position(i_r, :)');
    cdf_lo(isinf(log10_lo(i_r, :)) & log10_lo(i_r, :) < 0, :) = 0;
    cdf_hi(isinf(log10_hi(i_r, :)) & log10_hi(i_r, :) > 0, :) = 1;
    row = log(max(cdf_hi - cdf_lo, 1e-12));
    row(~known(i_r, :), :) = 0;
    log_p(i_r, :, :) = reshape(row, 1, n_az, n_d);
end
end

function values = interpolate_rows(table, position)
% Linear interpolation of the rows of TABLE at a fractional row POSITION per output row
i0 = min(floor(position), size(table, 1) - 1);
fraction = position - i0;
values = table(i0, :) .* (1 - fraction) + table(i0 + 1, :) .* fraction;
end
