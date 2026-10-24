function contrast_db = local_contrast(data, range_m, az_deg, window_m, floor_value)
%LOCAL_CONTRAST Excess (dB) of each gate over its local polar median.
%   DATA is linear (density or Z) on a [range az] sweep. The median is taken
%   over about WINDOW_M along range and WINDOW_M across azimuth (wrapping
%   around north); gates below FLOOR_VALUE count as FLOOR_VALUE and missing
%   gates are ignored. High values mark compact targets that stand out from
%   their surroundings (wind turbines, towers, isolated cells).

v = 10 * log10(max(data, floor_value));
range_m = range_m(:);
n_range = max(1, round(window_m / median(diff(range_m))));
local = smoothdata(v, 1, "movmedian", n_range, "omitnan");
n_az = max(1, min(size(v, 2), round(window_m ./ (deg2rad(median(diff(az_deg))) * max(range_m, 1)))));
[groups, sizes] = findgroups(n_az);
for u = 1:numel(sizes)
    if sizes(u) > 1
        rows = groups == u;
        wrapped = local(rows, [end-sizes(u)+1:end 1:end 1:sizes(u)]);
        wrapped = smoothdata(wrapped, 2, "movmedian", sizes(u), "omitnan");
        local(rows, :) = wrapped(:, sizes(u)+1:end-sizes(u));
    end
end
contrast_db = v - local;
end
