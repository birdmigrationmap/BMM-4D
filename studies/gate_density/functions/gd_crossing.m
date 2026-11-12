function x0 = gd_crossing(x, y, level)
%GD_CROSSING First value of X at which Y crosses LEVEL (linear), NaN if none.

x0 = nan;
k = find(diff(sign(y - level)) ~= 0 & ~isnan(diff(y)), 1);
if ~isempty(k)
    x0 = interp1(y(k:k + 1), x(k:k + 1), level);
end
end
