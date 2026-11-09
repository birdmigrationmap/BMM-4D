function ground_smooth = smooth_ground(g, length_m)
%SMOOTH_GROUND Grid ground smoothed over a horizontal length.
%   GROUND_SMOOTH = BMM4D.GRID.SMOOTH_GROUND(G, LENGTH_M) is [nx ny]: the
%   grid ground max(g.dem, 0) smoothed with a Gaussian of standard
%   deviation LENGTH_M (edges replicated). The advection links and the
%   radar-profile prior keep the height above this surface, so birds follow
%   plateaus and basins and fly level over ridges and valleys narrower than
%   LENGTH_M. An infinite LENGTH_M gives zeros: height above sea level.

if isinf(length_m)
    ground_smooth = zeros(g.sz(1:2));
    return
end
ground_smooth = imgaussfilt(max(g.dem, 0), length_m / g.dx, Padding="replicate");
end
