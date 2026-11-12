function best = gd_profile_max(profile, grid)
%GD_PROFILE_MAX Grid maximum of a profile log likelihood, refined by a parabola.
%   PROFILE is the log likelihood at the evenly spaced values GRID; BEST is
%   the maximising value, within one grid step of the grid maximum.

[~, i] = max(profile);
i = min(max(i, 2), numel(grid) - 1);
y = profile(i - 1:i + 1);
curvature = y(3) - 2 * y(2) + y(1);
shift = 0;
if curvature < 0
    shift = max(min(0.5 * (y(1) - y(3)) / curvature, 1), -1);
end
best = grid(i) + shift * (grid(2) - grid(1));
end
