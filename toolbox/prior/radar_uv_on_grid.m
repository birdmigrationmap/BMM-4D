function uv_grid = radar_uv_on_grid(scan, g)
%RADAR_UV_ON_GRID Velocity profile of one radar scan on the grid layers.
%   UV_GRID = RADAR_UV_ON_GRID(SCAN, G) returns a numel(g.z)-by-2 [u v]
%   profile. SCAN.uv rows are height bins above the antenna (epvvp bins of
%   g.dz starting at 0 m), while g.z is above sea level, so each profile is
%   shifted by SCAN.height. Layers below the lowest or above the highest bin
%   take the nearest bin. NaN bins stay NaN (no advection link there).
%   For ground-referenced synthetic scans (height 0 or no height field) this
%   is the identity.

uv = scan.uv;
if isfield(g, "dz")
    dz = g.dz;
else
    dz = median(diff(g.z(:)));
end
radar_height = 0;
if isfield(scan, "height")
    radar_height = scan.height;
end
z_bin = ((1:size(uv, 1))' - 0.5) .* dz;
z_above_antenna = g.z(:) - radar_height;
if isscalar(z_bin)
    uv_grid = repmat(uv, numel(z_above_antenna), 1);
    return
end
uv_grid = interp1(z_bin, uv, min(max(z_above_antenna, z_bin(1)), z_bin(end)), "linear");
end
