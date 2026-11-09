function steps = build_radar_advection_steps(g, radar_scans, dt_seconds, varargin)
%BUILD_RADAR_ADVECTION_STEPS Interpolate radar velocity and build steps one interval at a time.
% The global advection penalty is set later in solve_advection_lean.
% Here, shoreline_crossing_coeff only edits the local link coefficients in
% the sparse advection mapping before the inverse problem is assembled.
%
% Vertical rule over terrain (ground_smoothing_m, needs active_mask_3d):
% - Inf (default): a link keeps its layer, the height above sea level
%   (build_profile_advection_steps); a link into the ground is dropped.
% - finite: a link keeps its height above the ground smoothed over
%   ground_smoothing_m (bmm4d.grid.smooth_ground), so birds follow plateaus
%   and basins and fly level over narrower ridges and valleys. The source
%   cell is shifted by the rise of the smoothed ground and linked to every
%   destination layer it overlaps, with coefficient sqrt(share of that
%   cell's air covered): a shift of a fraction of a layer is split between
%   two layers. A link whose destination is under the ground of an active
%   column goes to the lowest air cell of that column (the birds climb over
%   the terrain instead of being lost into it). Steps then carry
%   climb_valid (true for those links).

p = inputParser;
addParameter(p, "idw_neighbors", 10);
addParameter(p, "idw_power", 1.6);
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
addParameter(p, "active_mask_3d", []);
addParameter(p, "ground_smoothing_m", Inf);
parse(p, varargin{:});
terrain = isfinite(p.Results.ground_smoothing_m);
if terrain
    assert(~isempty(p.Results.active_mask_3d) && p.Results.boundary_mode == "open", ...
        "build_radar_advection_steps:TerrainOptions", "ground_smoothing_m needs active_mask_3d and the open boundary.")
    ground_smooth = bmm4d.grid.smooth_ground(g, p.Results.ground_smoothing_m);
end

n_radar = size(radar_scans, 1);
n_time = size(radar_scans, 2);
n_interval = n_time - 1;
nz = g.sz(3);
if isscalar(dt_seconds)
    dt_seconds = repmat(dt_seconds, n_interval, 1);
else
    dt_seconds = dt_seconds(:);
end

[interpolation_weights, ~] = bmm4d.prior.build_horizontal_idw_weights( ...
    g, radar_scans, ...
    n_neighbors=p.Results.idw_neighbors, ...
    power=p.Results.idw_power);

uv_profiles_radar = nan(nz, 2, n_radar, n_time);
for i_r = 1:n_radar
    for i_t = 1:n_time
        uv_profiles_radar(:, :, i_r, i_t) = bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g);
    end
end

steps = cell(n_interval, 1);
for i_interval = 1:n_interval
    uv_interval = zeros([g.sz 2]);
    profiles = 0.5 .* ( ...
        uv_profiles_radar(:, :, :, i_interval) + ...
        uv_profiles_radar(:, :, :, i_interval + 1));
    for i_component = 1:2
        profile_matrix = reshape( ...
            permute(profiles(:, i_component, :), [3 1 2]), ...
            n_radar, nz);
        uv_interval(:, :, :, i_component) = reshape( ...
            interpolation_weights * profile_matrix, g.sz);
    end
    if terrain
        steps{i_interval} = terrain_step(g, uv_interval, dt_seconds(i_interval), logical(p.Results.active_mask_3d), ...
            ground_smooth, p.Results.shoreline_crossing_coeff);
    else
        interval_step = bmm4d.advect.build_profile_advection_steps( ...
            g, uv_interval, dt_seconds(i_interval), ...
            boundary_mode=p.Results.boundary_mode, ...
            shoreline_crossing_coeff=p.Results.shoreline_crossing_coeff, ...
            active_mask_3d=p.Results.active_mask_3d);
        steps{i_interval} = interval_step{1};
    end
end
end

function step = terrain_step(g, uv_interval, dt_seconds, active_mask_3d, ground_smooth, shoreline_crossing_coeff)
% Links that keep the height above the smoothed ground, climbing into the
% lowest air cell where the destination is under the ground
nxy = prod(g.sz(1:2));
nz = g.sz(3);
z = squeeze(g.z);
grid_bottom = z(1) - g.dz / 2;
ground = max(g.dem, 0);
[has_active, lowest] = max(active_mask_3d, [], 3);
source_id = find(active_mask_3d(:));
[sx, sy, sz] = ind2sub(g.sz, source_id);
source_column = sub2ind(g.sz(1:2), sx, sy);

% Horizontal shift of each active cell (open boundary)
dest_x = sx + round(uv_interval(source_id) .* dt_seconds ./ g.dx);
dest_y = sy + round(uv_interval(source_id + prod(g.sz)) .* dt_seconds ./ g.dy);
valid = isfinite(dest_x) & isfinite(dest_y) & dest_x >= 1 & dest_x <= g.sz(1) & dest_y >= 1 & dest_y <= g.sz(2);
valid_id = find(valid);
dest_column = sub2ind(g.sz(1:2), dest_x(valid), dest_y(valid));

% The source cell's interval, shifted by the rise of the smoothed ground
rise = ground_smooth(dest_column) - ground_smooth(source_column(valid));
low = z(sz(valid)) - g.dz / 2 + rise;
high = z(sz(valid)) + g.dz / 2 + rise;
first = floor((low - grid_bottom) / g.dz + 1e-9) + 1;
last = max(ceil((high - grid_bottom) / g.dz - 1e-9), first);
n_link = last - first + 1;
link = repelem((1:numel(valid_id))', n_link);
dest_z = repelem(first, n_link) + (1:numel(link))' - repelem(cumsum(n_link) - n_link, n_link) - 1;
link_column = dest_column(link);

% Share of each destination cell's air that the interval covers
cell_bottom = max(grid_bottom + (dest_z - 1) * g.dz, ground(link_column));
cell_top = grid_bottom + dest_z * g.dz;
overlap = min(high(link), cell_top) - max(low(link), cell_bottom);
share = min(max(overlap ./ max(cell_top - cell_bottom, eps), 0), 1);
inside = dest_z >= 1 & dest_z <= nz;
keep = inside & share > 1e-6;
keep(keep) = active_mask_3d(link_column(keep) + (dest_z(keep) - 1) * nxy);

% Intervals entirely under the ground of an active column climb to its
% lowest air cell
climb = high <= ground(dest_column) & has_active(dest_column);
climb_column = dest_column(climb);
link_source = [source_id(valid_id(link(keep))); source_id(valid_id(climb))];
link_dest = [link_column(keep) + (dest_z(keep) - 1) * nxy; climb_column + (lowest(climb_column) - 1) * nxy];
coeff = [sqrt(share(keep)); ones(nnz(climb), 1)];
if shoreline_crossing_coeff ~= 1
    crossing = bmm4d.grid.find_land_water_transitions(g, mod(link_source - 1, nxy) + 1, mod(link_dest - 1, nxy) + 1);
    coeff(crossing) = coeff(crossing) * shoreline_crossing_coeff;
    is_active = coeff ~= 0;
    link_source = link_source(is_active);
    link_dest = link_dest(is_active);
    coeff = coeff(is_active);
    is_climb = [false(nnz(keep), 1); true(nnz(climb), 1)];
    is_climb = is_climb(is_active);
else
    is_climb = [false(nnz(keep), 1); true(nnz(climb), 1)];
end
step = struct("dt_seconds", dt_seconds, "boundary_mode", "open", ...
    "shoreline_crossing_coeff", shoreline_crossing_coeff, "coeff_valid", coeff, ...
    "source_id_valid", uint32(link_source), "dest_id_valid", uint32(link_dest), "climb_valid", is_climb);
end
