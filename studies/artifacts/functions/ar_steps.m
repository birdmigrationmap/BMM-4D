function steps = ar_steps(g, radar_scans, dt_seconds, active_mask_3d, inv, vertical, options)
%AR_STEPS Advection steps with a choice of vertical rule over terrain (dev).
%   STEPS = AR_STEPS(G, RADAR_SCANS, DT_SECONDS, ACTIVE_MASK_3D, INV, VERTICAL)
%   builds one step per interval as bmm4d.advect.build_radar_advection_steps
%   (interval wind: mean of the two frames' radar profiles, IDW to the grid;
%   horizontal shift rounded to cells; open boundary; links between active
%   cells only; INV.shoreline_crossing_coeff on land-water crossings). With
%   VERTICAL "asl" a link keeps its layer (height above sea level) and the
%   steps equal the toolbox ones. With "agl" a link keeps its height above
%   the ground: a cell k layers above the lowest active cell of its column
%   goes to the cell k layers above the lowest active cell of the
%   destination column. Birds carried over a ridge then rise with it and
%   come down after it, and the column is carried whole, except the layers
%   that would leave the top of the grid.
%   With "sigma" (terrain-following coordinate) a link keeps its relative
%   position between the column's lowest and highest active cells: source
%   cell j of n_src layers covers [j, j+1)/n_src of its column and links to
%   every destination layer that overlaps the same interval of the
%   destination column (one or more). Each link carries the ratio
%   H_src / H_dst of the two columns' air thickness (ratio_valid, for
%   ar_solve_advection_lean), so a column squeezed over a ridge keeps its
%   birds: no link goes into the ground or out of the top.
%   With "smooth" a link keeps its height above the ground smoothed with a
%   Gaussian of OPTIONS.smooth_length_m (G_L): the source cell is shifted
%   by G_L(dest) - G_L(source). Birds follow the regional terrain
%   (plateaus, basins) and fly level over the ridges and valleys smaller
%   than that scale. Infinite length: the "asl" steps.
%   With "hybrid" a link follows the ground (grid ground G) near it and is
%   level aloft: a cell at height h above its ground, with the ground
%   rising by dG = G(dest) - G(source), goes to h (1 - dG / H) above the
%   destination ground for h <= H = OPTIONS.fade_height_m, and keeps its
%   height above sea level above H (the lift fades linearly from dG at the
%   ground to 0 at H). Where dG >= H the low layer cannot be lifted and the
%   link is level.
%   In "smooth" and "hybrid" the source cell's interval is mapped, and it
%   links to every destination layer it overlaps, with coefficient
%   sqrt(share of that cell's air covered): a shift of a fraction of a
%   layer splits the link between two layers instead of rounding it away,
%   and flat ground gives the "asl" links.
%   OPTIONS.detour (default false; "asl", "smooth", "hybrid"): a link whose
%   destination layer is below the ground of an active destination column
%   (terrain above the bird's flight level) goes sideways instead of being
%   dropped: to the nearest column, perpendicular to the displacement, in
%   which that layer is air, searched on both sides up to the length of
%   the displacement (at most a 45 degree detour). Equidistant columns on
%   both sides share the link (coefficient 1 / sqrt(2) each). Steps then
%   have the field detour_valid (true for these links).
%   OPTIONS.climb (default false): such a link goes instead to the lowest
%   air cell of its destination column (the birds climb over the terrain
%   rather than being lost into it). Steps then have the field climb_valid.

arguments
    g
    radar_scans
    dt_seconds
    active_mask_3d
    inv
    vertical (1, 1) string
    options.smooth_length_m (1, 1) double = 100e3
    options.fade_height_m (1, 1) double = 1000
    options.detour (1, 1) logical = false
    options.climb (1, 1) logical = false
end

n_radar = size(radar_scans, 1);
n_interval = size(radar_scans, 2) - 1;
nz = g.sz(3);
nxy = prod(g.sz(1:2));
z = squeeze(g.z);
weights = bmm4d.prior.build_horizontal_idw_weights(g, radar_scans, ...
    n_neighbors=inv.velocity_idw_neighbors, power=inv.velocity_idw_power);
uv_radar = nan(nz, 2, n_radar, n_interval + 1);
for i_r = 1:n_radar
    for i_t = 1:n_interval + 1
        uv_radar(:, :, i_r, i_t) = bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g);
    end
end

% Source cells and their height above the lowest active cell of the column
[has_active, lowest] = max(active_mask_3d, [], 3);
lowest(~has_active) = nz + 1;
source_id = find(active_mask_3d(:));
[sx, sy, sz] = ind2sub(g.sz, source_id);
source_column = sub2ind(g.sz(1:2), sx, sy);
height_layers = sz - lowest(source_column);
ground = max(g.dem, 0);
if vertical == "sigma"
    % Active layers and air thickness (in layers) of every column
    n_layers = reshape(sum(active_mask_3d, 3), [], 1);
    air = zeros(g.sz);
    air(source_id) = bmm4d.grid.air_fraction(g, source_id);
    thickness = reshape(sum(air, 3), [], 1);
elseif vertical == "smooth"
    % Ground smoothed over the regional scale (constant for an infinite length)
    if isinf(options.smooth_length_m)
        smooth_ground = zeros(g.sz(1:2));
    else
        smooth_ground = imgaussfilt(ground, options.smooth_length_m / g.dx, Padding="replicate");
    end
    edge_low = z(sz) - g.dz / 2;
    edge_high = z(sz) + g.dz / 2;
elseif vertical == "hybrid"
    % Edges of each source cell, from its ground up, as heights above it
    source_ground = ground(source_column);
    edge_low = max(z(sz) - g.dz / 2, source_ground) - source_ground;
    edge_high = z(sz) + g.dz / 2 - source_ground;
end
grid_bottom = z(1) - g.dz / 2;

steps = cell(n_interval, 1);
for i_interval = 1:n_interval
    profiles = 0.5 .* (uv_radar(:, :, :, i_interval) + uv_radar(:, :, :, i_interval + 1));
    shift = zeros(numel(source_id), 2);
    for i_c = 1:2
        uv = weights * reshape(permute(profiles(:, i_c, :), [3 1 2]), n_radar, nz);   % nxy x nz
        cell_spacing = g.dx * (i_c == 1) + g.dy * (i_c == 2);
        shift(:, i_c) = round(uv(source_column + (sz - 1) * nxy) .* dt_seconds(i_interval) ./ cell_spacing);
    end
    dest_x = sx + shift(:, 1);
    dest_y = sy + shift(:, 2);
    valid = isfinite(dest_x) & isfinite(dest_y) & dest_x >= 1 & dest_x <= g.sz(1) & dest_y >= 1 & dest_y <= g.sz(2);
    dest_column = zeros(size(source_id));
    dest_column(valid) = sub2ind(g.sz(1:2), dest_x(valid), dest_y(valid));
    link_source = (1:numel(source_id))';
    dest_z = sz;
    link_share = [];
    blocked = false(size(source_id));
    if vertical == "agl"
        dest_z(valid) = lowest(dest_column(valid)) + height_layers(valid);
    elseif vertical == "smooth" || vertical == "hybrid"
        % Map the source cell's interval [low, high] (heights above sea level)
        valid_id = find(valid);
        dest_ground = ground(dest_column(valid));
        if vertical == "smooth"
            rise = smooth_ground(dest_column(valid)) - smooth_ground(source_column(valid));
            low = edge_low(valid) + rise;
            high = edge_high(valid) + rise;
        else
            rise = dest_ground - ground(source_column(valid));
            low = hybrid_height(edge_low(valid), rise, options.fade_height_m) + dest_ground;
            high = hybrid_height(edge_high(valid), rise, options.fade_height_m) + dest_ground;
        end
        % Every layer the interval overlaps; the one holding its centre is
        % the detour candidate when the whole interval is under the ground
        first = floor((low - grid_bottom) / g.dz + 1e-9) + 1;
        last = max(ceil((high - grid_bottom) / g.dz - 1e-9), first);
        centre = min(max(floor((0.5 * (low + high) - grid_bottom) / g.dz) + 1, first), last);
        n_link = last - first + 1;
        link_source = repelem(valid_id, n_link);
        offset = (1:numel(link_source))' - repelem(cumsum(n_link) - n_link, n_link);
        dest_z = repelem(first, n_link) + offset - 1;
        dest_column = dest_column(link_source);
        blocked = dest_z == repelem(centre, n_link) & repelem(high <= dest_ground, n_link);
        % Share of each destination cell's air that the interval covers
        cell_bottom = max(grid_bottom + (dest_z - 1) * g.dz, ground(dest_column));
        cell_top = grid_bottom + dest_z * g.dz;
        overlap = min(repelem(high, n_link), cell_top) - max(repelem(low, n_link), cell_bottom);
        link_share = min(max(overlap ./ max(cell_top - cell_bottom, eps), 0), 1);
        valid = link_share > 1e-6;
    elseif vertical == "sigma"
        % Destination layers overlapping the source cell's share of its column
        valid(valid) = n_layers(dest_column(valid)) > 0;
        n_source = n_layers(source_column(valid));
        n_dest = n_layers(dest_column(valid));
        first = floor(height_layers(valid) .* n_dest ./ n_source);
        last = ceil((height_layers(valid) + 1) .* n_dest ./ n_source) - 1;
        valid_id = find(valid);
        link_source = repelem(valid_id, last - first + 1);
        offset = (1:numel(link_source))' - repelem(cumsum(last - first + 1) - (last - first + 1), last - first + 1);
        dest_z = lowest(dest_column(link_source)) + repelem(first, last - first + 1) + offset - 1;
        valid = true(size(link_source));
        dest_column = dest_column(link_source);
    end

    % Links into the ground of an active column: the bird's level is below
    % the terrain there (detour candidates)
    if isempty(link_share)
        blocked = valid;
    end
    blocked = blocked & dest_z >= 1 & dest_z <= nz;
    blocked(blocked) = has_active(dest_column(blocked)) & dest_z(blocked) < lowest(dest_column(blocked));
    valid = valid & dest_z >= 1 & dest_z <= nz;
    dest_id = zeros(size(dest_z));
    dest_id(valid) = dest_column(valid) + (dest_z(valid) - 1) * nxy;
    valid(valid) = active_mask_3d(dest_id(valid));
    link_dest = dest_id(valid);
    if isempty(link_share)
        link_weight = ones(nnz(valid), 1);
    else
        link_weight = sqrt(link_share(valid));
    end
    link_detour = false(nnz(valid), 1);
    if options.detour && any(blocked)
        blocked_source = link_source(blocked);
        [detour_link, detour_dest, detour_weight] = detour_links(g, active_mask_3d, dest_x(blocked_source), ...
            dest_y(blocked_source), dest_z(blocked), shift(blocked_source, :));
        link_source = [link_source(valid); blocked_source(detour_link)];
        link_dest = [link_dest; detour_dest];
        link_weight = [link_weight; detour_weight];
        link_detour = [link_detour; true(numel(detour_dest), 1)];
    elseif options.climb && any(blocked)
        % Into the lowest air cell of the destination column
        blocked_column = dest_column(blocked);
        link_source = [link_source(valid); link_source(blocked)];
        link_dest = [link_dest; blocked_column + (lowest(blocked_column) - 1) * nxy];
        link_weight = [link_weight; ones(nnz(blocked), 1)];
        link_detour = [link_detour; true(nnz(blocked), 1)];
    else
        link_source = link_source(valid);
    end
    coeff = 1;
    if inv.shoreline_crossing_coeff ~= 1 || options.detour || options.climb || ~isempty(link_share)
        coeff = link_weight;
        if inv.shoreline_crossing_coeff ~= 1
            crossing = bmm4d.grid.find_land_water_transitions(g, source_column(link_source), mod(double(link_dest) - 1, nxy) + 1);
            coeff(crossing) = coeff(crossing) * inv.shoreline_crossing_coeff;
        end
    end
    steps{i_interval} = struct("dt_seconds", dt_seconds(i_interval), "boundary_mode", "open", ...
        "shoreline_crossing_coeff", inv.shoreline_crossing_coeff, "coeff_valid", coeff, ...
        "source_id_valid", uint32(source_id(link_source)), "dest_id_valid", uint32(link_dest));
    if vertical == "sigma"
        steps{i_interval}.ratio_valid = single(thickness(source_column(link_source)) ./ thickness(mod(double(link_dest) - 1, nxy) + 1));
    end
    if options.detour
        steps{i_interval}.detour_valid = link_detour;
    elseif options.climb
        steps{i_interval}.climb_valid = link_detour;
    end
end
end

function h_dest = hybrid_height(h, rise, fade_height)
% Height above the destination ground of a bird h above the source ground
% when the ground rises by RISE: lifted by rise (1 - h / H) below H, level
% above H; level everywhere where rise >= H (the low layer cannot follow)
h_dest = h - rise;
follow = rise < fade_height & h <= fade_height;
h_dest(follow) = h(follow) .* (1 - rise(follow) / fade_height);
end

function [index, dest_id, weight] = detour_links(g, active_mask_3d, dest_x, dest_y, dest_z, shift)
% Sideways destinations of blocked links, perpendicular to the displacement
% and at the same layer, nearest first, up to the displacement length.
% INDEX: the blocked link of each detour link found.
nxy = prod(g.sz(1:2));
length_cells = hypot(shift(:, 1), shift(:, 2));
perpendicular = [-shift(:, 2) shift(:, 1)] ./ max(length_cells, eps);
pending = find(length_cells >= 1);
found_index = cell(0, 1);
found_dest = cell(0, 1);
found_weight = cell(0, 1);
for distance = 1:floor(max([length_cells; 0]))
    pending = pending(length_cells(pending) >= distance);
    if isempty(pending)
        break
    end
    hit = false(numel(pending), 2);
    candidate = zeros(numel(pending), 2);
    for side = 1:2
        sign_side = 2 * side - 3;
        cx = dest_x(pending) + round(sign_side * distance * perpendicular(pending, 1));
        cy = dest_y(pending) + round(sign_side * distance * perpendicular(pending, 2));
        inside = cx >= 1 & cx <= g.sz(1) & cy >= 1 & cy <= g.sz(2);
        candidate(inside, side) = cx(inside) + (cy(inside) - 1) * g.sz(1) + (dest_z(pending(inside)) - 1) * nxy;
        hit(inside, side) = active_mask_3d(candidate(inside, side));
    end
    n_hit = sum(hit, 2);
    for side = 1:2
        found_index{end + 1, 1} = pending(hit(:, side)); %#ok<AGROW>
        found_dest{end + 1, 1} = candidate(hit(:, side), side); %#ok<AGROW>
        found_weight{end + 1, 1} = 1 ./ sqrt(n_hit(hit(:, side))); %#ok<AGROW>
    end
    pending = pending(n_hit == 0);
end
index = vertcat(found_index{:});
dest_id = vertcat(found_dest{:});
weight = vertcat(found_weight{:});
end
