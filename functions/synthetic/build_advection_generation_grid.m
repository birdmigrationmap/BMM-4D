function [generation_grid, crop] = build_advection_generation_grid(g, uv_profile, time_offsets, dt_seconds, varargin)
%BUILD_ADVECTION_GENERATION_GRID Expand the inversion g for synthetic advection.

p = inputParser;
addParameter(p, "safety_cells", 1);
parse(p, varargin{:});

center_id = find(time_offsets == 0, 1);
if isempty(center_id)
    error("build_advection_generation_grid:MissingReferenceTime", ...
        "time_offsets must include 0 for the reference field.");
end

n_time = numel(time_offsets);
cumulative_x = zeros(g.sz(3), n_time);
cumulative_y = zeros(g.sz(3), n_time);
shift_profiles = get_interval_shift_profiles(uv_profile, g, n_time, dt_seconds);

for i_t = center_id + 1:n_time
    cumulative_x(:, i_t) = cumulative_x(:, i_t - 1) + shift_profiles(:, 1, i_t - 1);
    cumulative_y(:, i_t) = cumulative_y(:, i_t - 1) + shift_profiles(:, 2, i_t - 1);
end

for i_t = center_id - 1:-1:1
    cumulative_x(:, i_t) = cumulative_x(:, i_t + 1) - shift_profiles(:, 1, i_t);
    cumulative_y(:, i_t) = cumulative_y(:, i_t + 1) - shift_profiles(:, 2, i_t);
end

safety_cells = p.Results.safety_cells;
% The reference field is defined at the center time step. To generate
% earlier and later fields without clipping, the padded domain must cover
% the reference-time source support required by every cropped time slice.
pad_left = max(0, max(cumulative_x(:))) + safety_cells;
pad_right = max(0, -min(cumulative_x(:))) + safety_cells;
pad_bottom = max(0, max(cumulative_y(:))) + safety_cells;
pad_top = max(0, -min(cumulative_y(:))) + safety_cells;

generation_grid = g;
generation_grid.x = g.x(1) + ((-pad_left):(g.sz(1) + pad_right - 1)) .* g.dx;
generation_grid.y = g.y(1) + ((-pad_bottom):(g.sz(2) + pad_top - 1))' .* g.dy;
generation_grid.z = g.z;
generation_grid.sz = [numel(generation_grid.x) numel(generation_grid.y) numel(generation_grid.z)];
generation_grid.lat = generation_grid.f_lat(generation_grid.y);
generation_grid.lon = generation_grid.f_lon(generation_grid.x);
generation_grid.mask_water = false(generation_grid.sz(1:2));
generation_grid.dem = zeros(generation_grid.sz(1:2));

crop = struct( ...
    "x_idx", (1:g.sz(1)) + pad_left, ...
    "y_idx", (1:g.sz(2)) + pad_bottom, ...
    "pad_left", pad_left, ...
    "pad_right", pad_right, ...
    "pad_bottom", pad_bottom, ...
    "pad_top", pad_top);
end

function shift_profiles = get_interval_shift_profiles(uv_profile, g, n_time, dt_seconds)

if ismatrix(uv_profile) && size(uv_profile, 2) == 2
    if size(uv_profile, 1) ~= g.sz(3)
        error("build_advection_generation_grid:InvalidProfileSize", ...
            "uv_profile must have one row per g altitude.");
    end
    shift_profiles = repmat(reshape(uv_profile, [g.sz(3) 2 1]), 1, 1, n_time - 1);
elseif ndims(uv_profile) == 3 && size(uv_profile, 1) == g.sz(3) && size(uv_profile, 2) == 2
    if size(uv_profile, 3) ~= n_time - 1
        error("build_advection_generation_grid:InvalidProfileCount", ...
            "Time-varying uv_profile must have one nz-by-2 slice per time interval.");
    end
    shift_profiles = uv_profile;
else
    error("build_advection_generation_grid:InvalidProfileShape", ...
        "uv_profile must be nz-by-2 or nz-by-2-by-(n_time-1).");
end

shift_profiles(:, 1, :) = round(shift_profiles(:, 1, :) .* dt_seconds ./ g.dx);
shift_profiles(:, 2, :) = round(shift_profiles(:, 2, :) .* dt_seconds ./ g.dy);
end
