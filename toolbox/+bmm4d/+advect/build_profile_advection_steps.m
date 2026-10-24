function steps = build_profile_advection_steps(g, uv_profiles, dt_seconds, varargin)
%BUILD_PROFILE_ADVECTION_STEPS Build one advection step per time interval.
% Each step stores only valid source/destination index pairs and their
% coefficients. Inversion and synthetic generation reuse these mappings.
%
% Optional tuning:
% - boundary_mode="open" drops links advecting outside the g, while
%   "zero_gradient" clamps them to the nearest edge cell.
% - shoreline_crossing_coeff scales links whose source and destination
%   cross the land-water mask. Use 1 for full coupling, 0 to remove those
%   links, or an intermediate value to weaken cross-shore advection.

if iscell(uv_profiles)
    n_interval = numel(uv_profiles);
elseif ismatrix(uv_profiles) && size(uv_profiles, 2) == 2
    n_interval = 1;
elseif ndims(uv_profiles) == 3
    n_interval = size(uv_profiles, 3);
elseif ndims(uv_profiles) == 4 && size(uv_profiles, 4) == 2
    n_interval = 1;
elseif ndims(uv_profiles) == 5
    n_interval = size(uv_profiles, 5);
else
    error("build_profile_advection_steps:InvalidProfiles", ...
        "uv_profiles must be a cell array, a nz-by-2 profile or profile sequence, a single nx-by-ny-by-nz-by-2 field, or a nx-by-ny-by-nz-by-2-by-n_interval array.");
end

if isscalar(dt_seconds)
    dt_seconds = repmat(dt_seconds, n_interval, 1);
else
    dt_seconds = dt_seconds(:);
end

if numel(dt_seconds) ~= n_interval
    error("build_profile_advection_steps:InvalidDt", ...
        "dt_seconds must be scalar or have one value per interval.");
end

p = inputParser;
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
addParameter(p, "active_mask_3d", []);
parse(p, varargin{:});
options = p.Results;
options.boundary_mode = string(options.boundary_mode);
if ~any(options.boundary_mode == ["open" "zero_gradient"])
    error("build_profile_advection_steps:UnknownBoundaryMode", ...
        "Unknown boundary_mode '%s'. Use 'open' or 'zero_gradient'.", options.boundary_mode);
end

% Source cell indices are identical for every interval, so build them once.
nxy = g.sz(1) * g.sz(2);
[x_id, y_id] = ind2sub(g.sz(1:2), (1:nxy)');
source.x = repmat(x_id, g.sz(3), 1);
source.y = repmat(y_id, g.sz(3), 1);
source.z = repelem((1:g.sz(3))', nxy, 1);
source.id = uint32((1:prod(g.sz))');
source.horizontal_id = repmat((1:nxy)', g.sz(3), 1);

steps = cell(n_interval, 1);
for i_interval = 1:n_interval
    if iscell(uv_profiles)
        uv_interval = uv_profiles{i_interval};
    elseif ismatrix(uv_profiles)
        uv_interval = uv_profiles;
    elseif ndims(uv_profiles) == 4
        uv_interval = uv_profiles;
    elseif ndims(uv_profiles) == 5
        uv_interval = uv_profiles(:, :, :, :, i_interval);
    else
        uv_interval = uv_profiles(:, :, i_interval);
    end

    steps{i_interval} = build_interval_step( ...
        g, uv_interval, dt_seconds(i_interval), options, source);
end
end

function step = build_interval_step(g, uv_interval, dt_seconds, options, source)

boundary_mode = options.boundary_mode;
nx = g.sz(1);
ny = g.sz(2);
nxyz = prod(g.sz);
source_x = source.x;
source_y = source.y;
z_id = source.z;

if ismatrix(uv_interval) && size(uv_interval, 2) == 2
    shift_x_cells = round(uv_interval(:, 1) .* dt_seconds ./ g.dx);
    shift_y_cells = round(uv_interval(:, 2) .* dt_seconds ./ g.dy);
    dest_x = source_x + shift_x_cells(z_id);
    dest_y = source_y + shift_y_cells(z_id);
elseif ndims(uv_interval) == 4 && size(uv_interval, 1) == g.sz(1) && size(uv_interval, 2) == g.sz(2) && size(uv_interval, 3) == g.sz(3) && size(uv_interval, 4) == 2
    shift_x_cells = round(uv_interval(:, :, :, 1) .* dt_seconds ./ g.dx);
    shift_y_cells = round(uv_interval(:, :, :, 2) .* dt_seconds ./ g.dy);
    dest_x = source_x + shift_x_cells(:);
    dest_y = source_y + shift_y_cells(:);
else
    error("build_profile_advection_steps:InvalidIntervalVelocity", ...
        "Each velocity interval must be nz-by-2 or nx-by-ny-by-nz-by-2.");
end
id_edge = dest_x < 1 | dest_x > nx | dest_y < 1 | dest_y > ny;
id_nonfinite = ~isfinite(dest_x) | ~isfinite(dest_y);

if boundary_mode == "zero_gradient"
    dest_x = min(max(dest_x, 1), nx);
    dest_y = min(max(dest_y, 1), ny);
    is_valid = ~id_nonfinite;
elseif boundary_mode == "open"
    is_valid = ~(id_edge | id_nonfinite);
end

if ~isempty(options.active_mask_3d)
    active_mask_3d = logical(options.active_mask_3d);
    source_active = active_mask_3d(:);
    valid_dest_id = sub2ind(g.sz, dest_x(is_valid), dest_y(is_valid), z_id(is_valid));
    dest_active_valid = active_mask_3d(valid_dest_id);
    active_valid = false(nxyz, 1);
    active_valid(is_valid) = dest_active_valid;
    is_valid = is_valid & source_active & active_valid;
end

source_id_valid = source.id(is_valid);
dest_id_valid = uint32(sub2ind(g.sz, dest_x(is_valid), dest_y(is_valid), z_id(is_valid)));
coeff_valid = 1;
if options.shoreline_crossing_coeff ~= 1
    coeff_valid = ones(numel(source_id_valid), 1);
    dest_horizontal_id_valid = sub2ind([nx ny], dest_x(is_valid), dest_y(is_valid));
    is_shoreline_crossing = bmm4d.grid.find_land_water_transitions( ...
        g, source.horizontal_id(is_valid), dest_horizontal_id_valid);
    coeff_valid(is_shoreline_crossing) = options.shoreline_crossing_coeff;
    is_active = coeff_valid ~= 0;
    source_id_valid = source_id_valid(is_active);
    dest_id_valid = dest_id_valid(is_active);
    coeff_valid = coeff_valid(is_active);
end

step = struct( ...
    "dt_seconds", dt_seconds, ...
    "boundary_mode", boundary_mode, ...
    "shoreline_crossing_coeff", options.shoreline_crossing_coeff, ...
    "coeff_valid", coeff_valid, ...
    "source_id_valid", source_id_valid, ...
    "dest_id_valid", dest_id_valid);
end
