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
        g, uv_interval, dt_seconds(i_interval), varargin{:});
end
end

function step = build_interval_step(g, uv_interval, dt_seconds, varargin)

p = inputParser;
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
addParameter(p, "active_mask_3d", []);
parse(p, varargin{:});
boundary_mode = string(p.Results.boundary_mode);

nx = g.sz(1);
ny = g.sz(2);
nz = g.sz(3);
nxy = nx * ny;
nxyz = prod(g.sz);

[x_id, y_id] = ind2sub([nx ny], (1:nxy)');
source_x = repmat(x_id, nz, 1);
source_y = repmat(y_id, nz, 1);
z_id = repelem((1:nz)', nxy, 1);
source_id = uint32((1:nxyz)');
source_horizontal_id = repmat((1:nxy)', nz, 1);

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
else
    error("build_profile_advection_steps:UnknownBoundaryMode", ...
        "Unknown boundary_mode '%s'. Use 'open' or 'zero_gradient'.", boundary_mode);
end

if ~isempty(p.Results.active_mask_3d)
    active_mask_3d = logical(p.Results.active_mask_3d);
    source_active = active_mask_3d(:);
    valid_dest_id = sub2ind(g.sz, dest_x(is_valid), dest_y(is_valid), z_id(is_valid));
    dest_active_valid = active_mask_3d(valid_dest_id);
    active_valid = false(nxyz, 1);
    active_valid(is_valid) = dest_active_valid;
    is_valid = is_valid & source_active & active_valid;
end

source_id_valid = source_id(is_valid);
dest_id_valid = uint32(sub2ind(g.sz, dest_x(is_valid), dest_y(is_valid), z_id(is_valid)));
coeff_valid = 1;
if p.Results.shoreline_crossing_coeff ~= 1
    coeff_valid = ones(numel(source_id_valid), 1);
    dest_horizontal_id_valid = sub2ind([nx ny], dest_x(is_valid), dest_y(is_valid));
    is_shoreline_crossing = find_land_water_transitions( ...
        g, source_horizontal_id(is_valid), dest_horizontal_id_valid);
    coeff_valid(is_shoreline_crossing) = p.Results.shoreline_crossing_coeff;
    is_active = coeff_valid ~= 0;
    source_id_valid = source_id_valid(is_active);
    dest_id_valid = dest_id_valid(is_active);
    coeff_valid = coeff_valid(is_active);
end

step = struct( ...
    "dt_seconds", dt_seconds, ...
    "boundary_mode", boundary_mode, ...
    "shoreline_crossing_coeff", p.Results.shoreline_crossing_coeff, ...
    "coeff_valid", coeff_valid, ...
    "source_id_valid", source_id_valid, ...
    "dest_id_valid", dest_id_valid);
end
