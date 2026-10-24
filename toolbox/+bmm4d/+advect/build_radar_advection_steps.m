function steps = build_radar_advection_steps(g, radar_scans, dt_seconds, varargin)
%BUILD_RADAR_ADVECTION_STEPS Interpolate radar velocity and build steps one interval at a time.
% The global advection penalty is set later in solve_advection_lean.
% Here, shoreline_crossing_coeff only edits the local link coefficients in
% the sparse advection mapping before the inverse problem is assembled.

p = inputParser;
addParameter(p, "idw_neighbors", 10);
addParameter(p, "idw_power", 1.6);
addParameter(p, "boundary_mode", "open");
addParameter(p, "shoreline_crossing_coeff", 1);
addParameter(p, "active_mask_3d", []);
parse(p, varargin{:});

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
    interval_step = bmm4d.advect.build_profile_advection_steps( ...
        g, uv_interval, dt_seconds(i_interval), ...
        boundary_mode=p.Results.boundary_mode, ...
        shoreline_crossing_coeff=p.Results.shoreline_crossing_coeff, ...
        active_mask_3d=p.Results.active_mask_3d);
    steps{i_interval} = interval_step{1};
end
end
