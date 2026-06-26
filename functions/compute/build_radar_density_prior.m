function [density_prior, density_prior_info] = build_radar_density_prior(grid, radar_scans, observations, varargin)
%BUILD_RADAR_DENSITY_PRIOR Build a 4D density prior from radar vertical profiles.
% Each radar/time scan gets one vertical density profile estimated from the
% compact observations. The radar profiles are then interpolated horizontally
% to the full Cartesian grid, returning DENSITY_PRIOR with size
% [grid.sz n_time].

p = inputParser;
addParameter(p, "lambda_profile", 1);
addParameter(p, "lambda_time_profile", 0);
addParameter(p, "maxit", 100);
addParameter(p, "tol", 1e-6);
addParameter(p, "idw_neighbors", 10);
addParameter(p, "idw_power", 2);
addParameter(p, "split_land_water_profiles", false);
addParameter(p, "plot_diag", false);
parse(p, varargin{:});

n_radar = size(radar_scans, 1);
n_time = size(radar_scans, 2);
nz = grid.sz(3);
nxy = prod(grid.sz(1:2));
[interpolation_weights, radar_xy] = build_horizontal_idw_weights( ...
    grid, radar_scans, ...
    n_neighbors=p.Results.idw_neighbors, ...
    power=p.Results.idw_power);

if isfield(grid, "mask_water")
    water_mask = grid.mask_water;
else
    water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
    water_mask = logical(water_mask_data.Fonwater({grid.f_lat(grid.y), grid.f_lon(grid.x)}))';
end

% Map vertical profiles to field cells and then to observation rows.
vertical_expansion_base = repelem(speye(nz), nxy, 1);
water_mask_3d = repmat(water_mask, 1, 1, nz);
if p.Results.split_land_water_profiles
    land_expansion = vertical_expansion_base;
    water_expansion = vertical_expansion_base;
    land_expansion(water_mask_3d(:), :) = 0;
    water_expansion(~water_mask_3d(:), :) = 0;
    vertical_expansion = [land_expansion water_expansion];
    n_profile_type = 2;
    profile_labels = ["land" "water"];
else
    vertical_expansion = vertical_expansion_base;
    vertical_expansion(water_mask_3d(:), :) = 0;
    n_profile_type = 1;
    profile_labels = "land";
end
profile_width = nz * n_profile_type;
vp_operator_shared = observations.F * vertical_expansion;
vp_row_active = any(vp_operator_shared, 2);
[~, ~, obs_radar_id, ~] = ind2sub( ...
    observations.obs_grid_size, double(observations.row_ids));

profile_smoothness = get_l(nz, 1);
profile_smoothness_base = profile_smoothness ./ sqrt(size(profile_smoothness, 1));
if p.Results.split_land_water_profiles
    profile_smoothness_base = blkdiag(profile_smoothness_base, profile_smoothness_base);
end
profile_smoothness_base_norm = norm(profile_smoothness_base, "fro");
nonzero_obs_indices = cell(n_time, 1);
for i_t = 1:n_time
    nonzero_obs_indices{i_t} = uint32(find(observations.values(:, i_t)));
end

density_profiles_radar = zeros(nz, n_radar, n_time, n_profile_type);
profile_data_norms = zeros(n_radar, n_time);
profile_smoothness_scales = zeros(n_radar, n_time);
time_smoothness_scales = zeros(n_radar, 1);
for i_r = 1:n_radar
    A_parts = cell(n_time, 1);
    b_parts = cell(n_time, 1);
    for i_t = 1:n_time
        nonzero_obs_idx = nonzero_obs_indices{i_t};
        profile_obs_idx = nonzero_obs_idx(vp_row_active(nonzero_obs_idx) & obs_radar_id(nonzero_obs_idx) == i_r);
        profile_obs_values = full(observations.values(profile_obs_idx, i_t));
        data_scale = 1 / sqrt(max(numel(profile_obs_values), 1));
        block_ids = (i_t - 1) * profile_width + (1:profile_width);
        n_left = block_ids(1) - 1;
        n_right = profile_width * n_time - block_ids(end);
        A_data = [ ...
            sparse(numel(profile_obs_values), n_left) ...
            data_scale .* vp_operator_shared(profile_obs_idx, :) ...
            sparse(numel(profile_obs_values), n_right)];
        profile_data_norms(i_r, i_t) = norm(A_data, "fro");
        smoothness_scale = compute_regularization_scale( ...
            profile_data_norms(i_r, i_t), ...
            profile_smoothness_base_norm, ...
            p.Results.lambda_profile);
        profile_smoothness_scales(i_r, i_t) = smoothness_scale;
        A_smooth = [ ...
            sparse(size(profile_smoothness_base, 1), n_left) ...
            smoothness_scale .* profile_smoothness_base ...
            sparse(size(profile_smoothness_base, 1), n_right)];
        A_parts{i_t} = [A_data; A_smooth];
        b_parts{i_t} = [data_scale .* profile_obs_values; sparse(size(profile_smoothness_base, 1), 1)];
    end

    A_time = sparse(0, profile_width * n_time);
    b_time_regularization = sparse(0, 1);
    if p.Results.lambda_time_profile > 0 && n_time > 1
        temporal_difference = get_temporal_difference_matrix(profile_width, n_time);
        valid_data_norms = profile_data_norms(i_r, profile_data_norms(i_r, :) > 0);
        if isempty(valid_data_norms)
            reference_data_norm = 0;
        else
            reference_data_norm = median(valid_data_norms);
        end
        temporal_scale = compute_regularization_scale( ...
            reference_data_norm, ...
            norm(temporal_difference, "fro"), ...
            p.Results.lambda_time_profile);
        time_smoothness_scales(i_r) = temporal_scale;
        row_id = zeros(2 * profile_width * (n_time - 1), 1);
        col_id = zeros(2 * profile_width * (n_time - 1), 1);
        values = zeros(2 * profile_width * (n_time - 1), 1);
        for i_interval = 1:n_time - 1
            row_block = (i_interval - 1) * profile_width + (1:profile_width)';
            ids = (i_interval - 1) * 2 * profile_width + (1:(2 * profile_width))';
            row_id(ids) = [row_block; row_block];
            col_id(ids) = [(i_interval - 1) * profile_width + (1:profile_width)'; i_interval * profile_width + (1:profile_width)'];
            values(ids) = temporal_scale .* [-ones(profile_width, 1); ones(profile_width, 1)];
        end
        A_time = sparse(row_id, col_id, values, profile_width * (n_time - 1), profile_width * n_time);
        b_time_regularization = sparse(profile_width * (n_time - 1), 1);
    end

    A_profile = vertcat(A_parts{:}, A_time);
    b_profile = vertcat(b_parts{:}, b_time_regularization);
    profile_sequence = lsqr(A_profile, b_profile, p.Results.tol, p.Results.maxit);
    profile_sequence = reshape(profile_sequence, profile_width, n_time);
    for i_profile_type = 1:n_profile_type
        profile_ids = (i_profile_type - 1) * nz + (1:nz);
        density_profiles_radar(:, i_r, :, i_profile_type) = reshape(profile_sequence(profile_ids, :), nz, 1, n_time);
    end
end

% Interpolate all altitude/time profiles in one sparse matrix multiply.
if p.Results.split_land_water_profiles
    land_profile_matrix = reshape(permute(density_profiles_radar(:, :, :, 1), [2 1 3]), n_radar, []);
    water_profile_matrix = reshape(permute(density_profiles_radar(:, :, :, 2), [2 1 3]), n_radar, []);
    density_prior_land = reshape(full(interpolation_weights * land_profile_matrix), [grid.sz n_time]);
    density_prior_water = reshape(full(interpolation_weights * water_profile_matrix), [grid.sz n_time]);
    water_mask_4d = repmat(water_mask, 1, 1, nz, n_time);
    density_prior = density_prior_land;
    density_prior(water_mask_4d) = density_prior_water(water_mask_4d);
else
    profile_matrix = reshape(permute(density_profiles_radar, [2 1 3 4]), n_radar, []);
    interpolated_profiles = full(interpolation_weights * profile_matrix);
    density_prior = reshape(interpolated_profiles, [grid.sz n_time]);
end

if p.Results.plot_diag
    i_t_diag = max(1, ceil(n_time / 2));
    radar = table( ...
        cellfun(@(scan) scan.lon, radar_scans(:, 1)), ...
        cellfun(@(scan) scan.lat, radar_scans(:, 1)), ...
        VariableNames=["lon" "lat"]);
    plot_radar_interpolation_diag( ...
        density_prior(:, :, :, i_t_diag), density_profiles_radar(:, :, i_t_diag, 1), radar, ...
        grid=grid, ...
        field_time=density_prior, ...
        time_values=1:n_time, ...
        figure_name="Density prior diagnostics", ...
        figure_title=sprintf("Density prior | t = %s", string(radar_scans{1, i_t_diag}.datetime)), ...
        xy_title_prefix="Interpolated density prior", ...
        xz_title_prefix="XZ density prior", ...
        tz_title="Time-z density prior", ...
        tx_title="Time-x density prior");
end

density_prior_info = struct( ...
    "controls", orderfields(struct(p.Results)), ...
    "setup", struct( ...
        "n_radar", n_radar, ...
        "n_time", n_time, ...
        "n_altitude", nz, ...
        "n_water_cell", nnz(water_mask), ...
        "profile_labels", profile_labels), ...
    "radar_profiles", density_profiles_radar, ...
    "profile_data_norms", profile_data_norms, ...
    "profile_smoothness_scales", profile_smoothness_scales, ...
    "time_smoothness_scales", time_smoothness_scales, ...
    "radar_xy", radar_xy);
end

function scale = compute_regularization_scale(data_norm, regularization_norm, lambda_value)

if lambda_value == 0
    scale = 0;
    return
end

reference_norm = max(data_norm, 1);
scale = lambda_value .* reference_norm ./ max(regularization_norm, eps);
end

function temporal_difference = get_temporal_difference_matrix(n_profile, n_time)

row_id = zeros(2 * n_profile * (n_time - 1), 1);
col_id = zeros(2 * n_profile * (n_time - 1), 1);
values = zeros(2 * n_profile * (n_time - 1), 1);
for i_interval = 1:n_time - 1
    row_block = (i_interval - 1) * n_profile + (1:n_profile)';
    ids = (i_interval - 1) * 2 * n_profile + (1:(2 * n_profile))';
    row_id(ids) = [row_block; row_block];
    col_id(ids) = [(i_interval - 1) * n_profile + (1:n_profile)'; i_interval * n_profile + (1:n_profile)'];
    values(ids) = [-ones(n_profile, 1); ones(n_profile, 1)];
end
temporal_difference = sparse(row_id, col_id, values, n_profile * (n_time - 1), n_profile * n_time);
end
