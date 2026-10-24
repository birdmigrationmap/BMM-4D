function [density_prior, info] = build_density_prior_fast(g, radar_scans, observations, varargin)
%BUILD_DENSITY_PRIOR_FAST Radar-profile density prior, scaled for networks.
%   One vertical profile per radar and time is fitted to the observation
%   rows (regularized least squares, smooth in height and optionally in
%   time), then interpolated horizontally (IDW by default, as in the paper).
%   Returns DENSITY_PRIOR [g.sz n_time]. Same arithmetic as the paper's
%   build_radar_density_prior (removed, MIGRATION.md step 3), but:
%   - observation rows are taken from the transposed profile operator
%     (column slicing); row slicing a sparse matrix scans all of it, which
%     made the paper version scale with radars x times x rows;
%   - with active_id, the prior is returned only on those cells, single
%     precision, [numel(active_id) n_time], built one time at a time instead
%     of as a full-grid double field.
%   - interpolation "gaussian" (length_scale_m) uses smooth Gaussian radar
%     weights (build_gaussian_weights) instead of nearest-neighbour IDW;
%   - include_water true fits the profiles on water cells too (by default
%     only the land part of each beam enters the profile fit); the profiles
%     are applied to every cell in both cases.
%   - split_land_water true fits a land and a water profile per radar and
%     time and applies each to its own cells (paper Cape May case).
%   - plot_diag true draws the interpolation diagnostics (full-grid output
%     only).
%   Observation rows are assigned to radars by observations.row_radar_id,
%   or from observations.row_ids when it is missing (single-case systems).

p = inputParser;
addParameter(p, "lambda_profile", 1);
addParameter(p, "lambda_time_profile", 0);
addParameter(p, "maxit", 100);
addParameter(p, "tol", 1e-6);
addParameter(p, "idw_neighbors", 3);
addParameter(p, "idw_power", 2.25);
addParameter(p, "active_id", []);
addParameter(p, "interpolation", "idw");
addParameter(p, "length_scale_m", 150000);
addParameter(p, "include_water", false);
addParameter(p, "split_land_water", false);
addParameter(p, "plot_diag", false);
parse(p, varargin{:});
assert(~(p.Results.include_water && p.Results.split_land_water), ...
    "build_density_prior_fast:Options", "Use either include_water or split_land_water.")

n_radar = size(radar_scans, 1);
n_time = size(radar_scans, 2);
nz = g.sz(3);
nxy = prod(g.sz(1:2));
if p.Results.interpolation == "gaussian"
    interpolation_weights = build_gaussian_weights(g, radar_scans, p.Results.length_scale_m);
else
    interpolation_weights = build_horizontal_idw_weights(g, radar_scans, ...
        n_neighbors=p.Results.idw_neighbors, power=p.Results.idw_power);
end

% Map one profile per radar-time to land cells (all cells with
% include_water; a land and a water profile with split_land_water), then to
% observation rows for the profile fit.
water_mask = logical(get_water_mask(g));
cell_id = (1:nxy * nz)';
cell_layer = ceil(cell_id ./ nxy);
is_land = ~repmat(water_mask(:), nz, 1);
if p.Results.include_water
    is_land(:) = true;
end
vertical_expansion = sparse(cell_id(is_land), cell_layer(is_land), 1, nxy * nz, nz);
n_profile_type = 1;
if p.Results.split_land_water
    vertical_expansion = [vertical_expansion sparse(cell_id(~is_land), cell_layer(~is_land), 1, nxy * nz, nz)];
    n_profile_type = 2;
end
profile_width = nz * n_profile_type;
vp_operator_t = (observations.F * vertical_expansion)';   % profile_width x n_rows, sliced by column
vp_row_active = full(any(vp_operator_t, 1))';
if isfield(observations, "row_radar_id")
    obs_radar_id = double(observations.row_radar_id);
else
    [~, ~, obs_radar_id, ~] = ind2sub(observations.obs_grid_size, double(observations.row_ids));
end

profile_smoothness_base = difference_operator(nz, 1);
profile_smoothness_base = profile_smoothness_base ./ sqrt(size(profile_smoothness_base, 1));
if p.Results.split_land_water
    profile_smoothness_base = blkdiag(profile_smoothness_base, profile_smoothness_base);
end
profile_smoothness_base_norm = norm(profile_smoothness_base, "fro");
use_time_profile = p.Results.lambda_time_profile > 0 && n_time > 1;
if use_time_profile
    temporal_difference = temporal_difference_matrix(profile_width, n_time);
    temporal_difference_norm = norm(temporal_difference, "fro");
end

% Observation rows and values per radar and time, read once per time.
row_ids = cell(n_radar, n_time);
row_values = cell(n_radar, n_time);
for i_t = 1:n_time
    [rows_t, ~, values_t] = find(observations.values(:, i_t));
    keep = vp_row_active(rows_t);
    rows_t = rows_t(keep);
    values_t = full(values_t(keep));
    radar_t = obs_radar_id(rows_t);
    for i_r = 1:n_radar
        in_radar = radar_t == i_r;
        row_ids{i_r, i_t} = rows_t(in_radar);
        row_values{i_r, i_t} = values_t(in_radar);
    end
end

% One regularized least-squares problem per radar for its profile sequence.
profiles = zeros(profile_width, n_radar, n_time);
profile_data_norms = zeros(n_radar, n_time);
for i_r = 1:n_radar
    A_parts = cell(n_time, 1);
    b_parts = cell(n_time, 1);
    for i_t = 1:n_time
        n_obs = numel(row_ids{i_r, i_t});
        data_scale = 1 / sqrt(max(n_obs, 1));
        A_block = data_scale .* vp_operator_t(:, row_ids{i_r, i_t})';
        profile_data_norms(i_r, i_t) = norm(A_block, "fro");
        smoothness_scale = regularization_scale(profile_data_norms(i_r, i_t), profile_smoothness_base_norm, p.Results.lambda_profile);
        block = [A_block; smoothness_scale .* profile_smoothness_base];
        A_parts{i_t} = [sparse(size(block, 1), (i_t - 1) * profile_width) block ...
            sparse(size(block, 1), (n_time - i_t) * profile_width)];
        b_parts{i_t} = [data_scale .* row_values{i_r, i_t}; zeros(size(profile_smoothness_base, 1), 1)];
    end
    A_time = sparse(0, profile_width * n_time);
    b_time = zeros(0, 1);
    if use_time_profile
        valid_norms = profile_data_norms(i_r, profile_data_norms(i_r, :) > 0);
        reference_norm = 0;
        if ~isempty(valid_norms)
            reference_norm = median(valid_norms);
        end
        A_time = regularization_scale(reference_norm, temporal_difference_norm, p.Results.lambda_time_profile) .* temporal_difference;
        b_time = zeros(profile_width * (n_time - 1), 1);
    end
    profile_sequence = lsqr(vertcat(A_parts{:}, A_time), vertcat(b_parts{:}, b_time), p.Results.tol, p.Results.maxit);
    profiles(:, i_r, :) = reshape(profile_sequence, profile_width, 1, n_time);
end

% Interpolate the radar profiles horizontally; with split_land_water, water
% cells take the interpolated water profile.
water_cell = repmat(water_mask(:), nz, 1);
if isempty(p.Results.active_id)
    density_prior = zeros([g.sz n_time]);
    for i_t = 1:n_time
        density_prior(:, :, :, i_t) = reshape(interpolate_profiles(interpolation_weights, profiles(:, :, i_t), ...
            nz, n_profile_type, water_cell), g.sz);
    end
else
    active_id = p.Results.active_id(:);
    density_prior = zeros(numel(active_id), n_time, "single");
    for i_t = 1:n_time
        field_t = interpolate_profiles(interpolation_weights, profiles(:, :, i_t), nz, n_profile_type, water_cell);
        density_prior(:, i_t) = single(field_t(active_id));
    end
end

if p.Results.plot_diag && isempty(p.Results.active_id)
    i_t_diag = max(1, ceil(n_time / 2));
    radar = table(cellfun(@(scan) scan.lon, radar_scans(:, 1)), cellfun(@(scan) scan.lat, radar_scans(:, 1)), ...
        VariableNames=["lon" "lat"]);
    plot_radar_interpolation_diag(density_prior(:, :, :, i_t_diag), profiles(1:nz, :, i_t_diag), radar, ...
        g=g, field_time=density_prior, time_values=1:n_time, figure_name="Density prior diagnostics", ...
        figure_title=sprintf("Density prior | t = %s", string(radar_scans{1, i_t_diag}.datetime)), ...
        xy_title_prefix="Interpolated density prior", xz_title_prefix="XZ density prior", ...
        tz_title="Time-z density prior", tx_title="Time-x density prior");
end

info = struct("radar_profiles", profiles, "profile_data_norms", profile_data_norms);
end

function field = interpolate_profiles(interpolation_weights, profiles_t, nz, n_profile_type, water_cell)
% Field [nxy*nz 1] of one time from the radar profiles [profile_width n_radar]
field = full(interpolation_weights * profiles_t(1:nz, :)');   % nxy x nz
field = field(:);
if n_profile_type == 2
    field_water = full(interpolation_weights * profiles_t(nz + 1:end, :)');
    field(water_cell) = field_water(water_cell);
end
end

function scale = regularization_scale(data_norm, regularization_norm, lambda_value)
if lambda_value == 0
    scale = 0;
    return
end
scale = lambda_value .* max(data_norm, 1) ./ max(regularization_norm, eps);
end

function D = temporal_difference_matrix(n_profile, n_time)
n_pair = n_profile * (n_time - 1);
rows = [(1:n_pair)'; (1:n_pair)'];
cols = [(1:n_pair)'; (1:n_pair)' + n_profile];
D = sparse(rows, cols, [-ones(n_pair, 1); ones(n_pair, 1)], n_pair, n_profile * n_time);
end
