function [density_prior, info] = build_density_prior_fast(g, radar_scans, observations, varargin)
%BUILD_DENSITY_PRIOR_FAST Radar-profile density prior, scaled for networks.
%   Prod copy of toolbox/prior/build_radar_density_prior.m (land
%   profiles only, no diagnostics) with the same arithmetic but:
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
parse(p, varargin{:});

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
% include_water), then to observation rows for the profile fit.
water_mask = logical(get_water_mask(g));
cell_id = (1:nxy * nz)';
is_land = ~repmat(water_mask(:), nz, 1);
if p.Results.include_water
    is_land(:) = true;
end
vertical_expansion = sparse(cell_id(is_land), ceil(cell_id(is_land) ./ nxy), 1, nxy * nz, nz);
vp_operator_t = (observations.F * vertical_expansion)';   % nz x n_rows, sliced by column
vp_row_active = full(any(vp_operator_t, 1))';
obs_radar_id = double(observations.row_radar_id);

profile_smoothness_base = get_l(nz, 1);
profile_smoothness_base = profile_smoothness_base ./ sqrt(size(profile_smoothness_base, 1));
profile_smoothness_base_norm = norm(profile_smoothness_base, "fro");
use_time_profile = p.Results.lambda_time_profile > 0 && n_time > 1;
if use_time_profile
    temporal_difference = temporal_difference_matrix(nz, n_time);
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
profiles = zeros(nz, n_radar, n_time);
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
        A_parts{i_t} = [sparse(size(block, 1), (i_t - 1) * nz) block sparse(size(block, 1), (n_time - i_t) * nz)];
        b_parts{i_t} = [data_scale .* row_values{i_r, i_t}; zeros(size(profile_smoothness_base, 1), 1)];
    end
    A_time = sparse(0, nz * n_time);
    b_time = zeros(0, 1);
    if use_time_profile
        valid_norms = profile_data_norms(i_r, profile_data_norms(i_r, :) > 0);
        reference_norm = 0;
        if ~isempty(valid_norms)
            reference_norm = median(valid_norms);
        end
        A_time = regularization_scale(reference_norm, temporal_difference_norm, p.Results.lambda_time_profile) .* temporal_difference;
        b_time = zeros(nz * (n_time - 1), 1);
    end
    profile_sequence = lsqr(vertcat(A_parts{:}, A_time), vertcat(b_parts{:}, b_time), p.Results.tol, p.Results.maxit);
    profiles(:, i_r, :) = reshape(profile_sequence, nz, 1, n_time);
end

% Interpolate the radar profiles horizontally.
if isempty(p.Results.active_id)
    profile_matrix = reshape(permute(profiles, [2 1 3]), n_radar, []);
    density_prior = reshape(full(interpolation_weights * profile_matrix), [g.sz n_time]);
else
    active_id = p.Results.active_id(:);
    density_prior = zeros(numel(active_id), n_time, "single");
    for i_t = 1:n_time
        field_t = full(interpolation_weights * profiles(:, :, i_t)');   % nxy x nz
        density_prior(:, i_t) = single(field_t(active_id));
    end
end

info = struct("radar_profiles", profiles, "profile_data_norms", profile_data_norms);
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
