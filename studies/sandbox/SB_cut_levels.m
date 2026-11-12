%% Sandbox: levels per radar, cut type and distance in the 3D inversion
% Dev study, the 3D version of gate_density GD_M. Every observation is read as
%   y = c(radar, cut type, distance) * beta(lambda) * (F rho)
% - c: one level per radar, cut type (split: surveillance pass of a cut
%   scanned twice; batch: one pass) and distance knot (6, 15, 25, 37.5,
%   62.5 km, linear in between, constant beyond), within [P_out / P_all, 1.1],
%   the share of the echo the radar's clutter filter left on that radar's cuts
%   of that type and band over the night (gd_cfp_power). Split cuts at the
%   last knot are the anchor (1). Each level is pulled toward 1 by one
%   sweep's worth of pseudo-observations at the prediction.
% - beta: the sampling curve by expected birds per gate lambda = (F rho) x
%   counting volume, learned on the 124 reference scans (GD_M), fixed here.
% Solved in rounds: the field with y / f and weights w f^2 (f = c beta, w the
% production weights from the estimator's sd at the prior, floor 1.5), the
% radar prior rebuilt from y / f; then the levels from the new field.
% CASE_NAME (default "sandbox", kept if set before running, e.g.
% "night_20261021"); a case cleaned without the estimator's sd uses equal
% weights instead of the production weights.
% Outputs: figures/cut_levels<suffix>.csv (levels and bounds per radar),
% cut_levels_sections<suffix>_<HHmm>.png, cut_levels_sections<suffix>.gif,
% cut_levels<suffix>.mat (fields of the first and last round); suffix "" for
% the sandbox, "_<case>" otherwise.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"), ...
    fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("case_name", "var")
    case_name = "sandbox";
end
suffix = "";
if case_name ~= "sandbox"
    suffix = "_" + case_name;
end
S = sandbox_load(case_name);
c = case_settings(case_name);
manifest = load(fullfile(c.data_dir, "manifest.mat"));
gd = load(fullfile(project_root, "studies", "gate_density", "figures", "GD_M_cut_levels.mat"), "beta", "lambda_knots");
beam = bmm4d.radar.gate.beam();

%% Settings ---------------------------------------------------------------
band_edges_km = [2 10 20 30 45 80];
knots_km = (band_edges_km(1:end - 1) + band_edges_km(2:end)) / 2;
n_knot = numel(knots_km);
upper_level = 1.1;
n_round = 3;
n_radar = height(S.radar_selected);

%% 1. Reported and removed power per radar, cut type and band, over the night ---------------
power_file = fullfile(figure_dir, "cut_levels_power" + suffix + ".mat");
if isfile(power_file)
    load(power_file, "power")
else
    raw_files = manifest.raw_files(manifest.raw_files ~= "");
    power = cell(numel(raw_files), 1);
    pool = gcp("nocreate");
    if isempty(pool) || pool.NumWorkers > 4
        delete(pool);
        parpool("Processes", 4);
    end
    parfor k = 1:numel(raw_files)
        power{k} = gd_cfp_power(raw_files(k), band_edges_km, 6);
    end
    [~, stem] = fileparts(raw_files);
    station = extractBefore(string(stem), 5);
    save(power_file, "power", "station")
end
load(power_file, "station")
delete(gcp("nocreate"))
% Bounds per radar, type and band; the elevations of each radar's split cuts
bounds = ones(n_knot, 2, n_radar, 2);
split_elevations = cell(n_radar, 1);
for i_r = 1:n_radar
    volumes = power(station == S.radar_selected.name(i_r));
    p_out = zeros(n_knot, 2);
    p_all = zeros(n_knot, 2);
    elevations = [];
    for k = 1:numel(volumes)
        for g = 1:2
            of_type = volumes{k}.is_split == (g == 1);
            p_out(:, g) = p_out(:, g) + sum(volumes{k}.p_out(:, of_type), 2);
            p_all(:, g) = p_all(:, g) + sum(volumes{k}.p_all(:, of_type), 2);
        end
        elevations = [elevations; volumes{k}.elev(volumes{k}.is_split)]; %#ok<AGROW>
    end
    share = p_out ./ p_all;
    share(p_all == 0) = 1;
    bounds(:, :, i_r, 1) = share;
    bounds(:, :, i_r, 2) = upper_level;
    bounds(n_knot, 1, i_r, :) = 1;                   % anchor: split cuts at the last knot
    split_elevations{i_r} = unique(round(elevations, 1));
end

%% 2. Observations: radar, cut type, distance weights ------------------------------------
O = S.observations;
[row, col, y] = find(O.values);
y(y < 1e-6) = 0;                                     % zero density is stored as eps
radar_of_row = double(O.row_radar_id);
range_km = double(O.row_range_m) / 1000;
elevation = round(double(O.row_elev_deg), 1);
% A row is split if its elevation is one of that radar's split cuts
is_split_row = false(size(radar_of_row));
for i_r = 1:n_radar
    rows = radar_of_row == i_r;
    is_split_row(rows) = ismember(elevation(rows), split_elevations{i_r}) | ...
        ismember(elevation(rows) + 0.1, split_elevations{i_r}) | ismember(elevation(rows) - 0.1, split_elevations{i_r});
end
type_of_row = 2 - is_split_row;
H = sparse(interp1(knots_km(:), eye(n_knot), min(max(range_km, knots_km(1)), knots_km(end)), "linear"));
volume = beam.volume_sample_km3(double(O.row_range_m));
hat_lambda = @(lambda) interp1(log(gd.lambda_knots(:)), eye(numel(gd.lambda_knots)), ...
    min(max(log(lambda), log(gd.lambda_knots(1))), log(gd.lambda_knots(end))), "linear");
% Distinct elevations per radar, type and knot (one sweep's worth of pull)
sweeps_at_knot = zeros(n_knot, 2, n_radar);
for i_r = 1:n_radar
    for g = 1:2
        for k = 1:n_knot
            rows = radar_of_row == i_r & type_of_row == g & H(:, k) > 0;
            sweeps_at_knot(k, g, i_r) = numel(unique(elevation(rows)));
        end
    end
end

%% 3. Rounds: field, then levels ----------------------------------------------------------
prior_prediction = O.F(:, S.active_id) * double(sandbox_radar_prior(S));
if isfield(O, "sd")
    weights_production = bmm4d.solve.observation_weights(O, prior_prediction, c.inversion.observation_error_floor);
    w = weights_production(sub2ind(size(weights_production), row, col));
else
    w = ones(size(y));                               % volumes cleaned without the estimator's sd
end
level = ones(n_knot, 2, n_radar);
factor = ones(size(y));
runs = cell(n_round + 1, 1);
x_initial = [];
for i_round = 0:n_round
    % Field from the corrected observations, with the prior rebuilt from them
    S_round = S;
    S_round.observations.values = sparse(row, col, max(y ./ factor, 1e-9), size(O.values, 1), size(O.values, 2));
    S_round.prior_cache = containers.Map();
    weights = sparse(row, col, w .* factor .^ 2, size(O.values, 1), size(O.values, 2));
    runs{i_round + 1} = sandbox_run(S_round, struct("name", sprintf("cut levels, round %d", i_round), "video", false, ...
        "row_weight", @(r, e, z) weights, "x_initial", x_initial));
    x_initial = runs{i_round + 1}.field;
    fprintf("round %d: %.1f s, %d iterations\n", i_round, runs{i_round + 1}.seconds, runs{i_round + 1}.iterations)
    if i_round == n_round, break, end
    % Levels from this field: multiplicative update per radar, type and knot, pulled toward 1
    predicted = O.F(:, S.active_id) * double(runs{i_round + 1}.field);
    predicted = max(predicted(sub2ind(size(predicted), row, col)), 1e-6);
    sampling = hat_lambda(predicted .* volume(row)) * gd.beta;
    base = sampling .* predicted;
    group_of_row = sub2ind([2 n_radar], type_of_row, radar_of_row);
    level_at = @(level) level_of_rows(H, level, group_of_row);
    for inner = 1:20
        level_obs = level_at(level);
        ratio = y ./ max(level_obs(row) .* base, 1e-9);
        % Sums over times per row, then over the rows of each radar, type and knot
        update_row = accumarray(row, w .* base .* ratio, [numel(radar_of_row) 1]);
        weight_row = accumarray(row, w .* base, [numel(radar_of_row) 1]);
        update_sum = zeros(n_knot, 2 * n_radar);
        weight_sum = zeros(n_knot, 2 * n_radar);
        for k = 1:n_knot
            update_sum(k, :) = accumarray(group_of_row, full(H(:, k)) .* update_row, [2 * n_radar 1])';
            weight_sum(k, :) = accumarray(group_of_row, full(H(:, k)) .* weight_row, [2 * n_radar 1])';
        end
        flat_level = reshape(level, n_knot, []);
        pseudo = weight_sum ./ max(reshape(sweeps_at_knot, n_knot, []), 1);
        new_level = (flat_level .* update_sum + pseudo) ./ max(weight_sum + pseudo, eps);
        new_level(weight_sum == 0) = 1;
        level = min(max(reshape(new_level, n_knot, 2, n_radar), bounds(:, :, :, 1)), bounds(:, :, :, 2));
    end
    level_obs = level_at(level);
    factor = level_obs(row) .* sampling;
    fprintf("  batch levels at 6/15/25 km, median over radars: %s\n", sprintf("%5.2f", median(squeeze(level(1:3, 2, :)), 2)))
end

%% 4. Levels and bounds per radar ------------------------------------------------------------
batch_loss = 1 - squeeze(mean(bounds(2:3, 2, :, 1), 1));
level_table = table(S.radar_selected.name, batch_loss, squeeze(level(1, 2, :)), squeeze(level(2, 2, :)), squeeze(level(3, 2, :)), ...
    squeeze(bounds(1, 2, :, 1)), squeeze(bounds(2, 2, :, 1)), squeeze(level(1, 1, :)), squeeze(level(2, 1, :)), ...
    VariableNames=["radar" "batch_filter_loss" "batch_6km" "batch_15km" "batch_25km" "bound_6km" "bound_15km" "split_6km" "split_15km"]);
level_table = sortrows(level_table, "batch_filter_loss");
writetable(level_table, fullfile(figure_dir, "cut_levels" + suffix + ".csv"))
first_field = runs{1}.field;
last_field = runs{end}.field;
save(fullfile(bmm4d.util.local_data_dir(), "prod", case_name, "cut_levels.mat"), "first_field", "last_field", "level", "bounds", "-v7.3")
disp(level_table)
for k = [1 numel(runs)]
    r = runs{k};
    within = r.edges_km(2:end) <= 100;
    fprintf("%s: column / regional mean by radar distance %s, bias %.3f, detail %.3f, trusted misfit %.3f\n", r.name, ...
        sprintf("%5.2f", mean(r.ratio(within, :), 2, "omitnan")), r.bias_200km, r.detail, r.misfit_trusted)
end

%% 5. Sections through six radars, every frame ------------------------------------------------
% From low to high batch filter loss: KMPX, KLSX, KDVN, KDMX, KVNX, KDDC; rows: production
% weights (round 0), cut levels (last round), log2 ratio
g = S.g;
radar_selected = S.radar_selected;
section_radars = ["KMPX" "KLSX" "KDVN" "KDMX" "KVNX" "KDDC"];
row_names = ["round 0 (no levels)" "cut levels" "log_2 ratio"];
effective_radius = 6371200 * 4 / 3;
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
fields = {runs{1}.field, runs{end}.field};
n_time = size(fields{1}, 2);
gif_file = fullfile(figure_dir, "cut_levels_sections" + suffix + ".gif");
fig = bmm4d.plot.create_export_figure(4800, 2300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, 3, numel(section_radars), TileSpacing="tight", Padding="tight");
    title(t, case_name + ", " + string(c.times(i_t), "dd MMM yyyy HH:mm") + " UTC, west-east sections (dashed: 0.5° beam centre)", FontSize=18, FontWeight="bold")
    volumes = cell(1, 2);
    for k = 1:2
        volumes{k} = nan(g.sz, "single");
        volumes{k}(S.active_id) = fields{k}(:, i_t);
    end
    for i_row = 1:3
        for i_k = 1:numel(section_radars)
            i_r = find(radar_selected.name == section_radars(i_k), 1);
            [~, y_id] = min(abs(g.lat - radar_selected.lat(i_r)));
            if i_row == 3
                field = log2(max(volumes{2}, 1) ./ max(volumes{1}, 1));
            else
                field = volumes{i_row};
            end
            values = squeeze(field(:, y_id, :))';
            distance_km = (g.lon - radar_selected.lon(i_r)) * 111.32 * cosd(radar_selected.lat(i_r));
            ax = nexttile(t);
            imagesc(ax, distance_km, squeeze(g.z), values, AlphaData=~isnan(values))
            hold(ax, "on")
            beam_range = linspace(0, 200e3, 200);
            beam_height = sqrt(beam_range .^ 2 + effective_radius ^ 2 + 2 * beam_range * effective_radius * sind(0.5)) - effective_radius + radar_selected.height(i_r);
            plot(ax, beam_range / 1000, beam_height, "w--", -beam_range / 1000, beam_height, "w--", LineWidth=0.8)
            set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12, XLim=[-150 150], YLim=[0 3000])
            if i_row == 3
                clim(ax, [-1 1]); colormap(ax, signed_map)
                xlabel(ax, "Distance (km)")
            else
                clim(ax, [0 prctile(volumes{1}(:, y_id, :), 99, "all")]); colormap(ax, parula)
            end
            i_level = find(level_table.radar == section_radars(i_k));
            if i_row == 1
                title(ax, sprintf("%s, filter loss %.2f", section_radars(i_k), level_table.batch_filter_loss(i_level)), FontSize=14)
            end
            if i_row == 2
                text(ax, 0.02, 0.95, sprintf("batch c %.2f / %.2f / %.2f", level_table.batch_6km(i_level), ...
                    level_table.batch_15km(i_level), level_table.batch_25km(i_level)), Units="normalized", Color="w", FontSize=11)
            end
            if i_k == 1
                ylabel(ax, ["\bf" + row_names(i_row) "\rmHeight (m a.s.l.)"])
            end
            if i_k == numel(section_radars) && i_row >= 2
                cb = colorbar(ax);
                cb.FontSize = 11;
                if i_row == 3, cb.Label.String = "log_2 ratio"; else, cb.Label.String = "birds km^{-3}"; end
            end
        end
    end
    drawnow
    frame_image = getframe(fig);
    [indexed, color_map] = rgb2ind(frame_image.cdata, 256);
    if i_t == 1
        imwrite(indexed, color_map, gif_file, LoopCount=Inf, DelayTime=1.2)
    else
        imwrite(indexed, color_map, gif_file, WriteMode="append", DelayTime=1.2)
    end
    if any(i_t == round(n_time * [0.5 1]))
        bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "cut_levels_sections" + suffix + "_" + string(c.times(i_t), "HHmm") + ".png"), Resolution=150)
    end
end
close(fig)

%% Local functions ----------------------------------------------------------
function prior_active = sandbox_radar_prior(S)
% The radar-profile prior of the production inversion
inv = S.c.inversion;
evalc("prior_active = bmm4d.prior.build_density_prior_fast(S.g, S.radar_scans, S.observations, " + ...
    """lambda_profile"", inv.lambda_profile, ""lambda_time_profile"", inv.lambda_time_profile, " + ...
    """maxit"", inv.maxit, ""active_id"", S.active_id, ""interpolation"", inv.prior_interpolation, " + ...
    """length_scale_m"", inv.prior_length_scale_m, ""include_water"", true);");
end

function level_row = level_of_rows(H, level, group_of_row)
% Level of each row: its radar and type's knot levels, linear in distance
flat_level = reshape(level, size(level, 1), []);
level_row = full(sum(H .* flat_level(:, group_of_row)', 2));
end
