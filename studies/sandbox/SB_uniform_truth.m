%% Sandbox: a uniform bird density through the whole chain
% Truth: rho_true = 50 birds km^-3 everywhere from the ground to the grid
% top (3000 m a.s.l.), on the sandbox geometry (28 radars, 6 frames).
% Two simulations, both inverted like the sandbox (sandbox_run):
% 1. Perfect data: d = F rho_true on the sandbox rows. Tests the inversion.
% 2. Gates: every gate of every sweep gets the summed echo of a Poisson
%    number of birds (lambda = rho V_sample(r), from the bird-sum table;
%    gates drawn independently), rounded to 0.5 dB, with gates below a
%    detection limit set to no echo. Then the production estimator
%    (estimate_density), the gates kept by the production cleaning
%    (same mask as the real cleaned volume), build_observations and the
%    inversion. Tests estimator + rows + inversion.
% The gate density is rho_true times the share of the beam (two-way
% Gaussian, 4/3 Earth, as compute_F_elev) between the ground and 3000 m.
% Scores against the truth, by distance to the nearest radar in 1 km bins:
% reconstructed / true column, and the density by height.
% Outputs: studies/sandbox/figures/uniform_truth_*.png, .csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"), ...
    fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
if ~exist("S", "var")
    S = sandbox_load("sandbox");
end
g = S.g;
rho_true = 50;                                   % birds km^-3
top_m = 3000;                                    % grid top (m a.s.l.)
truth = rho_true * ones(numel(S.active_id), 1);
n_time = size(S.observations.values, 2);

%% 1. Perfect data -----------------------------------------------------------
S_perfect = S;
S_perfect.prior_cache = containers.Map();
[r, t, ~] = find(S.observations.values);
d_true = S.observations.F(:, S.active_id) * truth;
S_perfect.observations.values = sparse(r, t, d_true(r), size(S.observations.values, 1), n_time);
result_perfect = sandbox_run(S_perfect, struct("name", "uniform, perfect data", "video", false));

%% 2. Simulated gates through the estimator -----------------------------------
c = case_settings("sandbox");
synthetic_dir = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes", "uniform_" + rho_true);
c_syn = c;
c_syn.name = "sandbox_uniform";
c_syn.data_dir = fullfile(bmm4d.util.local_data_dir(), "prod", c_syn.name);
manifest = load(fullfile(c.data_dir, "manifest.mat"));
[~, stem] = fileparts(manifest.volume_files);
volume_files = manifest.volume_files;
has = volume_files ~= "";
volume_files(has) = fullfile(synthetic_dir, stem(has) + ".mat");
todo = find(has & ~isfile(volume_files));
if ~isempty(todo)
    if ~isfolder(synthetic_dir)
        mkdir(synthetic_dir);
    end
    beam = bmm4d.radar.gate.beam();
    table_s = bmm4d.radar.gate.bird_sum_table(beam, 2000, 20000);
    parms = struct("minRange", 3000, "maxRange", 150000, "maxAlt", 3500, "noEchoRange", 150000, "smoothing", 2000);
    source_files = manifest.volume_files(todo);
    out_files = volume_files(todo);
    parfor k = 1:numel(todo)
        loaded = load(source_files(k), "volume");
        volume = loaded.volume;
        volume.dens.sweeps = simulate_volume(volume, rho_true, top_m, beam, table_s, parms);
        parsave(out_files(k), volume)
    end
end
if ~isfile(fullfile(c_syn.data_dir, "case.mat"))
    if ~isfolder(c_syn.data_dir)
        mkdir(c_syn.data_dir);
    end
    raw_files = manifest.raw_files; radars = manifest.radars; times = manifest.times;
    save(fullfile(c_syn.data_dir, "manifest.mat"), "raw_files", "radars", "times", "volume_files")
    copyfile(fullfile(c.data_dir, "static_clutter.mat"), c_syn.data_dir)
    bmm4d.build_observations(c_syn);
end
loaded = load(fullfile(c_syn.data_dir, "case.mat"), "observations");
assert(isequal(size(loaded.observations.values), size(S.observations.values)), "Synthetic rows differ from the sandbox rows.")
S_gates = S;
S_gates.prior_cache = containers.Map();
S_gates.observations = loaded.observations;
result_gates = sandbox_run(S_gates, struct("name", "uniform, gates + estimator", "video", false));

%% Score against the truth ------------------------------------------------------
results = {result_perfect, result_gates};
[ix, iy, iz] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for i_r = 1:height(S.radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(i_r)) .* cosd(lon_grid - S.radar_selected.lon(i_r)))));
end
edges_km = 0:1:250;
center_km = edges_km(1:end - 1)' + 0.5;
column_bin = discretize(distance_km(:), edges_km);
true_column = accumarray(column_id, truth, [prod(g.sz(1:2)) 1]);
is_active = true_column > 0;
z_layer = squeeze(g.z);
height_edges = 0:250:3000;
column_ratio = nan(numel(center_km), numel(results));
height_ratio = nan(numel(center_km), numel(height_edges) - 1, numel(results));
cell_bin = column_bin(column_id);
height_bin = discretize(z_layer(iz), height_edges);
for i_run = 1:numel(results)
    field = double(results{i_run}.field);
    column = zeros(prod(g.sz(1:2)), n_time);
    for i_t = 1:n_time
        column(:, i_t) = accumarray(column_id, field(:, i_t), [prod(g.sz(1:2)) 1]);
    end
    ok = is_active & ~isnan(column_bin);
    column_ratio(:, i_run) = accumarray(column_bin(ok), mean(column(ok, :), 2) ./ true_column(ok), [numel(center_km) 1], @mean, nan);
    ok = ~isnan(cell_bin) & ~isnan(height_bin);
    height_ratio(:, :, i_run) = accumarray([cell_bin(ok) height_bin(ok)], mean(field(ok, :), 2) / rho_true, ...
        [numel(center_km) numel(height_edges) - 1], @mean, nan);
end
writetable(table(center_km, column_ratio(:, 1), column_ratio(:, 2), VariableNames=["distance_km" "perfect_data" "gates_estimator"]), ...
    fullfile(figure_dir, "uniform_truth.csv"))

%% Plot results -----------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
title(t, sprintf("Uniform truth (%d birds km^{-3}): reconstructed / true", rho_true), FontSize=14)
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
for i_run = 1:numel(results)
    field = double(results{i_run}.field);
    column = accumarray(column_id, mean(field, 2), [prod(g.sz(1:2)) 1]);
    ratio_map = nan(g.sz(1:2));
    ratio_map(is_active) = column(is_active) ./ true_column(is_active);
    ax = nexttile(t, (i_run - 1) * 3 + 1);
    bmm4d.plot.plot_night_map(ax, g, background, log2(ratio_map), [-0.5 0.5], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, signed_map); cb = colorbar(ax); cb.Label.String = "log_2(column / true), mean over frames";
    title(ax, results{i_run}.name, FontWeight="bold")
    ax = nexttile(t, (i_run - 1) * 3 + 2); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, center_km, column_ratio(:, i_run), "k-", LineWidth=1.5)
    yline(ax, 1, "k--"); xlim(ax, [0 250]); ylim(ax, [0.6 1.4]); set(ax, FontSize=12)
    xlabel(ax, "Distance to the nearest radar (km), 1 km bins"); ylabel(ax, "Column / true column")
    ax = nexttile(t, (i_run - 1) * 3 + 3);
    imagesc(ax, center_km, height_edges(1:end - 1) + 125, height_ratio(:, :, i_run)', AlphaData=~isnan(height_ratio(:, :, i_run)'))
    set(ax, YDir="normal", Color="w", FontSize=12); clim(ax, [0.5 1.5]); colormap(ax, signed_map)
    cb = colorbar(ax); cb.Label.String = "Density / true";
    xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Height (m a.s.l.)")
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "uniform_truth.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function sweeps = simulate_volume(volume, rho_true, top_m, beam, table_s, parms)
% Simulated gates of every sweep, through the production estimator, kept
% where the real cleaned sweep has a value (same rows as the sandbox).
sweeps = volume.dens.sweeps;
earth_radius_m = 4 / 3 * 6371200;
sigma_theta = deg2rad(0.925 / 2 / sqrt(2 * log(2))) / sqrt(2);   % compute_F_elev: exp(-dtheta^2 / sigma_tp^2)
ground_m = volume.height - 20;                                     % antenna about 20 m above ground
log10_lambda = table_s.log10_lambda(:);
for i_e = 1:numel(sweeps)
    keep = ~isnan(sweeps(i_e).data);
    if ~any(keep(:))
        continue
    end
    n_range = find(any(keep, 2), 1, "last") + 8;
    range_m = double(sweeps(i_e).range(1:n_range));
    range_m = range_m(:);
    n_az = numel(sweeps(i_e).az);
    % Beam share between the ground and the grid top
    centre_m = volume.height + range_m * sind(sweeps(i_e).elev) + range_m .^ 2 / (2 * earth_radius_m);
    spread_m = range_m * sigma_theta;
    fill = normcdf((top_m - centre_m) ./ spread_m) - normcdf((ground_m - centre_m) ./ spread_m);
    lambda = rho_true * fill .* beam.volume_sample_km3(range_m);
    % Draw each gate's summed echo from the table column of its lambda
    k = min(max(round((log10(lambda) - log10_lambda(1)) / 0.005) + 1, 1), numel(log10_lambda));
    s = zeros(n_range, n_az);
    for i_r = 1:n_range
        if lambda(i_r) < 10 ^ log10_lambda(1)
            continue
        end
        u = rand(1, n_az);
        id = sum(table_s.cdf(:, k(i_r)) < u, 1) + 1;
        id = min(id, numel(table_s.log10_s));
        s(i_r, :) = 10 .^ (table_s.log10_s(id)' - 0.02 * rand(1, n_az));
        s(i_r, id == 1) = 0;
    end
    % Radar output: 0.5 dB steps, no echo below the clear-air detection limit
    detection_dbz = -28 + 20 * log10(range_m / 5000);
    eta = s ./ beam.volume_eff_km3(range_m);
    dbz = gd_observe(eta, detection_dbz, -Inf);
    gate = struct("range", range_m, "az", sweeps(i_e).az, "elev", sweeps(i_e).elev, "data", dbz, ...
        "no_echo", isnan(dbz), "unit", "dBZ");
    scan = bmm4d.radar.filter.estimate_density(struct("gates", gate), parms);
    data = nan(size(sweeps(i_e).data), "like", sweeps(i_e).data);
    data(1:n_range, :) = scan.gates.data;
    data(~keep) = nan;
    sweeps(i_e).data = data;
end
end

function parsave(file, volume)
save(file, "volume")
end
