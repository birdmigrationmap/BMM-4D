%% Production test 5c: two-level inversion against the production run
% Compares the two-level inversion (P4c on the uncut case
% us_east_near0_zero150) with the production lean inversion (P4 on us_east,
% near-range cut 25 km, no-echo gates within 75 km). Both share the grid
% and active cells. Scores per frame: radar-distance bias (column /
% 150 km regional mean against distance to the nearest radar, RMS log
% within 200 km), domain total, and the data misfit on the rows both runs
% use alike (30-120 km from the radar, production case). Maps of the
% column density at three frames and their ratio. Figures go to
% studies/two_level/figures/us_east.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"), fullfile(project_root, "studies", "two_level", "functions"))
c_production = prod_case("us_east");
c_two_level = prod_case("us_east_near0_zero150");
figure_dir = fullfile(project_root, "studies", "two_level", "figures", "us_east");
% Observations from the uncut case: same grid, and its rows carry range and
% elevation (the production case file predates them)
load(fullfile(c_two_level.data_dir, "case.mat"), "observations", "g", "radar_selected")
production = load(fullfile(c_production.data_dir, "inversion.mat"), "active_field", "active_id");
two_level = load(fullfile(c_two_level.data_dir, "inversion_two_level.mat"), "active_field", "active_id", "inversion_info");
assert(isequal(production.active_id, two_level.active_id), "The two runs must share the active cells.")
active_id = production.active_id;
fields = {production.active_field, two_level.active_field};
names = ["production (lean, cuts)" "two-level (no cuts)"];
times = c_production.times(:);
n_time = numel(times);

%% Column density ------------------------------------------------------------
[ix, iy, ~] = ind2sub(g.sz, active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
is_active_column = false(g.sz(1:2));
is_active_column(column_id) = true;
column = cell(1, 2);
for i_run = 1:2
    column{i_run} = nan([g.sz(1:2) n_time], "single");
    for i_t = 1:n_time
        map = reshape(accumarray(column_id, double(fields{i_run}(:, i_t)), [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
        map(~is_active_column) = nan;
        column{i_run}(:, :, i_t) = max(map, 0);
    end
end
% Domain total in millions of birds (4 km x 4 km columns)
cell_area_km2 = g.dx * g.dy / 1e6;
total = cell2mat(cellfun(@(m) squeeze(sum(m, [1 2], "omitnan")) * cell_area_km2 / 1e6, column, UniformOutput=false));

%% Radar-distance bias ------------------------------------------------------------
ratio = cell(1, 2);
bias_200km = nan(n_time, 2);
for i_run = 1:2
    [ratio{i_run}, edges_km] = radar_distance_bias(fields{i_run}, active_id, g, radar_selected);
    within = edges_km(2:end) <= 200;
    ratio_within = ratio{i_run}(within, :) ./ mean(ratio{i_run}(within, :), 1, "omitnan");
    bias_200km(:, i_run) = sqrt(mean(log(ratio_within) .^ 2, 1, "omitnan"))';
end

%% Data misfit on the shared trusted rows ------------------------------------------
% Rows 30-120 km from the radar, which the production cuts leave as they
% are (except no-echo gates beyond 75 km); no-echo gates count as 0. The two-level prediction includes
% its near-range recovery a(r), which is about 1 there except above 2 deg.
observed = full(observations.values);
observed(observed > 0 & observed < 1e-6) = 0;
trusted = full(observations.values ~= 0) & observations.row_range_m >= 30000 & observations.row_range_m <= 120000;
curve = two_level.inversion_info.near_curve;
group = discretize(double(observations.row_elev_deg), [curve.elev_from; curve.elev_to(end)]);
group(isnan(group)) = height(curve);
recovery = 1 - exp(-(double(observations.row_range_m) ./ curve.r0_m(group)) .^ curve.k(group));
F_active = observations.F(:, active_id);
clear observations
misfit = nan(n_time, 2);
for i_t = 1:n_time
    rows = trusted(:, i_t);
    for i_run = 1:2
        predicted = F_active(rows, :) * double(fields{i_run}(:, i_t));
        if i_run == 2
            predicted = recovery(rows) .* predicted;
        end
        misfit(i_t, i_run) = norm(predicted - observed(rows, i_t)) / norm(observed(rows, i_t));
    end
end
clear F_active observed
summary = table(times, total(:, 1), total(:, 2), bias_200km(:, 1), bias_200km(:, 2), misfit(:, 1), misfit(:, 2), ...
    VariableNames=["time" "total_production_M" "total_two_level_M" "bias_production" "bias_two_level" ...
    "trusted_misfit_production" "trusted_misfit_two_level"]);
disp(summary)
fprintf("Mean over frames: bias %.3f / %.3f, trusted misfit %.3f / %.3f (production / two-level)\n", ...
    mean(bias_200km, 1, "omitnan"), mean(misfit, 1, "omitnan"))
writetable(summary, fullfile(figure_dir, "two_level_comparison.csv"))

%% Plot: maps at three frames ---------------------------------------------------------
[~, frame_id] = min(abs(times - (dateshift(times(1), "start", "day") + hours([25 29 33]))), [], 1);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
c_lim = [0 quantile(column{1}(:), 0.99)];
fig = create_export_figure(3600, 3000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, 3, TileSpacing="tight", Padding="tight");
for i_f = 1:3
    i_t = frame_id(i_f);
    for i_run = 1:2
        ax = nexttile(t);
        plot_night_map(ax, g, background, column{i_run}(:, :, i_t), c_lim, states, radar_selected, c_production.figure_bounds, ones(g.sz(1:2)));
        title(ax, sprintf("%s | %s UTC | bias %.3f", names(i_run), string(times(i_t), "HH:mm"), bias_200km(i_t, i_run)), ...
            FontSize=13, FontWeight="normal")
        if i_run == 2
            colorbar(ax);
        end
    end
    ax = nexttile(t);
    both = column{1}(:, :, i_t) > 5 & column{2}(:, :, i_t) > 5;
    log_ratio = log2(column{2}(:, :, i_t) ./ column{1}(:, :, i_t));
    log_ratio(~both) = nan;
    plot_night_map(ax, g, background, log_ratio, [-1 1], states, radar_selected, c_production.figure_bounds, ones(g.sz(1:2)));
    colormap(ax, signed_map); cb = colorbar(ax); cb.Label.String = "log_2(two-level / production)";
    title(ax, "Ratio (columns > 5 birds km^{-2})", FontSize=13, FontWeight="normal")
end
export_fixed_figure(fig, fullfile(figure_dir, "two_level_maps.png"), Resolution=150)
close(fig)

%% Plot: radar-distance curves and time series ----------------------------------------
fig = create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
colors = [0.4 0.4 0.4; 0 0.45 0.74];
center = (edges_km(1:end - 1) + edges_km(2:end)) / 2;
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_run = 1:2
    plot(ax, center, mean(ratio{i_run}, 2, "omitnan"), "-o", LineWidth=2, Color=colors(i_run, :), ...
        DisplayName=sprintf("%s, bias %.3f", names(i_run), mean(bias_200km(:, i_run), "omitnan")))
end
yline(ax, 1, "k--", HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean"); ylim(ax, [0.7 1.3])
legend(ax, Location="southwest"); set(ax, FontSize=13)
title(ax, "Radar-distance curve, mean over frames", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_run = 1:2
    plot(ax, times, bias_200km(:, i_run), "-o", LineWidth=2, Color=colors(i_run, :), DisplayName=names(i_run))
end
ylabel(ax, "Radar-distance bias"); legend(ax, Location="northwest"); set(ax, FontSize=13)
title(ax, "Bias per frame", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_run = 1:2
    plot(ax, times, total(:, i_run), "-o", LineWidth=2, Color=colors(i_run, :), DisplayName=names(i_run))
end
ylabel(ax, "Birds aloft in the domain (millions)"); legend(ax, Location="northwest"); set(ax, FontSize=13)
title(ax, "Domain total", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, "two_level_bias.png"), Resolution=150)
close(fig)
