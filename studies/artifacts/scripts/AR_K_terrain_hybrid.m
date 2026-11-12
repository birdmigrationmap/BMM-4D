%% Artefacts K: how should the reconstruction move birds over terrain? (dev)
% Follow-up of fix 2 / 2b. Three questions:
% 1. What terrain does the grid use? The real ground (30" DEM) against the
%    4 km solve-grid ground (nearest sampling at 2 km, block mean to 4 km)
%    and the cells switched off under it; the paper code did not use the
%    DEM at all.
% 2. What do the radar data show over terrain? Each observation row's beam
%    centre (F-weighted height) and the ground under its footprint; density
%    by height above ground and above sea level over low and high ground;
%    rows of the same beam (radar, frame, sweep, range) over different
%    ground; which height coordinate makes the rows of one radar collapse
%    onto one profile: height above ground (AGL), above sea level (ASL),
%    the hybrid "follow near the ground, level aloft" (fade height H), or
%    height above the ground smoothed over L km. Which rows to trust over
%    the ridges (blockage, lowest height seen).
% 3. The simplest advection that fits: ar_steps / ar_profile_prior
%    "smooth" (height above the ground smoothed over L) and "hybrid" (fade
%    height H), with or without detour (links into terrain sent sideways),
%    against ASL and AGL, on the solve and on held-out mountain radars
%    (prior rebuilt without the radar), in a batch job (ar_terrain_solves).
% Case ar_appalachians_bb (16 radars, 04:30-07:30 UTC, beam-blockage fix),
% production inversion settings with the white prior.
% Outputs: figures/AR_K_*.png, AR_K_*.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_appalachians_bb", [-89 -78 33.5 40], times);
S = ar_load(c);
g = S.g;
inv = c.inversion;
dt_seconds = seconds(diff(c.times));
results_file = fullfile(c.data_dir, "AR_K_results.mat");
section_time = datetime("2026-10-22 06:00");

%% What terrain the grid uses (figure K1) -------------------------------------------
% The DEM (data/reference/DEM/dem.mat) is GTOPO30 at 30" (0.93 x 0.75 km
% at 36 N), not smoothed. create_grid samples it at the 2 km grid cell
% centres (nearest, the Gaussian smoothing line is commented out),
% coarsen_grid averages 2 x 2 cells to the 4 km solve grid (g.dem), and
% build_active_mask switches off the cells whose top is below g.dem
% (air_fraction gives the lowest cell its air share). Heights are above
% sea level; there is no smoothed topography and no terrain-following
% coordinate. The paper code (BMM-within) loaded the DEM in create_grid but
% never used g.dem: its grid started at sea level everywhere.
dem_data = load(bmm4d.util.data_path("reference", "DEM", "dem.mat"), "dem", "dem_lat", "dem_lon");
dem_full = double(dem_data.dem);
dem_full(dem_full == -9999) = 0;                      % ocean
ground = max(g.dem, 0);
smooth_ground = imgaussfilt(ground, 100e3 / g.dx, Padding="replicate");
z = squeeze(g.z);
grid_top = z(end) + g.dz / 2;
terrain_map = interp1([0 0.45 1], [0.82 0.91 0.76; 0.86 0.77 0.58; 0.48 0.36 0.25], linspace(0, 1, 256));
map_box = [-85.5 -81.0 34.4 37.0];
section_lat = 35.6;                                   % through the Smokies (Clingmans Dome 35.56 N)
font = 13;

% How far the 4 km ground is from the real one, per 4 km cell
[lon_cells, lat_cells] = ndgrid(g.lon, g.lat);
in_box = lon_cells >= map_box(1) & lon_cells <= map_box(2) & lat_cells >= map_box(3) & lat_cells <= map_box(4);
lat_id = find(dem_data.dem_lat >= map_box(3) & dem_data.dem_lat <= map_box(4));
lon_id = find(dem_data.dem_lon >= map_box(1) & dem_data.dem_lon <= map_box(2));
cell_x = discretize(g.f_X(dem_data.dem_lon(lon_id)), [g.x(:) - g.dx / 2; g.x(end) + g.dx / 2]);
cell_y = discretize(g.f_Y(dem_data.dem_lat(lat_id)), [g.y(:) - g.dy / 2; g.y(end) + g.dy / 2]);
[cy_, cx_] = ndgrid(cell_y, cell_x);
pixel = dem_full(lat_id, lon_id);
cell_of_pixel = sub2ind(g.sz(1:2), cx_(:), cy_(:));
cell_max = accumarray(cell_of_pixel, pixel(:), [prod(g.sz(1:2)) 1], @max, NaN);
cell_mean = accumarray(cell_of_pixel, pixel(:), [prod(g.sz(1:2)) 1], @mean, NaN);
fprintf("4 km cells in the box: grid ground - cell mean %.0f m (median |.|), cell max - grid ground median %.0f m, max %.0f m\n", ...
    median(abs(ground(in_box) - cell_mean(in_box)), "omitnan"), median(cell_max(in_box) - ground(in_box), "omitnan"), ...
    max(cell_max(in_box) - ground(in_box)))

fig = bmm4d.plot.create_export_figure(3600, 2300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
maps = {dem_data.dem_lon(lon_id), dem_data.dem_lat(lat_id), pixel, "Real ground (30"" DEM)"; ...
    g.lon, g.lat, ground', "Grid ground (4 km, g.dem)"; ...
    g.lon, g.lat, smooth_ground', "Ground smoothed over 100 km"};   % [y x] for imagesc
for i_m = 1:3
    ax = nexttile(layout);
    imagesc(ax, maps{i_m, 1}, maps{i_m, 2}, maps{i_m, 3});
    hold(ax, "on")
    plot(ax, map_box(1:2), [section_lat section_lat], "k--", LineWidth=1.2)
    radar_in = S.radar_selected.lon >= map_box(1) & S.radar_selected.lon <= map_box(2) & ...
        S.radar_selected.lat >= map_box(3) & S.radar_selected.lat <= map_box(4);
    plot(ax, S.radar_selected.lon(radar_in), S.radar_selected.lat(radar_in), "k^", MarkerFaceColor="w", MarkerSize=8)
    text(ax, S.radar_selected.lon(radar_in) + 0.06, S.radar_selected.lat(radar_in) + 0.06, string(S.radar_selected.name(radar_in)), ...
        FontSize=font - 2, FontName="Helvetica")
    set(ax, YDir="normal", FontSize=font, DataAspectRatio=[1 cosd(35.7) 1])
    axis(ax, map_box)
    colormap(ax, terrain_map)
    clim(ax, [0 2000])
    title(ax, maps{i_m, 4}, FontSize=font + 1, FontWeight="bold")
    xlabel(ax, "Longitude (°)")
    if i_m == 1
        ylabel(ax, "Latitude (°)")
    end
end
cb = colorbar(ax, Color="k");
cb.Label.String = "Ground (m asl)";

% West-east section through the Smokies: real ground, grid ground, active
% cells shaded by their air share (switched-off cells white), smoothed ground
[~, y_id] = min(abs(g.lat - section_lat));
x_id = find(g.lon >= map_box(1) & g.lon <= map_box(2));
[~, dem_row] = min(abs(dem_data.dem_lat - section_lat));
air = zeros(g.sz);
air(S.active_id) = bmm4d.grid.air_fraction(g, S.active_id);
ax = nexttile(layout, [1 2]);
imagesc(ax, g.lon(x_id), z / 1000, squeeze(air(x_id, y_id, :))');   % imagesc wants [z x]
hold(ax, "on")
colormap(ax, interp1([0 1], [1 1 1; 0.70 0.82 0.95], linspace(0, 1, 256)))
clim(ax, [0 1])
h_real = area(ax, dem_data.dem_lon(lon_id), dem_full(dem_row, lon_id) / 1000, FaceColor=[0.55 0.55 0.55], EdgeColor="none", FaceAlpha=0.6);
edges_lon = [g.lon(x_id) - g.dx / g.dlon / 2, g.lon(x_id(end)) + g.dx / g.dlon / 2];
h_grid = stairs(ax, edges_lon, [ground(x_id, y_id); ground(x_id(end), y_id)] / 1000, Color=[0.75 0.1 0.1], LineWidth=1.8);
h_smooth = plot(ax, g.lon(x_id), smooth_ground(x_id, y_id) / 1000, "k--", LineWidth=1.5);
h_top = yline(ax, grid_top / 1000, "k:", LineWidth=1.2);
set(ax, YDir="normal", FontSize=font, Layer="top")
xlim(ax, map_box(1:2))
ylim(ax, [0 3.2])
xlabel(ax, "Longitude (°)")
ylabel(ax, "Height (km asl)")
title(ax, sprintf("Section at %.1f° N: active cells (blue, shade = air share)", section_lat), FontSize=font + 1, FontWeight="bold")
legend(ax, [h_real h_grid h_smooth h_top], ["Real ground (30"")" "Grid ground (4 km)" "Ground smoothed over 100 km" "Grid top (3000 m asl)"], ...
    Location="northeast", FontSize=font - 1)

% Rockies edge (39.75 N, Denver): the same sampling (nearest at 2 km, mean
% of 2 x 2 to 4 km) on a local grid; the us_east domain ends at -104
rockies_lat = 39.75;
rockies_lon = [-106.5 -102.5];
km_lon = 111.32 * cosd(rockies_lat);
fine_lon = rockies_lon(1):2 / km_lon:rockies_lon(2);
fine_lat = rockies_lat + [-1 1] / 111.32;
[~, fine_col] = min(abs(dem_data.dem_lon(:) - fine_lon), [], 1);
[~, fine_row] = min(abs(dem_data.dem_lat(:) - fine_lat), [], 1);
fine = dem_full(fine_row, fine_col);                  % [2 n_fine]
n_coarse = floor(numel(fine_lon) / 2);
coarse = mean(reshape(fine(:, 1:2 * n_coarse), 4, n_coarse), 1);
coarse_lon = mean(reshape(fine_lon(1:2 * n_coarse), 2, n_coarse), 1);
real_col = find(dem_data.dem_lon >= rockies_lon(1) & dem_data.dem_lon <= rockies_lon(2));
[~, real_row] = min(abs(dem_data.dem_lat - rockies_lat));
ax = nexttile(layout);
area(ax, dem_data.dem_lon(real_col), dem_full(real_row, real_col) / 1000, FaceColor=[0.55 0.55 0.55], EdgeColor="none", FaceAlpha=0.6)
hold(ax, "on")
coarse_edges = [coarse_lon - 2 / km_lon, coarse_lon(end) + 2 / km_lon];
stairs(ax, coarse_edges, [coarse coarse(end)] / 1000, Color=[0.75 0.1 0.1], LineWidth=1.8)
yline(ax, grid_top / 1000, "k:", LineWidth=1.2)
xline(ax, -104, "k--", "us\_east edge", FontSize=font - 2, LabelVerticalAlignment="bottom")
set(ax, FontSize=font, Layer="top")
xlim(ax, rockies_lon)
ylim(ax, [0 4.5])
xlabel(ax, "Longitude (°)")
ylabel(ax, "Height (km asl)")
title(ax, sprintf("Rockies edge, %.2f° N (same sampling)", rockies_lat), FontSize=font + 1, FontWeight="bold")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_K_grid.png"), Resolution=300);
close(fig)
plains = coarse_lon >= -104.2 & coarse_lon <= -103.8;
fprintf("Rockies edge: 4 km ground at -104 %.0f m (air below the grid top %.0f m); ground above the top west of %.2f\n", ...
    mean(coarse(plains)), grid_top - mean(coarse(plains)), max(coarse_lon(coarse >= grid_top)))

%% Radar rows over terrain (figures K2, K3) ---------------------------------------
% Each observation row: F-weighted beam centre (x, y, height asl) over the
% active cells, vertical spread, ground under the footprint (grid ground),
% and its radar's column ground. Rows mostly in active air (F weight on
% active cells > 0.5) at 30-150 km; every frame. Density relative to the
% radar-frame mean (removes the level of each radar and time, keeps the
% vertical structure).
F_active = S.observations.F(:, S.active_id);
[cell_x, cell_y, cell_z] = ind2sub(g.sz, S.active_id);
cell_ground = ground(sub2ind(g.sz(1:2), cell_x, cell_y));
row_weight = full(sum(F_active, 2));
row_x = full(F_active * g.x(cell_x)') ./ row_weight;
row_y = full(F_active * g.y(cell_y)) ./ row_weight;
row_z = full(F_active * z(cell_z)) ./ row_weight;
row_ground = full(F_active * cell_ground) ./ row_weight;
radar_x = g.f_X(S.radar_selected.lon);
radar_y = g.f_Y(S.radar_selected.lat);
radar_column = sub2ind(g.sz(1:2), round(interp1(g.x, 1:g.sz(1), radar_x)), round(interp1(g.y, 1:g.sz(2), radar_y)));
radar_ground = ground(radar_column);
row_use = row_weight > 0.5 & S.observations.row_range_m >= 30e3 & S.observations.row_range_m <= 150e3;
[obs_row, obs_time, obs_value] = find(S.observations.values);
keep = row_use(obs_row);
E = table(obs_row(keep), obs_time(keep), obs_value(keep), VariableNames=["row" "t" "density"]);
E.radar = S.observations.row_radar_id(E.row);
E.z = row_z(E.row);
E.ground = row_ground(E.row);
E.height = E.z - E.ground;                            % beam centre above the ground under it
E.rise = E.ground - radar_ground(E.radar);            % ground above the radar's ground
[radar_frame, ~, ~] = findgroups(E.radar, E.t);
frame_level = splitapply(@mean, E.density, radar_frame);
E.relative = E.density ./ frame_level(radar_frame);
E.log_density = log(max(E.density, 0.1));            % floored as ar_holdout_score
ground_edges = [-inf 300 600 1000 inf];
ground_names = ["< 300 m" "300-600 m" "600-1000 m" "> 1000 m"];
E.ground_class = discretize(E.ground, ground_edges);
fprintf("%d observations (%d rows), %.2f%% over ground > 1000 m\n", height(E), numel(unique(E.row)), 100 * mean(E.ground_class == 4))

% Median relative density by height above ground and above sea level, per
% ground class (bins with at least 300 observations)
height_edges = 0:200:3000;
height_centres = height_edges(1:end - 1) + 100;
profiles = table();
for i_c = 1:numel(ground_names)
    in_class = find(E.ground_class == i_c);
    for coordinate = ["agl" "asl"]
        h = E.height(in_class);
        if coordinate == "asl"
            h = E.z(in_class);
        end
        bin = discretize(h, height_edges);
        ok = ~isnan(bin);
        n = accumarray(bin(ok), 1, [numel(height_centres) 1]);
        value = accumarray(bin(ok), E.relative(in_class(ok)), [numel(height_centres) 1], @median, NaN);
        value(n < 300) = NaN;
        profiles = [profiles; table(repmat(ground_names(i_c), numel(height_centres), 1), ...
            repmat(coordinate, numel(height_centres), 1), height_centres', value, n, ...
            VariableNames=["ground" "coordinate" "height_m" "relative_density" "n"])]; %#ok<AGROW>
    end
end
writetable(profiles, fullfile(figure_dir, "AR_K_profiles.csv"))

% Same beam, different ground: rows of one radar, frame, sweep and 5 km
% range bin have the same beam height asl (within ~70 m) and differ only in
% azimuth. Density relative to the bin median against the ground relative
% to the bin median, in bands of beam height above the bin's median
% ground. Level flight: 1 at any ground; terrain following: denser over
% higher ground (closer to the dense low layer), as the AGL curve.
E.range_bin = floor(S.observations.row_range_m(E.row) / 5000);
beam = findgroups(E.radar, E.t, S.observations.row_elev_deg(E.row), E.range_bin);
beam_ground = splitapply(@median, E.ground, beam);
beam_relative = splitapply(@median, E.relative, beam);
beam_z = splitapply(@median, E.z, beam);
E.ground_anomaly = E.ground - beam_ground(beam);
E.beam_ratio = log(max(E.relative, 1e-3) ./ max(beam_relative(beam), 1e-3));
E.beam_height = beam_z(beam) - beam_ground(beam);
anomaly_edges = [-500 -300 -100 100 300 500 900 1500];
anomaly_centres = [-400 -200 0 200 400 700 1200];
band_edges = [0 600 1200 1800 3000];
band_names = ["0-600 m" "600-1200 m" "1200-1800 m" "1800-3000 m"];
lowland = profiles(profiles.ground == ground_names(1) & profiles.coordinate == "agl", :);
same_beam = table();
for i_b = 1:numel(band_names)
    in_band = find(E.beam_height >= band_edges(i_b) & E.beam_height < band_edges(i_b + 1));
    bin = discretize(E.ground_anomaly(in_band), anomaly_edges);
    ok = ~isnan(bin);
    n = accumarray(bin(ok), 1, [numel(anomaly_centres) 1]);
    value = exp(accumarray(bin(ok), E.beam_ratio(in_band(ok)), [numel(anomaly_centres) 1], @median, NaN));
    value(n < 300) = NaN;
    % AGL expectation from the lowland profile by height above ground
    band_height = median(E.beam_height(in_band));
    profile_at = @(h) interp1(lowland.height_m, lowland.relative_density, h, "linear", "extrap");
    agl_expected = profile_at(max(band_height - anomaly_centres', 100)) ./ profile_at(band_height);
    same_beam = [same_beam; table(repmat(band_names(i_b), numel(anomaly_centres), 1), anomaly_centres', value, ...
        agl_expected, n, VariableNames=["band" "ground_anomaly_m" "relative_density" "agl_expected" "n"])]; %#ok<AGROW>
end
writetable(same_beam, fullfile(figure_dir, "AR_K_same_beam.csv"))

fig = bmm4d.plot.create_export_figure(3600, 1350, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 1, 3, TileSpacing="compact", Padding="compact");   % "tight" lets the last legend overflow
class_colors = copper(numel(ground_names) + 1);
class_colors = class_colors(2:end, :);
for coordinate = ["agl" "asl"]
    ax = nexttile(layout);
    hold(ax, "on")
    for i_c = 1:numel(ground_names)
        p = profiles(profiles.ground == ground_names(i_c) & profiles.coordinate == coordinate, :);
        plot(ax, p.relative_density, p.height_m / 1000, "-o", Color=class_colors(i_c, :), ...
            MarkerFaceColor=class_colors(i_c, :), LineWidth=1.8, MarkerSize=5)
    end
    set(ax, XScale="log", FontSize=font, Box="on", XTick=[0.1 0.3 1 3 10], XTickLabel=["0.1" "0.3" "1" "3" "10"])
    xlim(ax, [0.05 15])
    ylim(ax, [0 3])
    grid(ax, "on")
    xlabel(ax, "Density / radar-frame mean")
    if coordinate == "agl"
        ylabel(ax, "Beam height above ground (km)")
        title(ax, "By height above ground", FontSize=font + 1, FontWeight="bold")
        legend(ax, "ground " + ground_names, Location="northeast", FontSize=font - 1)
    else
        ylabel(ax, "Beam height above sea level (km)")
        title(ax, "By height above sea level", FontSize=font + 1, FontWeight="bold")
    end
end
ax = nexttile(layout);
hold(ax, "on")
band_colors = winter(numel(band_names));
for i_b = 1:numel(band_names)
    p = same_beam(same_beam.band == band_names(i_b), :);
    plot(ax, p.ground_anomaly_m, p.relative_density, "-o", Color=band_colors(i_b, :), MarkerFaceColor=band_colors(i_b, :), ...
        LineWidth=1.8, MarkerSize=5)
end
for i_b = 1:numel(band_names)
    p = same_beam(same_beam.band == band_names(i_b), :);
    plot(ax, p.ground_anomaly_m(p.n >= 300), p.agl_expected(p.n >= 300), "--", Color=band_colors(i_b, :), LineWidth=1.2, ...
        HandleVisibility="off")
end
yline(ax, 1, "k:", LineWidth=1.5)
plot(ax, NaN, NaN, "k--", LineWidth=1.2)
set(ax, YScale="log", FontSize=font, Box="on", YTick=[0.5 0.7 1 1.5 2 3], XTick=-400:400:1200)
xlim(ax, [-500 1300])
ylim(ax, [0.13 4])
grid(ax, "on")
xlabel(ax, "Ground - beam's median ground (m)")
ylabel(ax, "Density / median of the same beam")
title(ax, "Same beam, different ground", FontSize=font + 1, FontWeight="bold")
band_legend = legend(ax, [band_names "level (ASL)" "following (AGL)"], Location="south", NumColumns=2, FontSize=font - 3);
band_legend.Title.String = "Beam above ground";
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_K_data.png"), Resolution=300);
close(fig)

% Which height coordinate makes each radar-frame's rows collapse onto one
% profile? Each row's height is mapped to the radar's column by the rule
% (the inverse of a link from the radar to the row's column), the
% radar-frame profile is the median log density per 100 m of that height,
% and the score is the median |log residual|, on rows whose ground differs
% from the radar's by >= 300 m (the rows where the rules differ) and on
% all rows. Smooth L: height above the ground smoothed over L (0 = AGL on
% the grid ground, infinite = ASL). Hybrid H: AGL near the ground fading
% to level at H (0 = ASL, infinite = AGL).
terrain_rows = abs(E.rise) >= 300;
lengths_km = [0 8 16 25 50 75 100 150 200 300 inf];
fades_m = [0 250 500 750 1000 1500 2000 3000 5000 inf];
coordinate_scores = table();
for length_km = lengths_km
    if length_km == 0
        reference = ground;
    elseif isinf(length_km)
        reference = zeros(size(ground));
    else
        reference = imgaussfilt(ground, length_km * 1e3 / g.dx, Padding="replicate");
    end
    row_reference = full(F_active * reference(sub2ind(g.sz(1:2), cell_x, cell_y))) ./ row_weight;
    mapped = E.z - row_reference(E.row) + reference(radar_column(E.radar));
    coordinate_scores = [coordinate_scores; collapse_score(E, radar_frame, mapped, terrain_rows, "smooth", length_km)]; %#ok<AGROW>
end
for fade = fades_m
    h = max(E.height, 0);
    if isinf(fade)
        mapped = h;
    else
        mapped = h + E.rise;                          % level
        follow = E.rise < fade & h < fade - E.rise;
        mapped(follow) = h(follow) ./ (1 - E.rise(follow) / fade);
    end
    mapped = mapped + radar_ground(E.radar);          % height asl in the radar's column, as above
    coordinate_scores = [coordinate_scores; collapse_score(E, radar_frame, mapped, terrain_rows, "hybrid", fade)]; %#ok<AGROW>
end
disp(coordinate_scores)
writetable(coordinate_scores, fullfile(figure_dir, "AR_K_coordinate.csv"))

fig = bmm4d.plot.create_export_figure(3600, 1250, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
rules = ["smooth" "hybrid"];
rule_titles = ["Height above the ground smoothed over L" "Terrain following, fading to level at H above ground"];
rule_labels = ["Smoothing length L (km)" "Fade height H (m)"];
end_positions = {[4 600], [125 10000]};               % where L or H = 0 and infinite are drawn
tick_values = {[8 25 50 100 200], [250 500 1000 2000 5000]};
end_names = {["AGL" "ASL"], ["ASL" "AGL"]};
for i_r = 1:2
    ax = nexttile(layout);
    hold(ax, "on")
    s = coordinate_scores(coordinate_scores.rule == rules(i_r), :);
    position = s.parameter;
    position(position == 0) = end_positions{i_r}(1);
    position(isinf(position)) = end_positions{i_r}(2);
    plot(ax, position, s.terrain_rows, "-o", Color=[0.80 0.25 0.10], MarkerFaceColor=[0.80 0.25 0.10], LineWidth=1.8, MarkerSize=6)
    plot(ax, position, s.all_rows, "-s", Color=[0.95 0.60 0.20], MarkerFaceColor=[0.95 0.60 0.20], LineWidth=1.8, MarkerSize=6)
    set(ax, XScale="log", FontSize=font, Box="on", XTick=[end_positions{i_r}(1) tick_values{i_r} end_positions{i_r}(2)], ...
        XTickLabel=[end_names{i_r}(1) string(tick_values{i_r}) end_names{i_r}(2)])
    xlim(ax, end_positions{i_r} .* [0.8 1.25])
    ylim(ax, [0.3 0.65])
    grid(ax, "on")
    xlabel(ax, rule_labels(i_r))
    if i_r == 2
        ax.YTickLabel = [];                           % shared y axis
    end
    if i_r == 1
        ylabel(ax, "Median |log(density / radar-frame profile)|")
        legend(ax, ["rows over ground >= 300 m from the radar's" "all rows"], Location="northeast", FontSize=font - 1)
    end
    title(ax, rule_titles(i_r), FontSize=font + 1, FontWeight="bold")
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_K_coordinate.png"), Resolution=300);
close(fig)


%% Which rows to trust over the ridges (figure K4) ------------------------------------
% The cleaning drops gates whose terrain horizon (beam-blockage map, 30"
% DEM) is above the beam's lower half-power edge (elevation - 0.475 deg),
% so a kept gate loses at most the part of the two-way beam below that
% edge (~5% for a Gaussian beam). Over the ridges this leaves rows only
% well above the ground. Left: lowest beam centre above the ground seen in
% each 4 km column (rows mostly in active air, any range, any frame).
% Middle: heights above ground of the rows over each ground class. Right:
% share of the 0.5 and 0.9 deg gates blocked, by range, for the mountain
% radars (map looked up at compass azimuth, as the fixed filter).
row_use_all = row_weight > 0.5 & any(S.observations.values, 2);
row_column = sub2ind(g.sz(1:2), min(max(round(row_x / g.dx) + 1, 1), g.sz(1)), min(max(round(row_y / g.dy) + 1, 1), g.sz(2)));
row_height = row_z - row_ground;
lowest_seen = accumarray(row_column(row_use_all), row_height(row_use_all), [prod(g.sz(1:2)) 1], @min, NaN);
lowest_seen = reshape(lowest_seen, g.sz(1:2));
ridge = ground > 1000 & any(S.active_mask_3d, 3);
fprintf("ridge columns (> 1000 m): %d, seen %.0f%%, lowest beam above ground: median %.0f m, 10%% %.0f m\n", nnz(ridge), ...
    100 * mean(isfinite(lowest_seen(ridge))), median(lowest_seen(ridge), "omitnan"), quantile(lowest_seen(ridge), 0.1))

bb = load(bmm4d.util.data_path("reference", "beam_blockage", "beam_blockage.mat"));
half_beam = 0.475;
mountain = ["KMRX" "KGSP" "KFCX" "KRLX" "KJKL" "KHTX"];
compass = 0:359;
blocked_share = table();
for radar_name = mountain
    horizon = bb.bb(:, :, bb.bb_radar == radar_name);
    horizon_compass = interp1(bb.bb_az, horizon', mod(90 - compass, 360), "nearest")';   % [range az]
    for elevation = [0.5 0.9]
        share = mean(horizon_compass > elevation - half_beam, 2);
        blocked_share = [blocked_share; table(repmat(radar_name, numel(bb.bb_range), 1), repmat(elevation, numel(bb.bb_range), 1), ...
            bb.bb_range(:), share, VariableNames=["radar" "elevation" "range_m" "blocked_share"])]; %#ok<AGROW>
    end
end
writetable(blocked_share(mod(blocked_share.range_m, 10000) == 0, :), fullfile(figure_dir, "AR_K_blockage.csv"))

fig = bmm4d.plot.create_export_figure(3600, 1350, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
% "compact" spacing: with "tight" the map's colorbar and the last panel overflow
layout = tiledlayout(fig, 1, 5, TileSpacing="compact", Padding="compact");
ax = nexttile(layout, [1 3]);
map_bounds = [-86.5 -79 34 38.8];
imagesc(ax, g.lon, g.lat, lowest_seen', AlphaData=isfinite(lowest_seen'));   % [y x] for imagesc
hold(ax, "on")
contour(ax, g.lon, g.lat, ground', [1000 1000], "Color", [0.8 0.1 0.1], LineWidth=1.2)
in_view = S.radar_selected.lon >= map_bounds(1) & S.radar_selected.lon <= map_bounds(2) & ...
    S.radar_selected.lat >= map_bounds(3) & S.radar_selected.lat <= map_bounds(4);
plot(ax, S.radar_selected.lon(in_view), S.radar_selected.lat(in_view), "k^", MarkerFaceColor="w", MarkerSize=7)
text(ax, S.radar_selected.lon(in_view) + 0.08, S.radar_selected.lat(in_view) + 0.08, string(S.radar_selected.name(in_view)), FontSize=font - 3)
colormap(ax, flipud(gray(256)))
clim(ax, [0 2000])
set(ax, YDir="normal", Color="w", FontSize=font, DataAspectRatio=[1 cosd(36.5) 1])
axis(ax, map_bounds)
xlabel(ax, "Longitude (°)")
ylabel(ax, "Latitude (°)")
title(ax, "Lowest beam centre above ground (m); red: ground 1000 m", FontSize=font + 1, FontWeight="bold")
colorbar(ax, Color="k");

ax = nexttile(layout);
hold(ax, "on")
height_bins = 0:100:3000;
for i_c = 1:numel(ground_names)
    in_class = row_use_all & discretize(row_ground, ground_edges) == i_c;
    counts = histcounts(row_height(in_class), height_bins);
    plot(ax, 100 * counts / sum(counts), height_bins(1:end - 1) / 1000 + 0.05, "-", Color=class_colors(i_c, :), LineWidth=2)
end
set(ax, FontSize=font, Box="on")
grid(ax, "on")
ylim(ax, [0 3])
xlabel(ax, "Rows (% per 100 m)")
ylabel(ax, "Beam centre above ground (km)")
title(ax, "Heights seen", FontSize=font + 1, FontWeight="bold")
class_legend = legend(ax, ground_names, Location="northeast", FontSize=font - 2);
class_legend.Title.String = "Ground";

ax = nexttile(layout);
hold(ax, "on")
radar_colors = lines(numel(mountain));
for i_r = 1:numel(mountain)
    for elevation = [0.5 0.9]
        b = blocked_share(blocked_share.radar == mountain(i_r) & blocked_share.elevation == elevation, :);
        style = "-";
        if elevation == 0.9
            style = "--";
        end
        plot(ax, b.range_m / 1000, 100 * b.blocked_share, style, Color=radar_colors(i_r, :), LineWidth=1.6, ...
            HandleVisibility=string(matlab.lang.OnOffSwitchState(elevation == 0.5)))
    end
end
set(ax, FontSize=font, Box="on")
grid(ax, "on")
xlim(ax, [0 200])
ylim(ax, [0 200])
yticks(ax, 0:20:100)
xlabel(ax, "Range (km)")
ylabel(ax, "Azimuths blocked (%)")
title(ax, ["Blocked gates" "0.5° solid, 0.9° dashed"], FontSize=font + 1, FontWeight="bold")
legend(ax, mountain, Location="northeast", FontSize=font - 3)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_K_trust.png"), Resolution=300);
close(fig)
disp(groupsummary(blocked_share(blocked_share.range_m >= 60e3 & blocked_share.range_m <= 120e3, :), ["radar" "elevation"], "mean", "blocked_share"))

%% Links of each rule over terrain ---------------------------------------------------
% As AR_G: uphill sources (destination ground > 300 m above the source's,
% from the sigma links, which give every source cell a link) and lee cells
% (in a column reached from ground > 300 m higher). Per rule: the share of
% uphill sources that keep a link, of lee cells that receive one, and the
% number of detour links.
nxy = prod(g.sz(1:2));
column_of = @(id) mod(double(id) - 1, nxy) + 1;
active_cell = false(g.sz);
active_cell(S.active_id) = true;
steps_sigma = ar_steps(g, S.radar_scans, dt_seconds, S.active_mask_3d, inv, "sigma");
link_rules = struct("name", {"ASL (production)", "AGL", "sigma", "hybrid H 1000 m", "smooth L 100 km", ...
    "smooth L 100 km + detour", "ASL + detour"}, "steps", {"asl", "agl", "sigma", "hybrid", "smooth", "smooth", "asl"}, ...
    "options", {{}, {}, {}, {"fade_height_m", 1000}, {"smooth_length_m", 100e3}, {"smooth_length_m", 100e3, "detour", true}, ...
    {"detour", true}});
links = table();
for rule = link_rules
    steps_rule = ar_steps(g, S.radar_scans, dt_seconds, S.active_mask_3d, inv, rule.steps, rule.options{:});
    counts = zeros(1, 6);
    for i_t = 1:numel(steps_rule)
        s = steps_sigma{i_t};
        rise = ground(column_of(s.dest_id_valid)) - ground(column_of(s.source_id_valid));
        uphill = unique(s.source_id_valid(rise > 300));
        lee_column = false(g.sz(1:2));
        lee_column(unique(column_of(s.dest_id_valid(rise < -300)))) = true;
        lee = find(active_cell & repmat(lee_column, 1, 1, g.sz(3)));
        m = steps_rule{i_t};
        n_detour = 0;
        if isfield(m, "detour_valid")
            n_detour = nnz(m.detour_valid);
        end
        counts = counts + [nnz(ismember(uphill, m.source_id_valid)) numel(uphill) nnz(ismember(lee, m.dest_id_valid)) ...
            numel(lee) numel(m.source_id_valid) n_detour];
    end
    links = [links; table(string(rule.name), counts(1) / counts(2), counts(3) / counts(4), counts(5), counts(6), ...
        VariableNames=["rule" "uphill_sources_linked" "lee_cells_reached" "links" "detour_links"])]; %#ok<AGROW>
end
disp(links)
writetable(links, fullfile(figure_dir, "AR_K_links.csv"))

%% Solves and held-out radars (batch job, ~70 min) -------------------------------------
% Vertical rules (links and their prior), each solved with the production
% settings, then each mountain radar left out in turn (as AR_G).
variants = struct( ...
    "name", {"ASL (production)", "AGL", "hybrid H 1000 m", "smooth L 25 km", "smooth L 50 km", "smooth L 100 km", ...
    "smooth L 100 km + detour", "smooth L 150 km + detour", "ASL + detour"}, ...
    "steps", {"asl", "agl", "hybrid", "smooth", "smooth", "smooth", "smooth", "smooth", "asl"}, ...
    "steps_options", {{}, {}, {"fade_height_m", 1000}, {"smooth_length_m", 25e3}, {"smooth_length_m", 50e3}, ...
    {"smooth_length_m", 100e3}, {"smooth_length_m", 100e3, "detour", true}, {"smooth_length_m", 150e3, "detour", true}, ...
    {"detour", true}}, ...
    "prior", {"asl", "agl", "hybrid", "smooth", "smooth", "smooth", "smooth", "smooth", "asl"}, ...
    "prior_options", {{}, {}, {"fade_height_m", 1000}, {"smooth_length_m", 25e3}, {"smooth_length_m", 50e3}, ...
    {"smooth_length_m", 100e3}, {"smooth_length_m", 100e3}, {"smooth_length_m", 150e3}, {}});
held_out = intersect(["KMRX" "KGSP" "KFCX" "KRLX" "KJKL" "KHTX"], string(S.radar_selected.name), "stable");
solve_job = batch(@ar_terrain_solves, 0, {c, variants, held_out, section_time, results_file, 6}, Pool=0, ...
    AutoAddClientPath=true, CurrentFolder=project_root);

%% Scores of the solves and the held-out radars ----------------------------------------
wait(solve_job)
results = load(results_file);
disp(results.scores)
writetable(results.scores, fullfile(figure_dir, "AR_K_solves.csv"))
writetable(results.holdout, fullfile(figure_dir, "AR_K_holdout.csv"))
holdout_mean = groupsummary(results.holdout, "solve", "mean", ["relative_error" "median_abs_log_ratio"]);
[~, order] = ismember(string({results.variants.name}), holdout_mean.solve);
disp(holdout_mean(order, :))
disp(unstack(results.holdout(:, ["radar" "solve" "relative_error"]), "relative_error", "solve"))

%% Section across the Smokies and maps (figures K5, K6) ------------------------------
% West-east section at 35.6 N, 06:00 UTC. Top: the observed rows whose beam
% centre is within 6 km of the section (any radar, any sweep), at their
% beam-centre height, coloured by density. Below: the reconstruction of
% each rule, density by height above sea level; real ground (30") in gray,
% grid ground (4 km) in red. imagesc wants [z x], so the [x z] slice is
% transposed here.
shown = ["ASL (production)" "AGL" "hybrid H 1000 m" "smooth L 50 km"];
[~, shown_id] = ismember(shown, string({results.variants.name}));
i_section = find(c.times == section_time);
[~, y_id] = min(abs(g.lat - section_lat));
x_id = find(g.lon >= map_box(1) & g.lon <= map_box(2));
sections = cell(numel(shown), 1);
for i_v = 1:numel(shown)
    frame = nan(g.sz);
    frame(S.active_id) = results.section_fields{shown_id(i_v)};
    sections{i_v} = squeeze(frame(x_id, y_id, :));
end
near_section = row_weight > 0.5 & abs(row_y - g.f_Y(section_lat)) <= 6e3 & ...
    g.f_lon(row_x) >= map_box(1) & g.f_lon(row_x) <= map_box(2) & S.observations.values(:, i_section) ~= 0;
observed = full(S.observations.values(near_section, i_section));
limit = quantile([cell2mat(cellfun(@(s) s(:), sections, UniformOutput=false)); observed], 0.99);

fig = bmm4d.plot.create_export_figure(3600, 2700, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, numel(shown) + 1, 1, TileSpacing="tight", Padding="tight");
for i_p = 0:numel(shown)
    ax = nexttile(layout);
    hold(ax, "on")
    if i_p == 0
        [~, order] = sort(observed);                  % densest on top
        lon_near = g.f_lon(row_x(near_section));
        z_near = row_z(near_section);
        scatter(ax, lon_near(order), z_near(order) / 1000, 14, observed(order), "filled")
        panel_title = sprintf("Observed rows within 6 km of %.1f° N (beam centre), %s UTC", section_lat, string(section_time, "HH:mm"));
    else
        imagesc(ax, g.lon(x_id), z / 1000, sections{i_p}', AlphaData=~isnan(sections{i_p}'));
        panel_title = shown(i_p);
    end
    area(ax, dem_data.dem_lon(lon_id), dem_full(dem_row, lon_id) / 1000, FaceColor=[0.55 0.55 0.55], EdgeColor="none")
    stairs(ax, edges_lon, [ground(x_id, y_id); ground(x_id(end), y_id)] / 1000, Color=[0.75 0.1 0.1], LineWidth=1.4)
    set(ax, YDir="normal", Color="w", FontSize=font, Box="on", Layer="top")
    colormap(ax, parula)
    clim(ax, [0 limit])
    xlim(ax, map_box(1:2))
    ylim(ax, [0 3])
    ylabel(ax, "Height (km asl)")
    title(ax, panel_title, FontSize=font + 1, FontWeight="bold")
    if i_p < numel(shown)
        ax.XTickLabel = [];
    end
end
xlabel(ax, "Longitude (°)")
cb = colorbar(ax, Color="k");
cb.Layout.Tile = "east";
cb.Label.String = "Birds km^{-3}";
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_K_section.png"), Resolution=300);
close(fig)

% Column maps at 04:30, 06:00 and 07:30 with the ground
map_names = ["ASL (production)" "AGL" "hybrid 1000 m" "smooth 50 km"];
map_results = cellfun(@(name, column) struct("name", name, "column", column), cellstr(map_names), ...
    results.columns(shown_id)', UniformOutput=false);
ar_plot_maps(S, map_results, [1 4 7], fullfile(figure_dir, "AR_K_maps.png"), ...
    bounds=[-88 -79 34 38.5], extra=ground, extra_title="Ground (m)", extra_limits=[0 1600]);

%% Check of the chosen rule: smoothed ground 50 km with climb (figures K7, K8) ----------
% The rule kept: links and prior by height above the ground smoothed over
% 50 km; a link into the real (4 km) ground goes to the lowest air cell of
% its destination column. Solved beside ASL and AGL on the Appalachian
% case and on the Plains case (to the -104 edge, ground to 1900 m), to look
% for artefacts: column ratio to ASL, sections, the lowest cell against the
% cell 300 m above it (a layer piled into the lowest cell would show), and
% the edge score of fix 3.
rules = struct("name", {"ASL (production)", "AGL", "smooth 50 km + climb"}, "steps", {"asl", "agl", "smooth"}, ...
    "steps_options", {{}, {}, {"smooth_length_m", 50e3, "climb", true}}, "prior", {"asl", "agl", "smooth"}, ...
    "prior_options", {{}, {}, {"smooth_length_m", 50e3}});
c_plains = ar_test_case("ar_plains_bb", [-104 -94 37.5 47.5], times);
check_cases = {S, ar_load(c_plains)};
check_names = ["appalachians" "plains"];
check_lats = {[35.0 35.6 36.2 37.0], [39.0 41.0 43.0 45.0]};
check_lons = {[-85.5 -80.5], [-106 -94]};
for i_c = 1:2
    S_check = check_cases{i_c};
    dt_check = seconds(diff(S_check.c.times));
    check_results = cell(numel(rules), 1);
    edges = zeros(numel(rules), 1);
    for i_r = 1:numel(rules)
        r = rules(i_r);
        steps_check = ar_steps(S_check.g, S_check.radar_scans, dt_check, S_check.active_mask_3d, inv, r.steps, r.steps_options{:});
        prior_check = ar_profile_prior(S_check.g, S_check.radar_scans, S_check.radar_profiles, S_check.active_mask_3d, ...
            inv.prior_length_scale_m, r.prior, r.prior_options{:});
        check_results{i_r} = ar_solve(S_check, r.name, steps=steps_check, prior=prior_check);
        edge = ar_edge_score(S_check, check_results{i_r}.column);
        edges(i_r) = edge.edge;
    end
    skin = ar_check_figures(S_check, check_results, check_lats{i_c}, check_lons{i_c}, section_time, ...
        string(fullfile(figure_dir, "AR_K_check_" + check_names(i_c))));
    skin.edge = edges;
    disp(skin)
    writetable(skin, fullfile(figure_dir, "AR_K_check_" + check_names(i_c) + ".csv"))
end

function score = collapse_score(E, radar_frame, mapped_height, terrain_rows, rule, parameter)
% Median |log residual| from each radar-frame's profile in the mapped height
bin = findgroups(radar_frame, floor(mapped_height / 100));
profile = splitapply(@median, E.log_density, bin);
residual = abs(E.log_density - profile(bin));
score = table(string(rule), parameter, median(residual(terrain_rows)), median(residual), ...
    VariableNames=["rule" "parameter" "terrain_rows" "all_rows"]);
end
