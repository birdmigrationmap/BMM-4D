%% Sandbox: dry, strong-migration nights over the sandbox domain
% Screens the Dark Ecology daily time series (2021-2025, one row per radar
% and night) for nights when the 28 sandbox radars are dry and migration is
% strong, to compare the radar imprint across nights: what repeats every
% night (calibration) and what changes (propagation, birds).
% A night qualifies if, in spring (15 Mar-31 May) or autumn (15 Aug-15 Nov):
% - at least 26 of the 28 radars report, each with less than 10% missing;
% - rain on less than 5% of the night at 24 radars or more (Dark Ecology
%   rain shares are rarely below 1%: median 6% over radar-nights);
% - median night-mean bird reflectivity of at least 300 cm^2 km^-2.
% Nights are ranked by the median over radars of the night-mean bird
% reflectivity; the 6 strongest per season are kept.
% Outputs: figures/night_screen.csv and night_screen.png (one map per
% night: night-mean bird reflectivity per radar, rain marked).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
radar = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
radar = radar(ismember(radar.name, c.radars), :);
n_radar = height(radar);

%% Read the Dark Ecology nights -------------------------------------------
daily_dir = fullfile(bmm4d.util.local_data_dir(), "gate_density", "darkeco", "daily");
nights = table();
for y = 2021:2025
    daily = readtable(fullfile(daily_dir, y + "-daily.csv"), TextType="string");
    nights = [nights; daily(daily.period == "night" & ismember(daily.station, radar.name), :)]; %#ok<AGROW>
end
nights.mean_refl = nights.reflectivity_hours ./ nights.period_length;     % cm^2 km^-2
day_of_year = day(nights.date, "dayofyear");
nights.season = repmat("", height(nights), 1);
nights.season(day_of_year >= 74 & day_of_year <= 151) = "spring";
nights.season(day_of_year >= 227 & day_of_year <= 319) = "autumn";
nights = nights(nights.season ~= "", :);

%% Score each night over the 28 radars ---------------------------------------
[date_id, dates] = findgroups(nights.date);
ok = nights.fraction_missing < 0.1;
n_ok = splitapply(@sum, ok, date_id);
n_dry = splitapply(@sum, ok & nights.fraction_rain < 0.05, date_id);
max_rain = splitapply(@max, nights.fraction_rain, date_id);
median_refl = splitapply(@(x) median(x, "omitnan"), nights.mean_refl, date_id);
season = splitapply(@(s) s(1), nights.season, date_id);
screen = table(dates, season, n_ok, n_dry, max_rain, median_refl, VariableNames=["date" "season" "n_radar" "n_dry" "max_rain" "median_refl"]);
screen = screen(screen.n_radar >= n_radar - 2 & screen.n_dry >= n_radar - 4 & screen.median_refl >= 300, :);
screen = sortrows(screen, "median_refl", "descend");
fprintf("%d dry nights (%d spring, %d autumn)\n", height(screen), nnz(screen.season == "spring"), nnz(screen.season == "autumn"))
selected = [screen(find(screen.season == "spring", 6), :); screen(find(screen.season == "autumn", 6), :)];
disp(selected)
writetable(selected, fullfile(figure_dir, "night_screen.csv"))

%% Plot results -----------------------------------------------------------
% One map per night: night-mean bird reflectivity at each radar (colour,
% log scale), rain share as a red ring when above 5%.
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, ceil(height(selected) / 3), TileSpacing="tight", Padding="tight");
title(t, "Dry nights over the sandbox radars (Dark Ecology, night means)", FontSize=14)
c_lim = log10([100 20000]);
for i_n = 1:height(selected)
    ax = nexttile(t); hold(ax, "on"); box(ax, "on")
    set(ax, Color=[0.85 0.92 0.97], FontSize=11)
    for s = states'
        patch(ax, s.Lon, s.Lat, [0.96 0.96 0.94], EdgeColor=[0.6 0.6 0.6])
    end
    one = nights(nights.date == selected.date(i_n), :);
    [~, id] = ismember(one.station, radar.name);
    scatter(ax, radar.lon(id), radar.lat(id), 120, log10(max(one.mean_refl, 1)), "filled", MarkerEdgeColor="k")
    wet = one.fraction_rain >= 0.05;
    scatter(ax, radar.lon(id(wet)), radar.lat(id(wet)), 260, "r", LineWidth=1.5)
    colormap(ax, parula); clim(ax, c_lim)
    xlim(ax, c.figure_bounds(1:2)); ylim(ax, c.figure_bounds(3:4))
    title(ax, sprintf("%s, %s (median %.0f)", string(selected.date(i_n), "yyyy-MM-dd"), selected.season(i_n), ...
        selected.median_refl(i_n)), FontWeight="bold", FontSize=12)
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Ticks = 2:4;
cb.TickLabels = ["100" "1000" "10000"];
cb.Label.String = "Night-mean bird reflectivity (cm^2 km^{-2}); red ring: rain > 5%";
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "night_screen.png"), Resolution=150)
close(fig)
