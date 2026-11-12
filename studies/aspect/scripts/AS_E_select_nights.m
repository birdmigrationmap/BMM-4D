%% Aspect E: choose the nights of the multi-night aspect fit (dev)
% One night (us_east, 21-22 Oct 2026, tracks mostly toward S) measures the
% aspect function precisely but cannot tell how it changes with season,
% region, speed or year, nor separate it from the fixed azimuthal patterns
% of each radar (blockage, calibration by sector): with one track
% direction both sit on the same axis. Nights with different track
% directions separate them (the aspect turns with the track, the radar
% patterns do not). The side-on term repeats every 180 deg, so spring and
% autumn (opposite tracks) do not help it; tracks turned by 30-90 deg do.
% Nights are taken from the Dark Ecology daily archive (2021-2025, one row
% per radar and night) for the 107 us_east radars:
% - spring (15 Mar-31 May) or autumn (15 Aug-15 Nov);
% - at least 90 radars report, each with less than 10% missing;
% - rain on less than 10% of the night at 70% of them or more;
% - median night-mean bird reflectivity of at least 200 cm^2 km^-2.
% The qualifying nights are split by season and domain-mean track (30 deg
% sectors); nights are taken in turn from each sector, strongest first, at
% least 4 days from a chosen night and at most 4 per year, so rare track
% directions come in early. No download: this only reads the archive.
% Frames: hourly, from 2 h after sunset at 37 N, 85 W, 6 frames (fewer
% if sunrise is less than 1 h after the last one).
% Outputs: <data>/aspect/multi_night/nights.csv (read by AS_F) and
% figures/AS_E_nights.csv / AS_E_nights.png.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
data_root = bmm4d.util.local_data_dir();
multi_dir = fullfile(data_root, "aspect", "multi_night");
figure_dir = fullfile(project_root, "studies", "aspect", "figures");
if ~isfolder(multi_dir)
    mkdir(multi_dir)
end
us = case_settings("us_east");
n_target = 16;
n_frames = 6;

%% Read the Dark Ecology nights of the us_east radars
daily_dir = fullfile(data_root, "gate_density", "darkeco", "daily");
nights = table();
for y = 2021:2025
    daily = readtable(fullfile(daily_dir, y + "-daily.csv"), TextType="string");
    nights = [nights; daily(daily.period == "night" & ismember(daily.station, us.radars), :)]; %#ok<AGROW>
end
nights.mean_refl = nights.reflectivity_hours ./ nights.period_length;     % cm^2 km^-2
day_of_year = day(nights.date, "dayofyear");
nights.season = repmat("", height(nights), 1);
nights.season(day_of_year >= 74 & day_of_year <= 151) = "spring";
nights.season(day_of_year >= 227 & day_of_year <= 319) = "autumn";
nights = nights(nights.season ~= "", :);

%% Score each night over the radars
[date_id, dates] = findgroups(nights.date);
ok = nights.fraction_missing < 0.1;
% Domain-mean track and its coherence, weighted by the night-mean reflectivity
has_track = ok & isfinite(nights.u) & isfinite(nights.v) & isfinite(nights.mean_refl);
weight = max(nights.mean_refl, 0) .* has_track;
u = nights.u;
v = nights.v;
u(~has_track) = 0;
v(~has_track) = 0;
axis_vector = exp(2i * deg2rad(nights.direction));
axis_vector(~has_track | ~isfinite(axis_vector)) = 0;
sum_weight = splitapply(@sum, weight, date_id);
screen = table(dates, splitapply(@(s) s(1), nights.season, date_id), splitapply(@sum, ok, date_id), ...
    splitapply(@sum, ok & nights.fraction_rain < 0.1, date_id), ...
    splitapply(@(x) median(x, "omitnan"), nights.mean_refl, date_id), ...
    mod(atan2d(splitapply(@sum, weight .* u, date_id), splitapply(@sum, weight .* v, date_id)), 360), ...
    abs(splitapply(@sum, weight .* axis_vector, date_id)) ./ sum_weight, ...
    splitapply(@sum, weight .* hypot(u, v), date_id) ./ sum_weight, ...
    VariableNames=["date" "season" "n_radar" "n_dry" "median_refl" "track_deg" "axis_coherence" "speed_ms"]);
qualify = screen.n_radar >= 90 & screen.n_dry >= 0.7 * screen.n_radar & screen.median_refl >= 200;
screen = screen(qualify, :);
screen.sector = discretize(screen.track_deg, 0:30:360);
fprintf("%d qualifying nights (%d spring, %d autumn)\n", height(screen), nnz(screen.season == "spring"), ...
    nnz(screen.season == "autumn"))

%% Take nights in turn from each season and track sector
[~, strata] = findgroups(screen(:, ["season" "sector"]));
screen = sortrows(screen, "median_refl", "descend");
[~, stratum] = ismember(screen(:, ["season" "sector"]), strata);
chosen = false(height(screen), 1);
while nnz(chosen) < n_target
    added = false;
    for s = 1:height(strata)
        candidates = find(stratum == s & ~chosen)';
        for i = candidates
            near = any(abs(days(screen.date(chosen) - screen.date(i))) < 4);
            full_year = nnz(year(screen.date(chosen)) == year(screen.date(i))) >= 4;
            if ~near && ~full_year
                chosen(i) = true;
                added = true;
                break
            end
        end
        if nnz(chosen) >= n_target
            break
        end
    end
    if ~added
        break
    end
end
selected = sortrows(screen(chosen, :), "date");

%% Frames, names and the run budget
% Night of local date D: sunset on D at 37 N, 85 W (getSunriseSunset takes
% a time within the night, as case_settings "night_YYYYMMDD")
selected.first_frame = NaT(height(selected), 1);
selected.n_frames = zeros(height(selected), 1);
for i = 1:height(selected)
    [sunset, sunrise] = bmm4d.util.getSunriseSunset(37, -85, selected.date(i) + hours(24 + 85 / 15));
    sunset.TimeZone = "";
    sunrise.TimeZone = "";
    first = dateshift(sunset + hours(2), "start", "hour", "next");
    selected.first_frame(i) = first;
    selected.n_frames(i) = min(n_frames, floor(hours(sunrise - hours(1) - first)) + 1);
end
selected.name = "aspect_" + string(selected.date, "yyyyMMdd");
% Budget from the us_east night: ~8.8 MB per raw volume, cleaning 2.6 h for
% 1869 volumes on 5 MistNet workers (6.3 s per volume on 4 workers), binning
% ~0.5 min per frame
selected.volumes = selected.n_frames * numel(us.radars);
selected.download_gb = selected.volumes * 8.8e-3;
selected.clean_h = selected.volumes * 6.3 / 3600;
disp(selected(:, ["name" "season" "track_deg" "axis_coherence" "median_refl" "n_dry" "first_frame" "n_frames"]))
fprintf("%d nights, %d volumes: about %.0f GB to download (less what is local), %.1f h of cleaning\n", ...
    height(selected), sum(selected.volumes), sum(selected.download_gb), sum(selected.clean_h))
writetable(selected, fullfile(multi_dir, "nights.csv"))
writetable(selected, fullfile(figure_dir, "AS_E_nights.csv"))

%% Plot results
% Track of every qualifying night against the day of year, chosen nights
% marked: the spread of track directions is what separates aspect from
% fixed radar patterns
fig = bmm4d.plot.create_export_figure(2400, 1000, DefaultAxesFontName="Helvetica", ...
    DefaultTextFontName="Helvetica");
ax = axes(fig);
hold(ax, "on");
box(ax, "on")
scatter(ax, day(screen.date, "dayofyear"), screen.track_deg, 30, [0.7 0.7 0.7], "filled")
scatter(ax, day(selected.date, "dayofyear"), selected.track_deg, 90, selected.axis_coherence, "filled", ...
    MarkerEdgeColor="k")
colormap(ax, copper(256))
cb = colorbar(ax);
cb.Label.String = "Track coherence over the radars (axis, 0-1)";
xlabel(ax, "Day of year")
ylabel(ax, "Domain-mean track (deg)")
yticks(ax, 0:45:360)
ylim(ax, [0 360])
set(ax, FontSize=13)
title(ax, sprintf("Qualifying nights 2021-2025 (grey) and the %d chosen", height(selected)), FontSize=14)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AS_E_nights.png"), Resolution=300)
close(fig)
