%% Night set A: choose the nights for the v2 studies (dev)
% One set of full nights over the us_east radars (107 radars east of
% -104 deg) serves the open questions of the v2 studies at once:
% - aspect (studies/aspect): tracks turned by 30-90 deg between nights
%   separate the aspect from fixed radar patterns; both seasons flip the
%   front/back term; a range of speeds and heights;
% - near-radar dip (studies/near_radar): many radar-nights, both seasons,
%   weak and strong nights (clear-air and precipitation scan strategies),
%   full nights for the gap-filling solves;
% - terrain (studies/artifacts AR_K): strong dry nights over the
%   Appalachians with tracks across the ridges and along them;
% - prior / coverage edges (studies/artifacts PRIOR_OPTIONS.md): strong
%   dry nights over the sparse Plains radars, and profile anomalies over
%   many nights for the horizontal range;
% - radar patterns without birds (aspect plan B1): a few days of
%   widespread rain, daytime, low sweeps;
% - the paper cases (Key West, Midwest), for the v2 rerun of the paper.
% Source: the Dark Ecology daily archive (2021-2025, one row per radar and
% night or day). No download.
% Core nights: dry over the domain (rain < 10% of the night at 70% of the
% reporting radars or more, at least 90 radars reporting), median
% night-mean bird reflectivity of at least 100 cm^2 km^-2. Three nights in
% each of six migration phases, chosen for diversity: the first is the
% strongest of the phase, each next one the night farthest (in track
% axis, track direction, speed and reflectivity) from every night chosen
% so far, at least 4 days from a chosen night and at most 5 per year.
% Then nights are added until the targeted questions have enough:
% 2 Appalachian nights across the ridges and 2 along them, 3 strong Plains
% nights, all dry in that region.
% Frames: every 30 min from 1 h after sunset at the east edge (37 N,
% 70 W) to 30 min before sunrise at the west edge (37 N, 102 W); each
% radar keeps only its night volumes (bmm4d.config radar_night).
% Outputs: figures/NS_A_nights.csv and figures/NS_A_nights.png.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
data_root = bmm4d.util.local_data_dir();
figure_dir = fullfile(project_root, "studies", "nights", "figures");
if ~isfolder(figure_dir)
    mkdir(figure_dir)
end
us = case_settings("us_east");
radar_file = bmm4d.util.data_path("reference", "radar.csv");
regions = struct("name", {"midwest", "plains", "appalachians"}, ...
    "bounds", {[-100 -85 36 46], [-104 -94 37.5 47.5], [-89 -78 33.5 40]});
for i_r = 1:numel(regions)
    regions(i_r).radars = bmm4d.radar.select_radars_in_bbox(radar_file, regions(i_r).bounds).name;
end
ridge_axis_deg = 40;     % Appalachian ridges run NE-SW
volume_gb = 8.8e-3;      % mean raw volume (us_east night)
cpu_s_per_volume = 26;   % MistNet cleaning, measured on merops-02 and the laptop

%% Read the Dark Ecology archive of the us_east radars
daily_dir = fullfile(data_root, "gate_density", "darkeco", "daily");
daily = table();
for y = 2021:2025
    one_year = readtable(fullfile(daily_dir, y + "-daily.csv"), TextType="string");
    daily = [daily; one_year(ismember(one_year.station, us.radars), :)]; %#ok<AGROW>
end
daily.mean_refl = daily.reflectivity_hours ./ daily.period_length;     % cm^2 km^-2
nights = daily(daily.period == "night", :);
day_of_year = day(nights.date, "dayofyear");
% Migration phases (day of year): spring 15 Mar-20 Apr, 21 Apr-20 May,
% 21 May-10 Jun; autumn 15 Aug-14 Sep, 15 Sep-20 Oct, 21 Oct-20 Nov
phase_edges = [74 111 141 162; 227 258 294 325];
phase_names = ["spring early" "spring peak" "spring late"; "autumn early" "autumn peak" "autumn late"];
nights.phase = repmat("", height(nights), 1);
for s = 1:2
    for p = 1:3
        nights.phase(day_of_year >= phase_edges(s, p) & day_of_year < phase_edges(s, p + 1)) = phase_names(s, p);
    end
end
nights = nights(nights.phase ~= "", :);

%% Score each night over the domain and the regions
[date_id, dates] = findgroups(nights.date);
screen = table(dates, splitapply(@(s) s(1), nights.phase, date_id), VariableNames=["date" "phase"]);
screen = [screen, night_metrics(nights, date_id, true(height(nights), 1), "")];
for i_r = 1:numel(regions)
    in_region = ismember(nights.station, regions(i_r).radars);
    screen = [screen, night_metrics(nights, date_id, in_region, regions(i_r).name + "_")]; %#ok<AGROW>
end
screen.season = extractBefore(screen.phase, " ");
qualify = screen.n_radar >= 90 & screen.dry_share >= 0.7 & screen.median_refl >= 100;
pool = screen(qualify, :);
fprintf("%d qualifying nights of %d\n", height(pool), height(screen))
disp(groupcounts(pool, "phase"))

%% Core nights: three per phase, as different as possible
% Features: the side-on axis (2 track, which repeats every 180 deg), the
% track direction (front/back), speed and log reflectivity
features = @(t) [cosd(2 * t.track_deg), sind(2 * t.track_deg), 0.5 * cosd(t.track_deg), ...
    0.5 * sind(t.track_deg), (t.speed_ms - 9) / 4, (log10(t.median_refl) - 2.7) / 0.4];
f_pool = features(pool);
chosen = false(height(pool), 1);
for phase = reshape(phase_names', 1, [])
    in_phase = pool.phase == phase;
    for k = 1:3
        allowed = in_phase & ~chosen & allowed_nights(pool, chosen);
        if ~any(allowed)
            break
        end
        if k == 1
            score = pool.median_refl;
        else
            score = min(pdist2(f_pool, f_pool(chosen, :)), [], 2);
        end
        score(~allowed) = -inf;
        [~, best] = max(score);
        chosen(best) = true;
    end
end
pool.purpose = repmat("core", height(pool), 1);

%% Targeted nights: terrain and coverage edges
% Appalachians: dry and strong there; across the ridges (|sin| of the
% track to the ridge axis >= 0.7) and along them (<= 0.3)
pool.crossing = abs(sind(pool.appalachians_track_deg - ridge_axis_deg));
good_app = pool.appalachians_dry_share >= 0.8 & pool.appalachians_median_refl >= 300;
[pool, chosen] = add_until(pool, chosen, good_app & pool.crossing >= 0.7, 2, "appalachians", "appalachians_median_refl");
[pool, chosen] = add_until(pool, chosen, good_app & pool.crossing <= 0.3, 2, "appalachians", "appalachians_median_refl");
% Plains: the sparse radars where the coverage edges show
good_plains = pool.plains_dry_share >= 0.8 & pool.plains_median_refl >= 300;
[pool, chosen] = add_until(pool, chosen, good_plains, 3, "plains", "plains_median_refl");
selected = sortrows(pool(chosen, :), "date");
% What each night can answer (every full night serves aspect, near-radar and prior)
selected.uses = repmat("aspect near-radar prior", height(selected), 1);
app = selected.appalachians_dry_share >= 0.8 & selected.appalachians_median_refl >= 300;
selected.uses(app & selected.crossing >= 0.7) = selected.uses(app & selected.crossing >= 0.7) + " terrain-across";
selected.uses(app & selected.crossing <= 0.3) = selected.uses(app & selected.crossing <= 0.3) + " terrain-along";
plains = selected.plains_dry_share >= 0.8 & selected.plains_median_refl >= 300;
selected.uses(plains) = selected.uses(plains) + " edges";

%% Frames and budget of the full nights
selected.first_frame = NaT(height(selected), 1);
selected.last_frame = NaT(height(selected), 1);
for i = 1:height(selected)
    in_night = selected.date(i) + hours(24 + 85 / 15);     % a time within the night
    [sunset_east, ~] = bmm4d.util.getSunriseSunset(37, -70, in_night);
    [~, sunrise_west] = bmm4d.util.getSunriseSunset(37, -102, in_night);
    sunset_east.TimeZone = "";
    sunrise_west.TimeZone = "";
    selected.first_frame(i) = dateshift(sunset_east + hours(1), "start", "hour") + ...
        minutes(30) * ceil(minutes(sunset_east + hours(1) - dateshift(sunset_east + hours(1), "start", "hour")) / 30);
    selected.last_frame(i) = dateshift(sunrise_west - minutes(30), "start", "hour") + ...
        minutes(30) * floor(minutes(sunrise_west - minutes(30) - dateshift(sunrise_west - minutes(30), "start", "hour")) / 30);
end
selected.n_frames = round(minutes(selected.last_frame - selected.first_frame) / 30) + 1;
selected.name = "night_us_" + string(selected.date, "yyyyMMdd");
selected.radar_set = repmat("us_east", height(selected), 1);
selected.step_min = repmat(30, height(selected), 1);
% About 86% of the radar-frames are at night (us_east, 21-22 Oct 2026)
selected.volumes = round(0.86 * selected.n_frames * numel(us.radars));

%% Rain days: radar patterns without birds
% Daytime, widespread rain: share of the reporting radars with rain over
% at least half of the day; April-May and September-October (no snow, a
% high melting layer); the three wettest days at least 30 days apart in
% different years. Six hourly frames 15-20 UTC.
day_rows = daily(daily.period == "day" & daily.fraction_missing < 0.1, :);
day_rows = day_rows(ismember(month(day_rows.date), [4 5 9 10]), :);
[day_id, day_dates] = findgroups(day_rows.date);
wet = table(day_dates, splitapply(@numel, day_rows.fraction_rain, day_id), ...
    splitapply(@(r) mean(r >= 0.5), day_rows.fraction_rain, day_id), VariableNames=["date" "n_radar" "wet_share"]);
wet = sortrows(wet(wet.n_radar >= 90, :), "wet_share", "descend");
rain_days = false(height(wet), 1);
for i = 1:height(wet)
    far = all(abs(days(wet.date(rain_days) - wet.date(i))) >= 30);
    new_year = ~any(year(wet.date(rain_days)) == year(wet.date(i)));
    if far && new_year
        rain_days(i) = true;
    end
    if nnz(rain_days) == 3
        break
    end
end
rain = wet(rain_days, :);
rain_rows = table(rain.date, repmat("rain day", height(rain), 1), "rain_us_" + string(rain.date, "yyyyMMdd"), ...
    repmat("us_east", height(rain), 1), rain.date + hours(15), rain.date + hours(20), repmat(6, height(rain), 1), ...
    repmat(60, height(rain), 1), round(6 * numel(us.radars) * ones(height(rain), 1)), rain.wet_share, ...
    VariableNames=["date" "purpose" "name" "radar_set" "first_frame" "last_frame" "n_frames" "step_min" "volumes" "wet_share"]);
disp(rain_rows)

%% The reference night and the paper cases (v2 rerun of the paper)
% The us_east night of 21-22 Oct 2026 (every dev study so far) is in the
% set too, so the cluster products include it.
reference = case_settings("us_east");
paper_rows = table(datetime("2026-10-21"), "reference", "night_us_20261021", "us_east", reference.times(1), ...
    reference.times(end), numel(reference.times), 30, round(0.86 * numel(reference.times) * numel(us.radars)), NaN, ...
    VariableNames=rain_rows.Properties.VariableNames);
paper_names = ["keywest" "midwest"];
for name = paper_names
    c = case_settings(name);
    step = minutes(c.times(2) - c.times(1));
    paper_rows = [paper_rows; table(dateshift(c.times(1) - hours(12), "start", "day"), "paper", name, name, ...
        c.times(1), c.times(end), numel(c.times), step, numel(c.radars) * numel(c.times), NaN, ...
        VariableNames=rain_rows.Properties.VariableNames)]; %#ok<AGROW>
end

%% The list and its budget
columns = ["date" "purpose" "name" "radar_set" "first_frame" "last_frame" "n_frames" "step_min" "volumes"];
rain_rows.uses = repmat("radar patterns", height(rain_rows), 1);
paper_rows.uses = repmat("paper rerun", height(paper_rows), 1);
paper_rows.uses(paper_rows.purpose == "reference") = "aspect near-radar prior terrain-across reference";
night_list = [selected(:, [columns "uses"]); rain_rows(:, [columns "uses"]); paper_rows(:, [columns "uses"])];
night_list.download_gb = night_list.volumes * volume_gb;
night_list.clean_cpu_h = night_list.volumes * cpu_s_per_volume / 3600;
detail = ["phase" "season" "track_deg" "axis_coherence" "speed_ms" "median_refl" "dry_share" ...
    "crossing" "appalachians_median_refl" "plains_median_refl"];
[~, i_sel] = ismember(night_list.name, selected.name);
for name = detail
    if isnumeric(selected.(name))
        night_list.(name) = nan(height(night_list), 1);
    else
        night_list.(name) = repmat("", height(night_list), 1);
    end
    night_list.(name)(i_sel > 0) = selected.(name)(i_sel(i_sel > 0));
end
disp(night_list(:, ["name" "phase" "track_deg" "speed_ms" "median_refl" "first_frame" "n_frames" "volumes" "uses"]))
fprintf("%d full nights, %d rain days, %d paper cases: %d volumes, %.0f GB to download, %.0f CPU-hours of cleaning\n", ...
    height(selected), height(rain_rows), height(paper_rows), sum(night_list.volumes), ...
    sum(night_list.download_gb), sum(night_list.clean_cpu_h))
writetable(night_list, fullfile(figure_dir, "NS_A_nights.csv"))

%% Plot results
% (a) Domain-mean track of every qualifying night against the day of
% year, the chosen nights by purpose; (b) speed against reflectivity.
fig = bmm4d.plot.create_export_figure(3000, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
terrain = contains(selected.uses, "terrain");
edges = contains(selected.uses, "edges");
groups = {~terrain & ~edges, "other studies", [0.20 0.20 0.20]; terrain & ~edges, "+ terrain", [0.65 0.40 0.15]; ...
    edges & ~terrain, "+ edges", [0.30 0.60 0.30]; terrain & edges, "+ terrain and edges", [0.45 0.30 0.60]};
panels = {"track_deg", "Domain-mean track (deg)"; "speed_ms", "Domain-mean speed (m/s)"};
for i_p = 1:2
    ax = nexttile(layout);
    hold(ax, "on");
    box(ax, "on")
    if i_p == 1
        x_all = day(pool.date, "dayofyear");
        x_sel = day(selected.date, "dayofyear");
        xlabel(ax, "Day of year")
    else
        x_all = pool.median_refl;
        x_sel = selected.median_refl;
        set(ax, XScale="log")
        xlabel(ax, "Median night-mean reflectivity (cm^2 km^{-2})")
    end
    scatter(ax, x_all, pool.(panels{i_p, 1}), 20, [0.75 0.75 0.75], "filled", DisplayName="qualifying nights")
    for i_g = 1:size(groups, 1)
        rows = groups{i_g, 1};
        scatter(ax, x_sel(rows), selected.(panels{i_p, 1})(rows), 90, groups{i_g, 3}, "filled", ...
            MarkerEdgeColor="w", DisplayName=groups{i_g, 2})
    end
    xlim(ax, "padded")
    ylabel(ax, panels{i_p, 2})
    if i_p == 1
        yticks(ax, 0:45:360)
        ylim(ax, [0 360])
        legend(ax, Location="northeast")
    end
    set(ax, FontSize=13)
end
title(layout, sprintf("Night set: %d full nights from %d qualifying (Dark Ecology 2021-2025)", height(selected), ...
    height(pool)), FontSize=14, FontName="Helvetica")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "NS_A_nights.png"), Resolution=300)
close(fig)

%% Local functions
function metrics = night_metrics(nights, date_id, rows, prefix)
% Per night, over the radars ROWS: reporting radars, dry share, median
% reflectivity, reflectivity-weighted track, its axis coherence and speed
n_date = max(date_id);
ok = rows & nights.fraction_missing < 0.1;
has_track = ok & isfinite(nights.u) & isfinite(nights.v) & isfinite(nights.mean_refl);
weight = max(nights.mean_refl, 0) .* has_track;
u = nights.u;
v = nights.v;
u(~has_track) = 0;
v(~has_track) = 0;
axis_vector = exp(2i * deg2rad(nights.direction));
axis_vector(~has_track | ~isfinite(axis_vector)) = 0;
n_radar = accumarray(date_id, ok, [n_date 1]);
refl = nights.mean_refl;
refl(~ok) = nan;
sum_weight = accumarray(date_id, weight, [n_date 1]);
metrics = table(n_radar, accumarray(date_id, ok & nights.fraction_rain < 0.1, [n_date 1]) ./ max(n_radar, 1), ...
    accumarray(date_id, refl, [n_date 1], @(x) median(x, "omitnan")), ...
    mod(atan2d(accumarray(date_id, weight .* u, [n_date 1]), accumarray(date_id, weight .* v, [n_date 1])), 360), ...
    abs(accumarray(date_id, weight .* axis_vector, [n_date 1])) ./ sum_weight, ...
    accumarray(date_id, weight .* hypot(u, v), [n_date 1]) ./ sum_weight, ...
    VariableNames=prefix + ["n_radar" "dry_share" "median_refl" "track_deg" "axis_coherence" "speed_ms"]);
end

function allowed = allowed_nights(pool, chosen)
% At least 4 days from every chosen night, at most 5 chosen nights per year
near = false(height(pool), 1);
for d = pool.date(chosen)'
    near = near | abs(days(pool.date - d)) < 4;
end
per_year = accumarray(year(pool.date(chosen)) - 2020, 1, [5 1]);
allowed = ~near & per_year(year(pool.date) - 2020) < 5;
end

function [pool, chosen] = add_until(pool, chosen, good, n_needed, purpose, rank_by)
% Add the strongest GOOD nights (by RANK_BY) until N_NEEDED chosen nights are GOOD
[~, order] = sort(pool.(rank_by), "descend");
for i = order'
    if nnz(chosen & good) >= n_needed
        break
    end
    if good(i) && ~chosen(i) && allowed_one(pool, chosen, i)
        chosen(i) = true;
        pool.purpose(i) = purpose;
    end
end
end

function ok = allowed_one(pool, chosen, i)
allowed = allowed_nights(pool, chosen);
ok = allowed(i);
end
