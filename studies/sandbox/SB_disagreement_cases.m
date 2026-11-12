%% Sandbox: radar-nights where the sweeps disagree near the radar
% KUEX on 21 Oct 2026 (SB_kuex) showed a low layer moving north under the
% southbound migrants, which the low sweeps sample over 0-30 km and the
% upper sweeps barely cross. Is that the usual story? The same diagnostics
% on eight radar-nights with a strongly negative score in SB_obs_sections
% (upper sweeps within 20 km read low) and two with a score >= 0:
%   1. motion profile of the volume (EP-VVP, direction moved to, speed);
%   2. cleaned density and 3. radial velocity of the lowest sweep, 30 km;
%   4. beam consistency: azimuth-mean density of the sweeps < 1.7 deg
%      (solid) against the profile of the narrow beams (>= 1.8 deg within
%      25 km) passed through their beam (0.925 deg, as the forward model;
%      dashed);
%   5. that narrow-beam density profile;
%   6. mean differential reflectivity ZDR by height (gates with echo within
%      30 km, all sweeps): insects are usually well above birds.
% One volume per radar-night (the one nearest 04 UTC for the full nights).
% Outputs: figures/disagreement_cases.png, disagreement_cases.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
cache_dir = bmm4d.radar.volume_cache_dir(c);
read_options = [c.read_options {"maxRange", 40000, "data_res", [0 0]}];
hpbw = 0.925;
cases = table( ...
    ["KGRB" "KABR" "KVNX" "KHPX" "KIWX" "KGRR" "KDDC" "KLOT" "KSGF" "KDMX"]', ...
    ["2023-05-20" "2025-09-27" "2021-10-21" "2023-05-20" "2023-05-16" "2021-05-11" "2021-10-21" "2026-10-21" ...
    "2025-10-07" "2023-05-16"]', VariableNames=["radar" "night"]);
scores = readtable(fullfile(figure_dir, "obs_sections_score.csv"), ReadRowNames=true, ReadVariableNames=true, VariableNamingRule="preserve");
cases.score = arrayfun(@(i) scores{cases.radar(i), cases.night(i)}, (1:height(cases))');

%% Read each case ---------------------------------------------------------------
profile_edges = (0:50:2500)';
profile_centre = (profile_edges(1:end - 1) + profile_edges(2:end)) / 2;
zdr_edges = (0:100:2500)';
data = cell(height(cases), 1);
for i_c = 1:height(cases)
    files = dir(fullfile(cache_dir, cases.radar(i_c) + "*.mat"));
    time = datetime(extractBetween(string({files.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
    night = string(dateshift(time - hours(12), "start", "day"), "yyyy-MM-dd");
    in_night = find(night == cases.night(i_c));
    [~, k] = min(abs(time(in_night) - (datetime(cases.night(i_c)) + days(1) + hours(4))));
    file = fullfile(files(in_night(k)).folder, files(in_night(k)).name);
    data{i_c} = case_data(file, read_options);
    data{i_c}.time = time(in_night(k));
    fprintf("%s %s %s\n", cases.radar(i_c), cases.night(i_c), string(data{i_c}.time))
end

%% Summary per case ----------------------------------------------------------------
for i_c = 1:height(cases)
    d = data{i_c};
    z = (50:100:2950)';
    uv = d.uv;
    direction = mod(atan2d(uv(:, 1), uv(:, 2)), 360);
    low = z < 300; aloft = z >= 1000 & z < 2000;
    cases.direction_low(i_c) = round(mod(atan2d(mean(uv(low, 1), "omitnan"), mean(uv(low, 2), "omitnan")), 360));
    cases.direction_aloft(i_c) = round(mod(atan2d(mean(uv(aloft, 1), "omitnan"), mean(uv(aloft, 2), "omitnan")), 360));
    cases.turning(i_c) = abs(mod(cases.direction_aloft(i_c) - cases.direction_low(i_c) + 180, 360) - 180);
    % Narrow-beam profile and beam test for the sweeps below 1.7 deg within 25 km
    [narrow_profile, measured, predicted] = beam_test(d.sweeps, profile_edges, hpbw);
    d.narrow_profile = narrow_profile; d.measured = measured; d.predicted = predicted;
    ratio = cellfun(@(m, p) sum(m(1:45), "omitnan") / sum(p(1:45), "omitnan"), measured, predicted);   % 3-25 km
    cases.low_sweeps_over_narrow(i_c) = round(median(ratio, "omitnan"), 2);
    % ZDR by height, and below 300 m against 1-2 km
    [zdr_profile, zdr_low, zdr_aloft] = zdr_by_height(d.sweeps, zdr_edges);
    d.zdr_profile = zdr_profile;
    cases.zdr_low(i_c) = round(zdr_low, 1);
    cases.zdr_aloft(i_c) = round(zdr_aloft, 1);
    data{i_c} = d;
end
disp(cases)
writetable(cases, fullfile(figure_dir, "disagreement_cases.csv"))

%% Plot results -----------------------------------------------------------
column_titles = ["Motion (direction moved to)" "Lowest sweep: cleaned density" "Lowest sweep: radial velocity" ...
    "Low sweeps: measured (solid) / from narrow beams (dashed)" "Narrow-beam density profile" "ZDR by height"];
fig = bmm4d.plot.create_export_figure(3600, 620 * height(cases), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, height(cases), 6, TileSpacing="compact", Padding="tight");
title(t, "Radar-nights with sweeps disagreeing near the radar (first eight) and two that agree (last two)", FontSize=17, Color="k")
sweep_colors = lines(8);
for i_c = 1:height(cases)
    d = data{i_c};
    z = (50:100:2950)';
    % 1. Motion profile: direction as dots, size by speed
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    direction = mod(atan2d(d.uv(:, 1), d.uv(:, 2)), 360);
    speed = hypot(d.uv(:, 1), d.uv(:, 2));
    ok = isfinite(direction);
    scatter(ax, direction(ok), z(ok), 8 + 6 * speed(ok), speed(ok), "filled")
    colormap(ax, copper(256)); clim(ax, [0 15])
    xlim(ax, [0 360]); ylim(ax, [0 2500]); set(ax, XTick=0:90:360, XTickLabel=["N" "E" "S" "W" "N"], FontSize=11)
    ylabel(ax, sprintf("%s %s\n%s UTC, score %.2f", cases.radar(i_c), cases.night(i_c), string(d.time, "HH:mm"), ...
        cases.score(i_c)), FontWeight="bold", FontSize=12)
    if i_c == 1, title(ax, column_titles(1), FontSize=12), end
    % 2-3. Maps of the lowest sweep (above 0.3 deg)
    lowest = find([d.sweeps.elev] > 0.3, 1);
    sweep = d.sweeps(lowest);
    [az_grid, range_grid] = meshgrid(sweep.az, sweep.range_m / 1000);
    x = range_grid .* sind(az_grid); y = range_grid .* cosd(az_grid);
    for i_q = 1:2
        ax = nexttile(t);
        if i_q == 1
            surface(ax, x, y, zeros(size(x)), log10(max(sweep.density, 0.1)), EdgeColor="none")
            colormap(ax, parula); clim(ax, [-1 2.5])
        else
            surface(ax, x, y, zeros(size(x)), sweep.vel, EdgeColor="none")
            colormap(ax, redblue_map()); clim(ax, [-20 20])
        end
        view(ax, 2); axis(ax, "equal"); xlim(ax, [-30 30]); ylim(ax, [-30 30]); set(ax, Color=[0.9 0.9 0.9], FontSize=10)
        hold(ax, "on")
        for r = 10:10:30, plot(ax, r * sind(0:360), r * cosd(0:360), "k:", LineWidth=0.5), end
        text(ax, -29, 28, sprintf("%.1f°", sweep.elev), FontSize=11, FontWeight="bold")
        if i_c == 1, title(ax, column_titles(1 + i_q), FontSize=12), end
    end
    % 4. Beam test
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    range_km = (3:0.5:40)';
    for i_l = 1:numel(d.measured)
        plot(ax, range_km, d.measured{i_l}, "-", Color=sweep_colors(i_l, :), LineWidth=2)
        plot(ax, range_km, d.predicted{i_l}, "--", Color=sweep_colors(i_l, :), LineWidth=2)
    end
    set(ax, YScale="log", FontSize=11); xlim(ax, [0 40]); ylim(ax, [0.3 1000])
    text(ax, 39, 0.5, sprintf("measured / narrow, 3-25 km: %.2f", cases.low_sweeps_over_narrow(i_c)), ...
        HorizontalAlignment="right", FontSize=10, BackgroundColor="w")
    if i_c == 1, title(ax, column_titles(4), FontSize=12), end
    % 5. Narrow-beam profile
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, d.narrow_profile, profile_centre, "k-", LineWidth=2)
    set(ax, XScale="log", FontSize=11); xlim(ax, [0.3 1000]); ylim(ax, [0 2500])
    if i_c == 1, title(ax, column_titles(5), FontSize=12), end
    % 6. ZDR
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    plot(ax, d.zdr_profile, (zdr_edges(1:end - 1) + zdr_edges(2:end)) / 2, "-", Color=[0.85 0.33 0.1], LineWidth=2)
    xlim(ax, [-2 10]); ylim(ax, [0 2500]); set(ax, FontSize=11)
    if i_c == 1, title(ax, column_titles(6), FontSize=12), end
    if i_c == height(cases)
        xlabel(nexttile(t, 6 * (i_c - 1) + 4), "Range (km)")
        xlabel(nexttile(t, 6 * (i_c - 1) + 5), "Birds km^{-3}")
        xlabel(ax, "ZDR (dB)")
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "disagreement_cases.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function d = case_data(cache_file, read_options)
% Raw moments (REF, VEL, ZDR, CFP) and the cleaned density at full gates,
% one pass per elevation (the surveillance pass of split cuts, with the
% velocity of its Doppler pass), within 40 km; and the volume's motion profile.
volume = load(cache_file).volume;
moments = bmm4d.radar.read_level2_moments(volume.source_file, ["REF" "VEL" "ZDR" "CFP"], 7);
cleaned = bmm4d.radar.read_rd(volume.source_file, read_options{:});
cleaned = cleaned{1};
sweeps = struct("elev", {}, "az", {}, "range_m", {}, "vel", {}, "zdr", {}, "density", {});
for e = unique(round([moments.elev] * 10) / 10)
    passes = moments(abs([moments.elev] - e) < 0.1 & arrayfun(@(m) ~isempty(m.REF), moments));
    if isempty(passes), continue, end
    [~, k] = max(arrayfun(@(m) numel(m.REF_range), passes));
    main = passes(k);
    n_gate = find(main.REF_range <= 40000, 1, "last");
    vel = nan(n_gate, numel(main.az));
    with_vel = passes(arrayfun(@(m) ~isempty(m.VEL), passes));
    if ~isempty(with_vel) && numel(with_vel(1).az) == numel(main.az)
        n_vel = min(n_gate, numel(with_vel(1).VEL_range));
        vel(1:n_vel, :) = with_vel(1).VEL(1:n_vel, :);
    end
    zdr = nan(n_gate, numel(main.az));
    if ~isempty(main.ZDR)
        n_zdr = min(n_gate, numel(main.ZDR_range));
        zdr(1:n_zdr, :) = main.ZDR(1:n_zdr, :);
    end
    [~, i_clean] = min(abs([cleaned.dens.sweeps.elev] - main.elev));
    clean_sweep = cleaned.dens.sweeps(i_clean);
    [~, az_id] = min(abs(mod(main.az(:) - clean_sweep.az(:)' + 180, 360) - 180), [], 2);
    range_id = interp1(clean_sweep.range(:), 1:numel(clean_sweep.range), main.REF_range(1:n_gate), "nearest", NaN);
    density = nan(n_gate, numel(main.az));
    ok = ~isnan(range_id);
    density(ok, :) = clean_sweep.data(range_id(ok), az_id);
    sweeps(end + 1) = struct("elev", main.elev, "az", main.az, "range_m", main.REF_range(1:n_gate), ... %#ok<AGROW>
        "vel", vel, "zdr", zdr, "density", density);
end
d = struct("sweeps", sweeps, "uv", volume.uv);
end

function [narrow_profile, measured, predicted] = beam_test(sweeps, profile_edges, hpbw)
% Profile from the sweeps >= 1.8 deg within 3-25 km (50 m bins); each sweep
% below 1.7 deg (above 0.3 deg): its azimuth-mean density by 0.5 km range
% and that profile through its Gaussian beam.
profile_centre = (profile_edges(1:end - 1) + profile_edges(2:end)) / 2;
heights = []; values = [];
for sweep = sweeps
    if sweep.elev < 1.7, continue, end
    near = sweep.range_m <= 25000 & sweep.range_m >= 3000;
    [~, height_m] = slant2ground(sweep.range_m(near), sweep.elev);
    heights = [heights; height_m]; %#ok<AGROW>
    values = [values; mean(sweep.density(near, :), 2, "omitnan")]; %#ok<AGROW>
end
bin = discretize(heights, profile_edges);
ok = ~isnan(bin) & ~isnan(values);
narrow_profile = accumarray(bin(ok), values(ok), [numel(profile_centre) 1], @mean, NaN);
filled = fillmissing(narrow_profile, "linear", EndValues="nearest");
range_m = (3000:500:40000)';
low = find([sweeps.elev] > 0.3 & [sweeps.elev] < 1.7);
measured = cell(1, numel(low)); predicted = cell(1, numel(low));
for i_l = 1:numel(low)
    sweep = sweeps(low(i_l));
    [~, centre_m] = slant2ground(range_m, sweep.elev);
    sd = range_m * deg2rad(hpbw) / (2 * sqrt(2 * log(2)));
    weight = exp(-(profile_centre' - centre_m) .^ 2 ./ (2 * sd .^ 2));
    predicted{i_l} = (weight * filled) ./ sum(weight, 2);
    range_bin = discretize(sweep.range_m, [range_m - 250; range_m(end) + 250]);
    ok = ~isnan(range_bin);
    measured{i_l} = accumarray(range_bin(ok), mean(sweep.density(ok, :), 2, "omitnan"), [numel(range_m) 1], @mean, NaN);
end
end

function [zdr_profile, zdr_low, zdr_aloft] = zdr_by_height(sweeps, zdr_edges)
% Mean ZDR (dB) of the gates with echo within 3-30 km, by height, all sweeps
heights = []; values = [];
for sweep = sweeps
    near = sweep.range_m >= 3000 & sweep.range_m <= 30000;
    [~, height_m] = slant2ground(sweep.range_m(near), sweep.elev);
    zdr = sweep.zdr(near, :);
    has = isfinite(zdr) & isfinite(sweep.density(near, :));
    [row, ~] = find(has);
    heights = [heights; height_m(row)]; %#ok<AGROW>
    values = [values; zdr(has)]; %#ok<AGROW>
end
bin = discretize(heights, zdr_edges);
ok = ~isnan(bin);
zdr_profile = accumarray(bin(ok), values(ok), [numel(zdr_edges) - 1 1], @median, NaN);
zdr_low = median(values(heights < 300), "omitnan");
zdr_aloft = median(values(heights >= 1000 & heights < 2000), "omitnan");
end

function map = redblue_map()
% Diverging blue-white-red for radial velocity
n = 128;
map = [linspace(0.2, 1, n)' linspace(0.3, 1, n)' ones(n, 1); ones(n, 1) linspace(1, 0.25, n)' linspace(1, 0.2, n)'];
end
