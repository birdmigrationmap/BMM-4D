%% Sandbox: raw reflectivity and estimator density by elevation near the radar
% Held out of the solve, the near rows of elevations >= 1.5 deg read low
% where gates hold few birds and those below 1.5 deg read high
% (SB_near_elevation). Is that in the radar data or in the gate estimator?
% The mean reflectivity of many gates does not depend on how many birds a
% gate holds (gates without echo count as zero), so per sweep and 100 m of
% height above the antenna, within 3-60 km, compare
%   raw:       mean reflectivity of the sweep, as density (/ 11 cm^2)
%   estimator: mean density of the cleaned volume (gate estimator)
% on the same gates (those kept by cleaning: no rain, no blocked beam). Each
% volume is cleaned again at full gate resolution with the production read
% settings except the rain filter (rhoHV instead of MistNet: both nights are
% dry, and the cached volumes are thinned to 2 km). Volumes of the two
% nights (night_20261021 and night_20210512, 28 radars), one per radar and
% hour.
% Adjacent sweeps of a volume at the same height differ in range (the higher
% one is closer) and may differ in waveform: the lowest cuts are long-PRT
% surveillance scans (PRF ~320 Hz, 0.5 deg radials), the upper ones short-PRT
% Doppler or batch scans (~925 Hz, 1 deg radials). The log ratio of a pair is
% fitted (robust least squares, all pairs where the lower sweep reads > 1
% bird km^-3) as f(r_upper) - f(r_lower) + w (short-PRT upper and long-PRT
% lower), with f a factor per range bin (mean 1 at 30-60 km) and w a factor
% per scan-strategy family, radar-night or hour.
% Outputs: figures/near_raw.png, near_raw.csv (one row per volume, sweep
% and height), near_raw_fit.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
nights = ["night_20261021" "night_20210512"];
height_edges = 0:100:3000;
range_limits_m = [3000 60000];
c = case_settings("sandbox");
read_options = [c.read_options {"removeRain", "rho", "mistnetScreen", false, "maxRange", 65000, "data_res", [0 0]}];

%% Volumes: one per radar and hour ------------------------------------------
volumes = table();
for night = nights
    c = case_settings(night);
    files = dir(fullfile(bmm4d.radar.volume_cache_dir(c), "K*.mat"));
    names = string({files.name})';
    radar = extractBefore(names, 5);
    time = datetime(extractBetween(names, 5, 19), InputFormat="yyyyMMdd_HHmmss");
    keep = ismember(radar, c.radars) & time >= c.times(1) - minutes(15) & time <= c.times(end) + minutes(15);
    hour_id = findgroups(radar(keep), dateshift(time(keep), "start", "hour"));
    [~, first] = unique(hour_id);
    kept = find(keep);
    volumes = [volumes; table(repmat(night, numel(first), 1), radar(kept(first)), time(kept(first)), ...
        string(fullfile({files(kept(first)).folder}, {files(kept(first)).name}))', ...
        VariableNames=["night" "radar" "time" "file"])]; %#ok<AGROW>
end
fprintf("%d volumes\n", height(volumes))

%% Sums per sweep and height ------------------------------------------------
% One file per volume, so an interrupted run resumes (about 10 s a volume)
part_dir = fullfile(fileparts(c.data_dir), "near_raw");
if ~isfolder(part_dir), mkdir(part_dir), end
part_file = fullfile(part_dir, volumes.radar + "_" + string(volumes.time, "yyyyMMdd_HHmmss") + ".mat");
todo = find(~isfile(part_file));
parfor i_t = 1:numel(todo)
    part = sweep_sums(volumes.file(todo(i_t)), read_options, height_edges, range_limits_m);
    parsave(part_file(todo(i_t)), part)
end
rows = table();
for i_v = 1:height(volumes)
    part = load(part_file(i_v)).part;
    if isempty(part), continue, end
    part.night = repmat(volumes.night(i_v), height(part), 1);
    part.radar = repmat(volumes.radar(i_v), height(part), 1);
    part.time = repmat(volumes.time(i_v), height(part), 1);
    rows = [rows; part]; %#ok<AGROW>
end
writetable(rows, fullfile(figure_dir, "near_raw.csv"))

%% Pairs of adjacent sweeps at the same height ---------------------------------
rows.e = round(rows.elev * 10) / 10;
rows = sortrows(rows, ["night" "radar" "time" "height_m" "e"]);
below = rows(1:end - 1, :);
above = rows(2:end, :);
same = below.radar == above.radar & below.time == above.time & below.height_m == above.height_m;
below = below(same, :);
above = above(same, :);
pairs = table(below.night, below.radar, below.time, below.e, above.e, below.range_m, above.range_m, ...
    below.raw_sum ./ below.n_gate, ...
    log2((above.raw_sum ./ above.n_gate) ./ (below.raw_sum ./ below.n_gate)), ...
    log2((above.estimator_sum ./ above.n_gate) ./ (below.estimator_sum ./ below.n_gate)), ...
    double(above.prf > 500) - double(below.prf > 500), ...
    VariableNames=["night" "radar" "time" "e_lower" "e_upper" "r_lower" "r_upper" "raw_lower" "log_raw" ...
    "log_estimator" "waveform_step"]);
pairs = pairs(pairs.raw_lower > 1 & isfinite(pairs.log_raw) & isfinite(pairs.log_estimator), :);
% Scan-strategy family: VCP 31/32 (0.5, 1.5, 2.5, 3.5, 4.5) or VCP 212/215/35 (0.5, 0.9, 1.3, 1.8, ...)
pairs.vcp32 = ismember(pairs.e_upper, [2.5 3.5 4.5]) | ismember(pairs.e_lower, [1.5 2.5 3.5]);
fprintf("%d pairs\n", height(pairs))

%% Range and waveform factors ----------------------------------------------------
range_edges_km = [3 5 7.5 10 15 20 30 45 60];
family = ["VCP 212/215/35" "VCP 31/32"];
[raw_range, raw_family] = fit_steps(pairs, "log_raw", 1 + pairs.vcp32, range_edges_km);
[estimator_range, estimator_family] = fit_steps(pairs, "log_estimator", 1 + pairs.vcp32, range_edges_km);
[radar_night, radar_night_name] = findgroups(pairs.night + " " + pairs.radar);
[~, radar_night_factor] = fit_steps(pairs, "log_estimator", radar_night, range_edges_km);
radar_night_pairs = accumarray(radar_night, abs(pairs.waveform_step));
radar_night_vcp32 = accumarray(radar_night, pairs.vcp32, [], @mean) > 0.5;
[hour_group, hour_night, hour_time] = findgroups(pairs.night, dateshift(pairs.time, "start", "hour"));
[~, hour_factor] = fit_steps(pairs, "log_estimator", hour_group, range_edges_km);
hour_pairs = accumarray(hour_group, abs(pairs.waveform_step));
range_label = string(range_edges_km(1:end - 1)) + "-" + string(range_edges_km(2:end)) + " km";
fit_table = table([range_label'; "short-PRT, " + family'], round([raw_range; raw_family], 2), ...
    round([estimator_range; estimator_family], 2), VariableNames=["term" "raw" "estimator"]);
disp(fit_table)
writetable(fit_table, fullfile(figure_dir, "near_raw_fit.csv"))

%% Plot results -----------------------------------------------------------
gd = readtable(fullfile(project_root, "studies", "gate_density", "figures", "GD_I_range_bias.csv"));
range_centre = (range_edges_km(1:end - 1) + range_edges_km(2:end)) / 2;
fig = bmm4d.plot.create_export_figure(3600, 1800, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "Adjacent sweeps at the same height, 571 volumes (28 radars, 21 Oct 2026 and 12 May 2021)", FontSize=15, Color="k")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
plot(ax, range_centre, raw_range, "-o", Color=[0.45 0.45 0.45], LineWidth=2, MarkerFaceColor=[0.45 0.45 0.45], ...
    DisplayName="raw reflectivity (pairs)")
plot(ax, range_centre, estimator_range, "-o", Color=[0 0.45 0.74], LineWidth=2, MarkerFaceColor=[0 0.45 0.74], ...
    DisplayName="estimator (pairs)")
plot(ax, gd.range_km, gd.All124Scans_radarOutput, "--", Color=[0.45 0.45 0.45], LineWidth=1.5, ...
    DisplayName="radar output, GD\_I (124 scans)")
plot(ax, gd.range_km, gd.All124Scans_estimator, "--", Color=[0 0.45 0.74], LineWidth=1.5, DisplayName="estimator, GD\_I")
yline(ax, 1, "k:", HandleVisibility="off")
xlim(ax, [0 60]); ylim(ax, [0 1.4]); set(ax, FontSize=13)
xlabel(ax, "Range (km)"); ylabel(ax, "Factor (1 at 30-60 km)")
legend(ax, Location="southeast", FontSize=12)
title(ax, "Range factor f(r)", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
show = radar_night_pairs >= 20;
[~, order] = sort(radar_night_factor);
order = order(show(order));
colors = [0.85 0.33 0.10; 0 0.45 0.74];
for i_f = 1:2
    in_family = order(radar_night_vcp32(order) == (i_f == 2));
    [~, position] = ismember(in_family, order);
    scatter(ax, radar_night_factor(in_family), position, 60, colors(i_f, :), "filled", DisplayName=family(i_f))
end
set(ax, YTick=1:numel(order), YTickLabel=extractAfter(radar_night_name(order), "night_"), FontSize=9, TickLabelInterpreter="none")
xline(ax, 1, "k:", HandleVisibility="off")
xlim(ax, [0.5 1.3]); ylim(ax, [0 numel(order) + 1])
xlabel(ax, "Short-PRT / long-PRT sweep (estimator)", FontSize=13)
legend(ax, Location="southeast", FontSize=12)
title(ax, "Waveform factor by radar and night", FontWeight="normal", FontSize=13)
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
night_colors = [0.47 0.67 0.19; 0.49 0.18 0.56];
for i_n = 1:numel(nights)
    h = hour_night == nights(i_n) & hour_pairs >= 30;
    plot(ax, hour(hour_time(h)) + 24 * (hour(hour_time(h)) < 12), hour_factor(h), "-o", Color=night_colors(i_n, :), ...
        LineWidth=2, MarkerFaceColor=night_colors(i_n, :), DisplayName=nights(i_n))
end
yline(ax, 1, "k:", HandleVisibility="off")
ylim(ax, [0.5 1.2]); set(ax, FontSize=13, XTick=24:2:34, XTickLabel=string(mod(24:2:34, 24)))
xlabel(ax, "Hour (UTC)"); ylabel(ax, "Short-PRT / long-PRT sweep (estimator)")
legend(ax, Location="southwest", FontSize=12, Interpreter="none")
title(ax, "Waveform factor by hour", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "near_raw.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function sums = sweep_sums(cache_file, read_options, height_edges, range_limits_m)
% Per sweep and height bin: gate count, gates with echo, sum of the raw
% density (no echo = 0) and of the estimator density, mean range, and the
% sweep's radials, PRF and gate length.
raw_file = load(cache_file).volume.source_file;
cleaned = bmm4d.radar.read_rd(raw_file, read_options{:});
cleaned = cleaned{1};
raw = bmm4d.radar.rsl2mat_wrap(raw_file, struct("max_elev", 7));
raw_sweeps = unique_elev_sweeps(raw, "dz");
sums = table();
for i_s = 1:numel(cleaned.dens.sweeps)
    sweep = cleaned.dens.sweeps(i_s);
    [~, i_raw] = min(abs([raw_sweeps.elev] - sweep.elev));
    [dbz, range_m, az] = sweep2mat(raw_sweeps(i_raw));
    if ~isequal(size(dbz), size(sweep.data)) || max(abs(az(:) - sweep.az(:))) > 0.01
        continue
    end
    raw_density = z_to_refl(idb(dbz)) / 11;
    raw_density(isnan(dbz)) = 0;
    use = ~isnan(sweep.data);
    [~, height_m] = slant2ground(range_m, sweep.elev);
    height_bin = discretize(height_m, height_edges);
    in_range = range_m >= range_limits_m(1) & range_m <= range_limits_m(2) & ~isnan(height_bin);
    range_all = repmat(range_m, 1, size(dbz, 2));
    for i_h = unique(height_bin(in_range))'
        g = use & in_range & height_bin == i_h;
        n = nnz(g);
        if n < 200, continue, end
        sums = [sums; table(sweep.elev, height_edges(i_h) + 50, n, nnz(g & ~isnan(dbz)), sum(raw_density(g)), ...
            sum(sweep.data(g)), mean(range_all(g)), numel(az), raw_sweeps(i_raw).prf, median(diff(range_m)), ...
            VariableNames=["elev" "height_m" "n_gate" "n_echo" "raw_sum" "estimator_sum" "range_m" "n_az" "prf" ...
            "gate_m"])]; %#ok<AGROW>
    end
end
end

function [range_factor, group_factor] = fit_steps(pairs, column, group, range_edges_km)
% Robust fit (Huber, 6 passes) of the pairs' log2 ratio COLUMN as
% f(r_upper) - f(r_lower) + w(group) * waveform_step; f has mean 0 at
% 30-60 km. Returns 2^f per range bin and 2^w per group.
n = height(pairs);
n_bin = numel(range_edges_km) - 1;
n_group = max(group);
upper_bin = discretize(pairs.r_upper / 1000, range_edges_km);
lower_bin = discretize(pairs.r_lower / 1000, range_edges_km);
A = [sparse([1:n 1:n]', [upper_bin; lower_bin], [ones(n, 1); -ones(n, 1)], n, n_bin), ...
    sparse(1:n, group, pairs.waveform_step, n, n_group)];
anchor = zeros(1, n_bin + n_group);
anchor(range_edges_km(1:end - 1) >= 30) = 1 / nnz(range_edges_km(1:end - 1) >= 30);
y = pairs.(column);
weight = ones(n, 1);
for i_pass = 1:6
    x = [A .* sqrt(weight); 100 * anchor] \ [y .* sqrt(weight); 0];
    residual = y - A * x;
    scale = 1.4826 * median(abs(residual));
    weight = min(1, 1.345 * scale ./ abs(residual));
end
range_factor = 2 .^ x(1:n_bin);
group_factor = 2 .^ x(n_bin + 1:end);
end

function parsave(file, part)
save(file, "part")
end
