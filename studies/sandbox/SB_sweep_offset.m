%% Sandbox: why sweeps of one volume disagree at the same height
% At the same height the upper sweeps read lower than the lowest ones
% (SB_near_raw): a range part (the upper sweep is closer) and a part on
% the short-PRT sweeps of VCP 212/215/35 (x0.84). Candidate causes, each
% with its own test:
%   - waveform (calibration, noise subtraction, few pulses): the lowest
%     elevations are scanned twice, a long-PRT surveillance pass (CS) and a
%     short-PRT Doppler pass (CD) that also records reflectivity. Same
%     beam, same gates, ~30 s apart: compare them gate by gate.
%   - clutter filter: it removes echo near zero radial velocity, i.e. where
%     the beam is perpendicular to the birds' flight. A filter loss depends
%     on the angle between beam and flight direction; a calibration offset
%     does not. Also the share of gates the filter touched (CFP >= 3 dB).
%   - noise and detection: censoring of weak echoes cuts the low tail of
%     the dBZ distribution; an offset shifts the whole distribution.
% The sweep types of VCP 212: 0.5/0.9/1.3 deg split (CS for reflectivity,
% CD for velocity), 1.8/2.4 deg batch (reflectivity from long-PRT pulses,
% velocity from short-PRT pulses, both in one rotation), 3.1 deg and up CD
% (everything from short-PRT pulses).
% Volumes: VCP 212/215/35 volumes of the two nights (night_20261021,
% night_20210512), 02-07 UTC, one per radar and hour.
% Outputs: figures/sweep_offset_example.png, sweep_offset_split.png,
% sweep_offset_azimuth.png, sweep_offset.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
nights = ["night_20261021" "night_20210512"];
range_edges_m = 3000:1000:60000;
sector_edges = -180:15:180;              % angle between beam and flight direction (deg)
dbz_edges = -35:1:45;
difference_edges = -15.25:0.5:15.25;          % REF is coded in 0.5 dB steps

%% Volumes ------------------------------------------------------------------
% VCP 212/215/35 volumes have a 0.9 deg sweep; VCP 31/32 do not (near_raw.csv)
sweeps_table = readtable(fullfile(figure_dir, "near_raw.csv"), TextType="string");
sweeps_table.e = round(sweeps_table.elev * 10) / 10;
volume_id = findgroups(sweeps_table.radar, sweeps_table.time);
has_09 = splitapply(@(e) any(e == 0.9), sweeps_table.e, volume_id);
[~, first] = unique(volume_id);
volumes = sweeps_table(first, ["night" "radar" "time"]);
volumes = volumes(has_09 & hour(volumes.time) >= 2 & hour(volumes.time) <= 7, :);
c = case_settings("sandbox");
cache_dir = bmm4d.radar.volume_cache_dir(c);
volumes.file = fullfile(cache_dir, volumes.radar + string(volumes.time, "yyyyMMdd_HHmmss") + "_V06.mat");
volumes = volumes(isfile(volumes.file), :);
fprintf("%d volumes\n", height(volumes))

%% Statistics per volume ------------------------------------------------------
stats = cell(height(volumes), 1);
parfor i_v = 1:height(volumes)
    stats{i_v} = volume_stats(volumes.file(i_v), range_edges_m, sector_edges, dbz_edges, difference_edges);
end
volumes.stats = stats;

%% Same beam, two waveforms: CD / CS on the split sweeps ------------------------
% Mean of linear Z on the gates valid in both passes (no range folding),
% below-threshold gates as zero, by range; and the gate-by-gate dB difference
split_elevations = [0.5 0.9 1.3];
split_z = zeros(numel(range_edges_m) - 1, 2, numel(split_elevations));
split_difference = zeros(numel(difference_edges) - 1, numel(split_elevations));
for i_v = 1:height(volumes)
    s = volumes.stats{i_v};
    for i_e = 1:numel(split_elevations)
        k = find(abs(s.split_elev - split_elevations(i_e)) < 0.15, 1);
        if isempty(k), continue, end
        split_z(:, :, i_e) = split_z(:, :, i_e) + s.split_z(:, :, k);
        split_difference(:, i_e) = split_difference(:, i_e) + s.split_difference(:, k);
    end
end
range_centre_km = (range_edges_m(1:end - 1) + range_edges_m(2:end)) / 2000;
split_ratio = squeeze(split_z(:, 2, :) ./ split_z(:, 1, :));
difference_centre = (difference_edges(1:end - 1) + difference_edges(2:end)) / 2;
for i_e = 1:numel(split_elevations)
    cdf = cumsum(split_difference(:, i_e)) / sum(split_difference(:, i_e));
    fprintf("%.1f deg CD - CS: median %.2f dB, mean linear ratio 3-60 km %.3f\n", split_elevations(i_e), ...
        difference_centre(find(cdf >= 0.5, 1)), sum(split_z(:, 2, i_e)) / sum(split_z(:, 1, i_e)))
end

%% Sweep types against the angle to the flight direction --------------------------
% Per sweep type, within 3-30 km: mean Z by sector, each volume and range bin
% normalized by its own mean over sectors (so range and height cancel), and
% the share of gates with CFP >= 3 dB
types = ["CS" "batch" "CD"];
type_elevations = {[0.5 0.9 1.3], [1.8 2.4], [3.1 4.0 5.1 6.4]};
near = range_edges_m(1:end - 1) < 30000;
n_sector = numel(sector_edges) - 1;
sector_shape = zeros(n_sector, numel(types));
sector_weight = zeros(n_sector, numel(types));
cfp_share = zeros(n_sector, numel(types));
cfp_gates = zeros(n_sector, numel(types));
for i_v = 1:height(volumes)
    s = volumes.stats{i_v};
    for i_s = 1:numel(s.sweeps)
        sweep = s.sweeps(i_s);
        i_type = find(types == sweep.type);
        if isempty(i_type) || ~any(abs(type_elevations{i_type} - sweep.elev) < 0.15), continue, end
        mean_z = sweep.sum_z(near, :) ./ max(sweep.n_valid(near, :), 1);
        ring_mean = sum(sweep.sum_z(near, :), 2) ./ max(sum(sweep.n_valid(near, :), 2), 1);
        use = ring_mean > 0 & all(sweep.n_valid(near, :) > 5, 2);
        sector_shape(:, i_type) = sector_shape(:, i_type) + sum(mean_z(use, :) ./ ring_mean(use), 1)';
        sector_weight(:, i_type) = sector_weight(:, i_type) + nnz(use);
        cfp_share(:, i_type) = cfp_share(:, i_type) + sum(sweep.n_cfp(near, :), 1)';
        cfp_gates(:, i_type) = cfp_gates(:, i_type) + sum(sweep.n_cfp_known(near, :), 1)';
    end
end
sector_shape = sector_shape ./ sector_weight;
cfp_share = cfp_share ./ cfp_gates;

%% dBZ distributions at the same height: 1.3 deg CS against 1.8 deg batch ----------
% Echo gates in the range bins where both sweeps are at 600-900 m above the antenna
dbz_cs = zeros(numel(dbz_edges) - 1, 1);
dbz_batch = zeros(numel(dbz_edges) - 1, 1);
n_cs = 0; n_batch = 0;
for i_v = 1:height(volumes)
    s = volumes.stats{i_v};
    for i_s = 1:numel(s.sweeps)
        sweep = s.sweeps(i_s);
        [~, height_m] = slant2ground(range_centre_km' * 1000, sweep.elev);
        band = height_m >= 600 & height_m <= 900;
        if sweep.type == "CS" && abs(sweep.elev - 1.3) < 0.15
            dbz_cs = dbz_cs + sum(sweep.dbz_histogram(:, band), 2);
            n_cs = n_cs + sum(sweep.n_valid(band, :), "all");
        elseif sweep.type == "batch" && abs(sweep.elev - 1.8) < 0.15
            dbz_batch = dbz_batch + sum(sweep.dbz_histogram(:, band), 2);
            n_batch = n_batch + sum(sweep.n_valid(band, :), "all");
        end
    end
end
band_range_cs = range_centre_km(abs(slant2ground_height(range_centre_km, 1.3) - 750) <= 150);
band_range_batch = range_centre_km(abs(slant2ground_height(range_centre_km, 1.8) - 750) <= 150);

%% Save the summary ---------------------------------------------------------------
summary = table(range_centre_km', split_ratio(:, 1), split_ratio(:, 2), split_ratio(:, 3), ...
    VariableNames=["range_km" "cd_over_cs_05" "cd_over_cs_09" "cd_over_cs_13"]);
writetable(summary, fullfile(figure_dir, "sweep_offset.csv"))

%% Plot: one volume, every sweep ------------------------------------------------------
% KUEX 22 Oct 2026 04 UTC (waveform factor 0.65 in SB_near_raw): mean
% density against height for every sweep and pass, 3-60 km
example = volumes(volumes.radar == "KUEX" & hour(volumes.time) == 4 & volumes.night == "night_20261021", :);
s = example.stats{1};
fig = bmm4d.plot.create_export_figure(3600, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "KUEX " + string(example.time, "dd MMM yyyy HH:mm") + " UTC: every sweep of one volume, 3-60 km", FontSize=16, Color="k")
type_colors = struct("CS", [0.85 0.33 0.10], "CD_split", [0.93 0.69 0.13], "batch", [0 0.45 0.74], "CD", [0.49 0.18 0.56]);
type_names = struct("CS", "long-PRT surveillance (split)", "CD_split", "short-PRT Doppler pass (split)", ...
    "batch", "batch (1.8, 2.4°)", "CD", "short-PRT Doppler (3.1° up)");
for i_panel = 1:2
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    shown = strings(0);
    for i_s = 1:numel(s.sweeps)
        sweep = s.sweeps(i_s);
        [~, height_m] = slant2ground(range_centre_km' * 1000, sweep.elev);
        density = sum(sweep.sum_z, 2) ./ max(sum(sweep.n_valid, 2), 1);
        density = z_to_refl(density) / 11;
        ok = sum(sweep.n_valid, 2) > 100;
        if i_panel == 1
            x = density(ok); y = height_m(ok);
        else
            x = density(ok); y = range_centre_km(ok);
        end
        name = sweep.type;
        handle_visibility = "off";
        if ~ismember(name, shown), handle_visibility = "on"; shown(end + 1) = name; end %#ok<SAGROW>
        line_style = "-";
        if sweep.type == "CD_split", line_style = "--"; end
        plot(ax, x, y, line_style, Color=type_colors.(name), LineWidth=2, DisplayName=type_names.(name), ...
            HandleVisibility=handle_visibility)
        text(ax, x(end), y(end), sprintf(" %.1f°", sweep.elev), Color=type_colors.(name), FontSize=11)
    end
    set(ax, XScale="log", FontSize=13)
    xlabel(ax, "Mean raw density (birds km^{-3})")
    if i_panel == 1
        ylabel(ax, "Height above antenna (m)"); ylim(ax, [0 3000])
        title(ax, "Against height", FontWeight="normal")
        legend(ax, Location="northeast", FontSize=12)
    else
        ylabel(ax, "Range (km)"); ylim(ax, [0 60])
        title(ax, "Against range", FontWeight="normal")
    end
end
% Third panel: the same, but each sweep's density divided by the 0.5-1.3 deg
% CS profile interpolated to its heights
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
cs = s.sweeps([s.sweeps.type] == "CS");
reference_height = []; reference_density = [];
for i_s = 1:numel(cs)
    [~, height_m] = slant2ground(range_centre_km' * 1000, cs(i_s).elev);
    density = sum(cs(i_s).sum_z, 2) ./ max(sum(cs(i_s).n_valid, 2), 1);
    ok = sum(cs(i_s).n_valid, 2) > 100 & density > 0;
    reference_height = [reference_height; height_m(ok)]; %#ok<AGROW>
    reference_density = [reference_density; density(ok)]; %#ok<AGROW>
end
height_grid = (50:100:2950)';
reference_profile = accumarray(discretize(reference_height, 0:100:3000), reference_density, [30 1], @mean, NaN);
for i_s = 1:numel(s.sweeps)
    sweep = s.sweeps(i_s);
    if sweep.type == "CS", continue, end
    [~, height_m] = slant2ground(range_centre_km' * 1000, sweep.elev);
    density = sum(sweep.sum_z, 2) ./ max(sum(sweep.n_valid, 2), 1);
    ok = sum(sweep.n_valid, 2) > 100 & height_m < 3000;
    reference = interp1(height_grid, reference_profile, height_m(ok));
    line_style = "-";
    if sweep.type == "CD_split", line_style = "--"; end
    plot(ax, range_centre_km(ok), density(ok) ./ reference, line_style, Color=type_colors.(sweep.type), LineWidth=2)
    last = find(isfinite(density(ok) ./ reference), 1, "last");
    r_ok = range_centre_km(ok); ratio_ok = density(ok) ./ reference;
    if ~isempty(last), text(ax, r_ok(last), ratio_ok(last), sprintf(" %.1f°", sweep.elev), Color=type_colors.(sweep.type), FontSize=11), end
end
yline(ax, 1, "k:")
set(ax, YScale="log", FontSize=13); ylim(ax, [0.1 3]); xlim(ax, [0 60])
xlabel(ax, "Range (km)"); ylabel(ax, "Sweep / long-PRT sweeps at the same height")
title(ax, "Ratio to the 0.5-1.3° long-PRT sweeps", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sweep_offset_example.png"), Resolution=150)
close(fig)

%% Plot: same beam, two waveforms ------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 1150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, sprintf("Same beam, two waveforms: the two passes of the split sweeps (%d volumes, VCP 212/215/35)", height(volumes)), ...
    FontSize=16, Color="k")
split_colors = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.47 0.67 0.19];
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_e = 1:numel(split_elevations)
    plot(ax, range_centre_km, split_ratio(:, i_e), "-", Color=split_colors(i_e, :), LineWidth=2, ...
        DisplayName=sprintf("%.1f°", split_elevations(i_e)))
end
yline(ax, 1, "k:", HandleVisibility="off"); yline(ax, 0.84, "--", Color=[0 0.45 0.74], LineWidth=1.5, ...
    DisplayName="upper sweeps vs CS (SB\_near\_raw)")
ylim(ax, [0.5 1.3]); xlim(ax, [0 60]); set(ax, FontSize=13)
xlabel(ax, "Range (km)"); ylabel(ax, "Mean Z, Doppler pass / surveillance pass")
legend(ax, Location="southeast", FontSize=12)
title(ax, "Mean reflectivity ratio on the same gates", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_e = 1:numel(split_elevations)
    plot(ax, difference_centre, split_difference(:, i_e) / sum(split_difference(:, i_e)) / 0.5, "-", ...
        Color=split_colors(i_e, :), LineWidth=2, DisplayName=sprintf("%.1f°", split_elevations(i_e)))
end
xline(ax, 0, "k:", HandleVisibility="off"); xline(ax, 10 * log10(0.84), "--", Color=[0 0.45 0.74], LineWidth=1.5, ...
    DisplayName="x0.84 = -0.76 dB")
xlim(ax, [-10 10]); set(ax, FontSize=13)
xlabel(ax, "Doppler - surveillance, same gate (dB)"); ylabel(ax, "Density of gates")
legend(ax, Location="northwest", FontSize=12)
title(ax, "Gate by gate (both with echo, 3-60 km)", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
dbz_centre = (dbz_edges(1:end - 1) + dbz_edges(2:end)) / 2;
stairs(ax, dbz_edges(1:end - 1), dbz_cs / n_cs, Color=[0.85 0.33 0.10], LineWidth=2, ...
    DisplayName=sprintf("1.3° surveillance, %.0f-%.0f km", min(band_range_cs), max(band_range_cs)))
stairs(ax, dbz_edges(1:end - 1), dbz_batch / n_batch, Color=[0 0.45 0.74], LineWidth=2, ...
    DisplayName=sprintf("1.8° batch, %.0f-%.0f km", min(band_range_batch), max(band_range_batch)))
set(ax, FontSize=13); xlim(ax, [-30 40])
xlabel(ax, "Gate reflectivity (dBZ)"); ylabel(ax, "Share of gates (per dB)")
legend(ax, Location="northeast", FontSize=12)
title(ax, "Gate values at 600-900 m height", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sweep_offset_split.png"), Resolution=150)
close(fig)

%% Plot: angle to the flight direction --------------------------------------------------
sector_centre = (sector_edges(1:end - 1) + sector_edges(2:end)) / 2;
type_line_colors = [0.85 0.33 0.10; 0 0.45 0.74; 0.49 0.18 0.56];
fig = bmm4d.plot.create_export_figure(3600, 1150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
title(t, "Against the angle between beam and flight direction (3-30 km); birds' radial velocity is ~0 at ±90°", ...
    FontSize=16, Color="k")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_type = 1:numel(types)
    plot(ax, sector_centre, sector_shape(:, i_type), "-o", Color=type_line_colors(i_type, :), LineWidth=2, ...
        MarkerFaceColor=type_line_colors(i_type, :), DisplayName=types(i_type) + " (" + strjoin(string(type_elevations{i_type}), ", ") + "°)")
end
xline(ax, [-90 90], "k--", HandleVisibility="off"); yline(ax, 1, "k:", HandleVisibility="off")
set(ax, FontSize=13, XTick=-180:45:180); xlim(ax, [-180 180])
xlabel(ax, "Beam azimuth - flight heading (°)"); ylabel(ax, "Mean Z / mean over the ring")
legend(ax, Location="south", FontSize=12)
title(ax, "Reflectivity around the ring", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_type = [1 2 3]
    if all(isnan(cfp_share(:, i_type))), continue, end
    plot(ax, sector_centre, cfp_share(:, i_type), "-o", Color=type_line_colors(i_type, :), LineWidth=2, ...
        MarkerFaceColor=type_line_colors(i_type, :), DisplayName=types(i_type))
end
xline(ax, [-90 90], "k--", HandleVisibility="off")
set(ax, FontSize=13, XTick=-180:45:180); xlim(ax, [-180 180]); ylim(ax, [0 1])
xlabel(ax, "Beam azimuth - flight heading (°)"); ylabel(ax, "Share of gates with clutter filter power ≥ 3 dB")
legend(ax, Location="south", FontSize=12)
title(ax, "Gates touched by the clutter filter", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "sweep_offset_azimuth.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function s = volume_stats(cache_file, range_edges_m, sector_edges, dbz_edges, difference_edges)
% Per sweep (pass) of one volume, by 1 km range bin and sector of the angle
% between beam and flight direction: gates valid (not range folded), sum
% of linear Z (below threshold = 0), gates with CFP >= 3 dB and with CFP
% recorded; dBZ histogram of all valid gates by range bin (no echo in the
% lowest bin). For split sweeps, the two passes compared on common gates.
volume = load(cache_file).volume;
% Flight heading: mean of the velocity profile where it is defined
uv = volume.uv(all(isfinite(volume.uv), 2), :);
heading = atan2d(mean(uv(:, 1)), mean(uv(:, 2)));
sweeps = bmm4d.radar.read_level2_moments(volume.source_file, ["REF" "VEL" "CFP"], 7);
elev = [sweeps.elev];
n_range = numel(range_edges_m) - 1;
n_sector = numel(sector_edges) - 1;
result = struct("elev", {}, "type", {}, "n_valid", {}, "sum_z", {}, "n_cfp", {}, "n_cfp_known", {}, "dbz_histogram", {});
split_elev = [];
split_z = zeros(n_range, 2, 0);
split_difference = zeros(numel(difference_edges) - 1, 0);
for i_s = 1:numel(sweeps)
    sweep = sweeps(i_s);
    if isempty(sweep.REF), continue, end
    twin = find(abs(elev - sweep.elev) < 0.1 & (1:numel(sweeps)) ~= i_s, 1);
    if isempty(sweep.VEL)
        type = "CS";
    elseif ~isempty(twin)
        type = "CD_split";
    elseif numel(sweep.REF_range) > numel(sweep.VEL_range)
        type = "batch";
    else
        type = "CD";
    end
    range_m = sweep.REF_range;
    code = sweep.REF_code;
    z = 10 .^ (sweep.REF / 10);
    z(code == 0) = 0;
    valid = code ~= 1;
    range_bin = discretize(range_m, range_edges_m);
    angle = mod(sweep.az - heading + 180, 360) - 180;
    sector = discretize(angle, sector_edges);
    [range_id, sector_id] = ndgrid(range_bin, sector);
    use = valid & ~isnan(range_id) & ~isnan(sector_id);
    subs = [range_id(use) sector_id(use)];
    n_valid = accumarray(subs, 1, [n_range n_sector]);
    sum_z = accumarray(subs, z(use), [n_range n_sector]);
    n_cfp = zeros(n_range, n_sector);
    n_cfp_known = zeros(n_range, n_sector);
    if ~isempty(sweep.CFP)
        cfp = sweep.CFP(1:numel(range_m), :);
        known = use & sweep.CFP_code(1:numel(range_m), :) ~= 0;
        n_cfp_known = accumarray([range_id(known) sector_id(known)], 1, [n_range n_sector]);
        filtered = known & cfp >= 3;
        n_cfp = accumarray([range_id(filtered) sector_id(filtered)], 1, [n_range n_sector]);
    end
    dbz = sweep.REF;
    dbz(code == 0) = dbz_edges(1);
    dbz_bin = discretize(min(max(dbz, dbz_edges(1)), dbz_edges(end) - 0.01), dbz_edges);
    hist_use = use & ~isnan(dbz_bin);
    dbz_histogram = accumarray([dbz_bin(hist_use) range_id(hist_use)], 1, [numel(dbz_edges) - 1 n_range]);
    result(end + 1) = struct("elev", sweep.elev, "type", type, "n_valid", n_valid, "sum_z", sum_z, ... %#ok<AGROW>
        "n_cfp", n_cfp, "n_cfp_known", n_cfp_known, "dbz_histogram", dbz_histogram);
    % Split sweep: this is the CS pass and the twin the CD pass
    if type == "CS" && ~isempty(twin) && numel(sweeps(twin).az) == numel(sweep.az) && ~isempty(sweeps(twin).REF)
        doppler = sweeps(twin);
        n_gate = min(numel(range_m), numel(doppler.REF_range));
        code_cs = code(1:n_gate, :);
        code_cd = doppler.REF_code(1:n_gate, :);
        both_valid = code_cs ~= 1 & code_cd ~= 1 & ~isnan(range_bin(1:n_gate));
        z_cs = z(1:n_gate, :);
        z_cd = 10 .^ (doppler.REF(1:n_gate, :) / 10);
        z_cd(code_cd == 0) = 0;
        bins = repmat(range_bin(1:n_gate), 1, numel(sweep.az));
        split_elev(end + 1) = sweep.elev; %#ok<AGROW>
        split_z(:, :, end + 1) = [accumarray(bins(both_valid), z_cs(both_valid), [n_range 1]), ...
            accumarray(bins(both_valid), z_cd(both_valid), [n_range 1])];
        both_echo = both_valid & code_cs > 1 & code_cd > 1;
        difference = doppler.REF(1:n_gate, :) - sweep.REF(1:n_gate, :);
        split_difference(:, end + 1) = histcounts(difference(both_echo), difference_edges)'; %#ok<AGROW>
    end
end
s = struct("sweeps", result, "split_elev", split_elev, "split_z", split_z, "split_difference", split_difference);
end

function height_m = slant2ground_height(range_km, elev)
[~, height_m] = slant2ground(range_km * 1000, elev);
end
