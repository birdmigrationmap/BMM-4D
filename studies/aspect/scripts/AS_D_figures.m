%% Aspect D: figures (dev)
% From the results of AS_B_overlap.m:
% - AS_aspect_function.png: (a) the aspect function measured on the radar
%   overlap without a shape (24 bins of psi = azimuth - track) and the
%   shared model at four heights; (b) side-on (b2) and front/back (b1)
%   amplitudes against height: overlap, single-radar population means, and
%   the track turning in time.
% - AS_kdmx_maps.png: KDMX 06:30 UTC 22 Oct 2026, mean density of the
%   0.5-1.3 deg sweeps within 90 km, uncorrected and with the three
%   corrections (top), and the factor each correction applied (bottom;
%   corrected / uncorrected), with the azimuthal profile at 10-60 km.
% - AS_overlap_validation.png: disagreement of overlapping radars and the
%   aspect harmonics left after each correction (out of fold), the aspect
%   function before and after the recommended correction, and the
%   single-radar harmonics at 10-60 km.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"))
data_root = bmm4d.util.local_data_dir();
figure_dir = fullfile(project_root, "studies", "aspect", "figures");
results = load(fullfile(data_root, "aspect", "overlap_results.mat"));
model = results.models.height_ramp;
side_color = [0 0.447 0.741];
front_color = [0.850 0.325 0.098];
font = {"FontName", "Helvetica", "FontSize", 13};

%% Figure 1: the aspect function
fig = bmm4d.plot.create_export_figure(3600, 1800, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
tiledlayout(fig, 1, 2, "TileSpacing", "tight", "Padding", "tight");

% (a) T(psi): overlap without a shape, and the shared model at four heights
ax = nexttile;
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
curve = results.curve_raw;
psi_line = (0:2:360)';
heights_km = [0.5 1 2 3];
model_colors = parula(numel(heights_km) + 1);
straight_profile = repmat([0 10], 30, 1);   % track 0 deg: azimuth = psi
for i = 1:numel(heights_km)
    template = as_aspect_template(psi_line, heights_km(i) * 1000 * ones(size(psi_line)), straight_profile, model);
    plot(ax, psi_line, template, "-", "Color", model_colors(i, :), "LineWidth", 2, ...
        "DisplayName", sprintf("model, %.1f km", heights_km(i)))
end
errorbar(ax, curve.psi, exp(curve.log_t), exp(curve.log_t) - exp(curve.log_t - 2 * curve.se), ...
    exp(curve.log_t + 2 * curve.se) - exp(curve.log_t), "o", "Color", [0.25 0.25 0.25], ...
    "MarkerFaceColor", [0.25 0.25 0.25], "MarkerSize", 6, "LineWidth", 1, "CapSize", 0, ...
    "DisplayName", "overlap (\pm2 se)")
set(ax, "YScale", "log", "XLim", [0 360], "XTick", 0:45:360, "YLim", [0.6 1.6], "YTick", [0.6 0.7 0.8 0.9 1 1.2 1.4 1.6])
xlabel(ax, "\psi = azimuth from radar - track (deg)")
ylabel(ax, "Relative reflectivity T(\psi)")
text(ax, [0 90 180 270] + 4, [0.62 0.62 0.62 0.62], ["tail-on" "side-on" "head-on" "side-on"], font{:}, "Color", [0.3 0.3 0.3])
legend(ax, "Location", "southoutside", "Box", "off", font{:}, "NumColumns", 3)
title(ax, "a  Aspect function", font{:}, "FontWeight", "bold")
set(ax, font{:})

% (b) Amplitudes against height
ax = nexttile;
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
overlap = results.overlap_by_height;
single = results.single_by_height;
h_line = (0:50:3500)';
plot(ax, max(0, model.b2_per_km * h_line / 1000), h_line / 1000, "-", "Color", side_color, "LineWidth", 2, ...
    "DisplayName", "model b_2 = 0.092 h")
plot(ax, model.b1 * min(1, h_line / model.b1_ramp_m), h_line / 1000, "-", "Color", front_color, "LineWidth", 2, ...
    "DisplayName", "model b_1 (ramp < 1 km)")
errorbar(ax, overlap.c2, overlap.height_m / 1000, 2 * overlap.c2_se, "horizontal", "o", "Color", side_color, ...
    "MarkerFaceColor", side_color, "MarkerSize", 8, "CapSize", 0, "LineWidth", 1.2, "DisplayName", "b_2 overlap")
errorbar(ax, overlap.c1, overlap.height_m / 1000, 2 * overlap.c1_se, "horizontal", "s", "Color", front_color, ...
    "MarkerFaceColor", front_color, "MarkerSize", 8, "CapSize", 0, "LineWidth", 1.2, "DisplayName", "b_1 overlap")
plot(ax, single.c2, single.height_m / 1000, "o--", "Color", side_color, "MarkerSize", 7, "LineWidth", 1, ...
    "DisplayName", "b_2 single radar")
plot(ax, single.c1, single.height_m / 1000, "s--", "Color", front_color, "MarkerSize", 7, "LineWidth", 1, ...
    "DisplayName", "b_1 single radar")
plot(ax, single.time_b2, single.height_m / 1000, "x:", "Color", side_color, "MarkerSize", 9, "LineWidth", 1, ...
    "DisplayName", "b_2 track turning")
plot(ax, single.time_b1, single.height_m / 1000, "x:", "Color", front_color, "MarkerSize", 9, "LineWidth", 1, ...
    "DisplayName", "b_1 track turning")
xline(ax, 0, "k-", "HandleVisibility", "off")
set(ax, "XLim", [-0.45 0.35], "YLim", [0 3.5])
xlabel(ax, "Amplitude (log units): b_2 side-on, b_1 tail-on minus head-on")
ylabel(ax, "Height above antenna (km)")
legend(ax, "Location", "southoutside", "Box", "off", "FontName", "Helvetica", "FontSize", 12, "NumColumns", 2)
title(ax, "b  Amplitudes against height", font{:}, "FontWeight", "bold")
set(ax, font{:})
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AS_aspect_function.png"), Resolution=300)
close(fig)

%% Figure 2: KDMX 06:30 UTC
volume_root = fullfile(data_root, "prod", "volumes");
file_name = "KDMX20261022_063022_V06.mat";
raw = as_read_volume(fullfile(volume_root, "us_east_fix", file_name), Inf);
versions = {raw, load(fullfile(volume_root, "9c2b4ef2", file_name)).volume, ...
    load(fullfile(volume_root, "us_east_fix", file_name)).volume, as_apply_aspect(raw, model)};
names = ["no correction" "ellipse (paper)" "per-volume symmetric" "shared aspect function"];
sweep_ids = 1:3;
max_range_m = 90000;
maps = cell(1, 4);
for i = 1:4
    [maps{i}, x_km, y_km] = as_volume_map(versions{i}, sweep_ids, max_range_m, 3000);
end
c_lim = [0 prctile(maps{1}(:), 99)];
track_uv = raw.uv(6, :);   % 500-600 m above the antenna
fig = bmm4d.plot.create_export_figure(3600, 1900, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 2, 4, "TileSpacing", "tight", "Padding", "tight");
ratio_map = [linspace(0.13, 1, 128)' linspace(0.40, 1, 128)' linspace(0.70, 1, 128)'; ...
    linspace(1, 0.70, 128)' linspace(1, 0.25, 128)' linspace(1, 0.15, 128)'];
for i = 1:4
    ax = nexttile(layout, i);
    % Transpose: maps are [x y], images are [row = y, column = x]
    imagesc(ax, x_km, y_km, maps{i}', "AlphaData", isfinite(maps{i}'))
    axis(ax, "xy", "equal", "tight")
    colormap(ax, parula(256))
    clim(ax, c_lim)
    set(ax, "Color", [0.96 0.96 0.96])
    hold(ax, "on")
    quiver(ax, 0, 0, 45 * track_uv(1) / norm(track_uv), 45 * track_uv(2) / norm(track_uv), 0, "k", "LineWidth", 2, ...
        "MaxHeadSize", 0.6)
    title(ax, names(i), font{:}, "FontWeight", "bold")
    xlabel(ax, "x (km)")
    if i == 1
        ylabel(ax, "y (km)")
    end
    set(ax, font{:})
    if i == 4
        cb = colorbar(ax);
        cb.Label.String = "Bird density (birds km^{-3})";
        cb.FontSize = 13;
        cb.Color = "k";
    end
end
% Azimuthal profile at 10-60 km (sweeps 0.5-1.3 deg)
ax = nexttile(layout, 5);
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
profile_colors = [0.2 0.2 0.2; 0.85 0.33 0.10; 0.47 0.67 0.19; 0 0.45 0.74];
centre = 5:10:355;
for i = 1:4
    [total, count] = deal(zeros(36, 1));
    for i_elev = sweep_ids
        sweep = versions{i}.dens.sweeps(i_elev);
        rows = sweep.range >= 10000 & sweep.range <= 60000;
        [i_range, i_az] = find(isfinite(sweep.data(rows, :)));
        values = sweep.data(rows, :);
        values = values(sub2ind(size(values), i_range, i_az));
        sector = min(36, floor(sweep.az(i_az) / 10) + 1);
        total = total + accumarray(sector(:), values, [36 1]);
        count = count + accumarray(sector(:), 1, [36 1]);
    end
    plot(ax, centre, total ./ count, "-", "Color", profile_colors(i, :), "LineWidth", 2, "DisplayName", names(i))
end
track_deg = mod(atan2d(track_uv(1), track_uv(2)), 360);
xline(ax, track_deg, "k--", "HandleVisibility", "off")
xline(ax, mod(track_deg + 180, 360), "k:", "HandleVisibility", "off")
set(ax, "XLim", [0 360], "XTick", 0:90:360, "YLim", [0 1.8 * max(ylim(ax))])
xlabel(ax, "Azimuth (deg; -- track, : upstream)")
ylabel(ax, "Mean density, 10-60 km")
legend(ax, "Location", "northwest", "Box", "off", "FontSize", 10)
title(ax, "Azimuthal profile", font{:}, "FontWeight", "bold")
set(ax, font{:})
for i = 2:4
    ax = nexttile(layout, i + 4);
    factor = maps{i} ./ maps{1};
    imagesc(ax, x_km, y_km, log2(factor)', "AlphaData", isfinite(factor'))
    axis(ax, "xy", "equal", "tight")
    colormap(ax, ratio_map)
    clim(ax, [-1 1])
    set(ax, "Color", [0.96 0.96 0.96])
    title(ax, names(i) + " / none", font{:}, "FontWeight", "bold")
    xlabel(ax, "x (km)")
    set(ax, font{:})
    if i == 4
        cb = colorbar(ax, "Ticks", log2([0.5 0.7 1 1.4 2]), "TickLabels", ["0.5" "0.7" "1" "1.4" "2"]);
        cb.Label.String = "Correction factor";
        cb.FontSize = 13;
        cb.Color = "k";
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AS_kdmx_maps.png"), Resolution=300)
close(fig)

%% Figure 3: validation on the overlap
scores = results.overlap_scores;
single_scores = results.single_scores;
labels = ["none" "ellipse (paper)" "per-volume symmetric" "shared symmetric" "shared b_2 + b_1" ...
    "shared b_2(h) + b_1" "shared b_2(h, v) + b_1" "shared b_2(h) + b_1 ramp"];
n_method = height(scores);
rows = 1:n_method;
fig = bmm4d.plot.create_export_figure(3600, 2600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "tight");
gray = [0.55 0.55 0.55];
method_axis = @(ax) set(ax, "YDir", "reverse", "YTick", rows, "YTickLabel", labels, "YLim", [0.5 n_method + 0.5]);

ax = nexttile(layout);
reduction = 100 * (1 - scores.disagreement / scores.disagreement(1));
barh(ax, rows, reduction, 0.6, "FaceColor", gray, "EdgeColor", "none")
method_axis(ax)
xlabel(ax, "Reduction of the overlap disagreement (%, RMS log ratio)")
title(ax, "a  Disagreement of overlapping radars (out of fold)", font{:}, "FontWeight", "bold")
grid(ax, "on")
set(ax, font{:})

ax = nexttile(layout);
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
errorbar(ax, scores.c2_left, rows' - 0.12, 2 * scores.c2_se, "horizontal", "o", "Color", side_color, ...
    "MarkerFaceColor", side_color, "MarkerSize", 8, "CapSize", 0, "LineWidth", 1.2, "DisplayName", "side-on b_2 left")
errorbar(ax, scores.c1_left, rows' + 0.12, 2 * scores.c1_se, "horizontal", "s", "Color", front_color, ...
    "MarkerFaceColor", front_color, "MarkerSize", 8, "CapSize", 0, "LineWidth", 1.2, "DisplayName", "front/back b_1 left")
xline(ax, 0, "k-", "HandleVisibility", "off")
method_axis(ax)
xlabel(ax, "Aspect harmonic left on the overlap (\pm2 se)")
legend(ax, "Location", "southwest", "Box", "off", font{:})
title(ax, "b  Aspect left on the overlap (out of fold)", font{:}, "FontWeight", "bold")
set(ax, font{:})

ax = nexttile(layout);
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
before = results.curve_raw;
after = results.curve_corrected;
errorbar(ax, before.psi, exp(before.log_t), exp(before.log_t) - exp(before.log_t - 2 * before.se), ...
    exp(before.log_t + 2 * before.se) - exp(before.log_t), "o-", "Color", gray, "MarkerFaceColor", gray, ...
    "CapSize", 0, "LineWidth", 1.2, "DisplayName", "none")
errorbar(ax, after.psi, exp(after.log_t), exp(after.log_t) - exp(after.log_t - 2 * after.se), ...
    exp(after.log_t + 2 * after.se) - exp(after.log_t), "o-", "Color", side_color, "MarkerFaceColor", side_color, ...
    "CapSize", 0, "LineWidth", 1.2, "DisplayName", "shared b_2(h) + b_1 ramp (out of fold)")
set(ax, "YScale", "log", "XLim", [0 360], "XTick", 0:45:360, "YLim", [0.6 1.6], "YTick", [0.6 0.7 0.8 0.9 1 1.2 1.4 1.6])
text(ax, [0 90 180 270] + 4, [0.62 0.62 0.62 0.62], ["tail-on" "side-on" "head-on" "side-on"], font{:}, "Color", [0.3 0.3 0.3])
xlabel(ax, "\psi = azimuth from radar - track (deg)")
ylabel(ax, "Relative reflectivity on the overlap")
legend(ax, "Location", "northwest", "Box", "off", font{:})
title(ax, "c  Aspect function on the overlap, before and after", font{:}, "FontWeight", "bold")
set(ax, font{:})

ax = nexttile(layout);
hold(ax, "on")
box(ax, "on")
grid(ax, "on")
plot(ax, single_scores.first_median, rows, "d", "Color", "k", "MarkerFaceColor", gray, "MarkerSize", 10, ...
    "DisplayName", "first harmonic, median (gradient + front/back)")
plot(ax, single_scores.c2_mean, rows, "o", "Color", side_color, "MarkerFaceColor", side_color, ...
    "MarkerSize", 9, "DisplayName", "side-on b_2, mean")
plot(ax, single_scores.c1_mean, rows, "s", "Color", front_color, "MarkerFaceColor", front_color, ...
    "MarkerSize", 9, "DisplayName", "tail-on minus head-on b_1, mean")
xline(ax, 0, "k-", "HandleVisibility", "off")
method_axis(ax)
set(ax, "XLim", [-0.2 0.32])
xlabel(ax, "Single-radar azimuthal harmonics at 10-60 km")
legend(ax, "Location", "southoutside", "Box", "off", font{:})
title(ax, "d  Each radar, 10-60 km (gradients included)", font{:}, "FontWeight", "bold")
set(ax, font{:})
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AS_overlap_validation.png"), Resolution=300)
close(fig)
