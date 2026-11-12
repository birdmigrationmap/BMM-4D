%% Gate density K: does the clutter filter remove birds?
% Zero-isodop test on the GD_J reference scans, which does not use the range
% bias. Birds whose radial velocity is near zero fall in the filter's notch;
% ground clutter does not depend on the birds' flight direction. A VAD fit to
% the unfiltered gates of each 1 km ring gives every gate the radial velocity
% the birds should have there. If the filter removes birds, gates near zero
% expected radial velocity lose more: more gates filtered, a lower output,
% and an output plus CFP that is not lower. Without the filter, side-on birds
% at the zero isodop would give a higher output, not a lower one.
% Control: the same with the fitted flight direction turned by 90 degrees,
% which keeps the clutter and the scan but moves the zero isodop.
% Two groups of cuts at 4-40 km: the lowest cuts (-0.2-1.3 deg, split: the
% surveillance pass has REF and CFP, the Doppler pass VEL, matched by azimuth
% and range) and the batch cuts (1.8-5.1 deg, REF, VEL and CFP on the same
% radials). Gates that lost more than MAX_CFP_DB are taken as strong clutter
% and left out of the densities.
% Figures: one example scan (output and CFP maps at the lowest cut and at
% 2.4 deg with the zero isodop, and the share of filtered gates by azimuth
% relative to the flight direction), a summary over all scans, and the full
% ratios with the control.
% Outputs: studies/gate_density/figures/GD_K_cfp_*.png, GD_K_cfp_birds.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");
scans = load(fullfile(bmm4d.util.local_data_dir(), "gate_density", "post2020_scans", "scans.mat")).scans;
n_scan = height(scans);

%% Settings ---------------------------------------------------------------
max_elevation_deg = 5.5;
group_edges_deg = [-1 1.6 max_elevation_deg];        % lowest (split) cuts, some radars below 0 deg; batch cuts
group_name = ["≤ 1.3°" "1.8-5.1°"];
ring_edges_km = 4:1:40;
range_edges_km = [4 10 20 40];
n_range_band = numel(range_edges_km) - 1;
% One band per group of cuts and range band
band_group = repelem(1:2, n_range_band);
band_range = repmat(1:n_range_band, 1, 2);
band_label = group_name(band_group) + ", " + string(range_edges_km(band_range)) + "-" + string(range_edges_km(band_range + 1)) + " km";
vr_edges = [0 1 2 3 4 6 8 Inf];                      % |expected radial velocity| (m/s)
vr_center = [0.5 1.5 2.5 3.5 5 7 10];
n_vr = numel(vr_center);
n_band = numel(band_label);
min_speed = 5;                                       % VAD speed (m/s) for a clear zero isodop
max_rms = 5;                                         % VAD residual (m/s)
min_fit_gates = 100;
max_cfp_db = 20;

%% Read the scans and accumulate by expected radial velocity ----------------------------
% Per |v_r| bin, band, scan and geometry (fitted, turned 90 deg): all gates, gates
% that lost >= 3 dB, and output and output plus CFP summed over the gates that
% lost at most MAX_CFP_DB
geometry_name = ["fitted flight direction" "control: turned 90°"];
n_all = zeros(n_vr, n_band, n_scan, 2);
n_lost = zeros(n_vr, n_band, n_scan, 2);
n_kept = zeros(n_vr, n_band, n_scan, 2);
sum_output = zeros(n_vr, n_band, n_scan, 2);
sum_before = zeros(n_vr, n_band, n_scan, 2);
% CFP strength: per CFP class, band and scan, over azimuths where the birds' expected
% |v_r| >= 4 m/s (outside the notch): gates, their output and output plus CFP, and the
% mean output of the unfiltered gates of the same ring times the number of gates
cfp_class_edges = [-Inf -0.5 0.5 2.5 5.5 10.5 20.5 30.5 Inf];   % -1: filter not applied
cfp_class_name = ["not applied" "0" "1-2" "3-5" "6-10" "11-20" "21-30" ">30"];
n_class = numel(cfp_class_name);
class_count = zeros(n_class, n_band, n_scan);
class_output = zeros(n_class, n_band, n_scan);
class_before = zeros(n_class, n_band, n_scan);
class_reference = zeros(n_class, n_band, n_scan);
tic
for i_s = 1:n_scan
    level2 = bmm4d.radar.read_level2_moments(scans.raw_file(i_s), ["REF" "VEL" "CFP"], max_elevation_deg);
    has_cfp = arrayfun(@(s) ~isempty(s.REF) && ~isempty(s.CFP), level2);
    for i_e = reshape(find(has_cfp), 1, [])
        cut = level2(i_e);
        % Velocity from the same cut, or from the Doppler pass of a split cut
        vel = gd_cut_velocity(level2, i_e);
        if isempty(vel)
            continue
        end
        n_gate = min(size(cut.REF, 1), size(cut.CFP_code, 1));
        range_km = cut.REF_range(1:n_gate) / 1000;
        dbz = cut.REF(1:n_gate, :);
        code = double(cut.CFP_code(1:n_gate, :));
        vel = vel(1:n_gate, :);
        i_group = discretize(cut.elev, group_edges_deg);
        filtered = code >= 8;
        cfp_db = max(code - 8, 0);
        output = z_to_refl(idb(dbz)) / 11;
        output(isnan(dbz)) = 0;
        before = output .* 10 .^ (cfp_db / 10);
        az = cut.az(:)';

        for i_ring = 1:numel(ring_edges_km) - 1
            rows = range_km >= ring_edges_km(i_ring) & range_km < ring_edges_km(i_ring + 1);
            % VAD fit on the unfiltered gates with a velocity
            use = ~filtered(rows, :) & ~isnan(vel(rows, :));
            if nnz(use) < min_fit_gates
                continue
            end
            [coef, fit_rms] = gd_vad_fit(vel(rows, :), az, use);
            speed = hypot(coef(2), coef(3));
            if fit_rms > max_rms || speed < min_speed
                continue
            end
            % Expected |v_r| from the fit, and with (u, v) turned by 90 deg for the control
            expected_vr = [abs(coef(1) + coef(2) * sind(az) + coef(3) * cosd(az)); ...
                abs(coef(1) + coef(3) * sind(az) - coef(2) * cosd(az))];
            i_band = (i_group - 1) * n_range_band + discretize(mean(range_km(rows)), range_edges_km);
            kept = cfp_db(rows, :) <= max_cfp_db;
            lost = code(rows, :) >= 11;
            o = output(rows, :);
            b = before(rows, :);
            for i_g = 1:2
                i_vr = discretize(repmat(expected_vr(i_g, :), nnz(rows), 1), vr_edges);
                n_all(:, i_band, i_s, i_g) = n_all(:, i_band, i_s, i_g) + accumarray(i_vr(:), 1, [n_vr 1]);
                n_lost(:, i_band, i_s, i_g) = n_lost(:, i_band, i_s, i_g) + accumarray(i_vr(:), lost(:), [n_vr 1]);
                n_kept(:, i_band, i_s, i_g) = n_kept(:, i_band, i_s, i_g) + accumarray(i_vr(kept), 1, [n_vr 1]);
                sum_output(:, i_band, i_s, i_g) = sum_output(:, i_band, i_s, i_g) + accumarray(i_vr(kept), o(kept), [n_vr 1]);
                sum_before(:, i_band, i_s, i_g) = sum_before(:, i_band, i_s, i_g) + accumarray(i_vr(kept), b(kept), [n_vr 1]);
            end
            % CFP strength, outside the notch
            fast = expected_vr(1, :) >= 4;
            c = code(rows, fast);
            unfiltered = c < 8;
            if nnz(unfiltered) < 50
                continue
            end
            o_fast = o(:, fast);
            b_fast = b(:, fast);
            i_class = discretize(max(c - 8, -1), cfp_class_edges);
            n = accumarray(i_class(:), 1, [n_class 1]);
            class_count(:, i_band, i_s) = class_count(:, i_band, i_s) + n;
            class_output(:, i_band, i_s) = class_output(:, i_band, i_s) + accumarray(i_class(:), o_fast(:), [n_class 1]);
            class_before(:, i_band, i_s) = class_before(:, i_band, i_s) + accumarray(i_class(:), b_fast(:), [n_class 1]);
            class_reference(:, i_band, i_s) = class_reference(:, i_band, i_s) + n * mean(o_fast(unfiltered));
        end
    end
    fprintf("[cfp birds] %d/%d %s\n", i_s, n_scan, scans.station(i_s))
end
toc
save(fullfile(figure_dir, "GD_K_cfp_birds.mat"), "n_all", "n_lost", "n_kept", "sum_output", "sum_before", ...
    "class_count", "class_output", "class_before", "class_reference", "cfp_class_name", ...
    "geometry_name", "vr_edges", "vr_center", "band_label", "scans")

%% CFP strength: output of filtered gates relative to the unfiltered gates of the ring ---------
% Outside the notch (|v_r| >= 4 m/s). If the filter removes only clutter and the
% filtered gates hold as many birds, the output ratio stays near 1 at every CFP; if it
% removes birds with the clutter, it falls as CFP grows; if the filtered gates hold
% fewer birds, it is low already at small CFP. Ratios per scan, median over scans
% with at least 50 gates in the class.
output_ratio = class_output ./ class_reference;
before_ratio = class_before ./ class_reference;
output_ratio(class_count < 50) = nan;
before_ratio(class_count < 50) = nan;
strength = table(repelem(band_label', n_class), repmat(cfp_class_name', n_band, 1), VariableNames=["band" "cfp_db"]);
strength.share_of_gates = reshape(sum(class_count, 3) ./ sum(class_count, [1 3]), [], 1);
strength.n_scan = reshape(sum(class_count >= 50, 3), [], 1);
strength.output_ratio = reshape(median(output_ratio, 3, "omitnan"), [], 1);
strength.output_cfp_ratio = reshape(median(before_ratio, 3, "omitnan"), [], 1);
strength{:, 3:end} = round(strength{:, 3:end}, 3);
disp(strength)
writetable(strength, fullfile(figure_dir, "GD_K_cfp_strength.csv"))

fig = bmm4d.plot.create_export_figure(2800, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
strength_bands = [1 2 4 5];                          % <= 1.3 deg and 1.8-5.1 deg, 4-10 and 10-20 km
strength_colors = [0.15 0.25 0.6; 0.15 0.25 0.6; 0.85 0.5 0.1; 0.85 0.5 0.1];
strength_styles = ["-o" "--o" "-o" "--o"];
quantities = {output_ratio, before_ratio};
quantity_title = ["Output" "Output + CFP"];
for i_q = 1:2
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for k = 1:numel(strength_bands)
        plot(ax, 1:n_class, median(squeeze(quantities{i_q}(:, strength_bands(k), :)), 2, "omitnan"), strength_styles(k), ...
            Color=strength_colors(k, :), LineWidth=2, MarkerSize=5, DisplayName=band_label(strength_bands(k)))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, FontSize=13, YScale="log"); xlim(ax, [0.5 n_class + 0.5]); xticks(ax, 1:n_class); xticklabels(ax, cfp_class_name)
    xlabel(ax, "Power removed by the filter, CFP (dB)")
    title(ax, quantity_title(i_q), FontWeight="bold", FontSize=14)
    if i_q == 1
        ylabel(ax, "Relative to unfiltered gates of the ring"); ylim(ax, [0.01 3])
        legend(ax, Location="southwest", FontSize=12)
    else
        ylim(ax, [0.3 1000])
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_K_cfp_strength.png"), Resolution=150)
close(fig)

%% Zero isodop: per scan, each bin relative to the bins at |v_r| >= 4 m/s ----------------
% Ratios per scan (so that dense scans do not dominate), then the median over scans
reference = vr_center >= 4;
share_lost = n_lost ./ n_all;
ratio = {share_lost ./ (sum(n_lost(reference, :, :, :), 1) ./ sum(n_all(reference, :, :, :), 1)), ...
    (sum_output ./ n_kept) ./ (sum(sum_output(reference, :, :, :), 1) ./ sum(n_kept(reference, :, :, :), 1)), ...
    (sum_before ./ n_kept) ./ (sum(sum_before(reference, :, :, :), 1) ./ sum(n_kept(reference, :, :, :), 1))};
ratio_name = ["share lost ≥ 3 dB" "output" "output + CFP"];
ratio_column = ["share_lost_ratio" "output_ratio" "output_cfp_ratio"];
% Scans with at least 200 gates in every bin of the band, for both geometries
enough = all(n_all >= 200, [1 4]);
for i_q = 1:numel(ratio)
    ratio{i_q}(~repmat(enough, n_vr, 1, 1, 2)) = nan;
end
share_lost(~repmat(enough, n_vr, 1, 1, 2)) = nan;
summary = table();
for i_g = 1:2
    for i_b = 1:n_band
        rows = table(repmat(geometry_name(i_g), n_vr, 1), repmat(band_label(i_b), n_vr, 1), ...
            vr_center', repmat(nnz(enough(1, i_b, :)), n_vr, 1), median(share_lost(:, i_b, :, i_g), 3, "omitnan"), ...
            VariableNames=["geometry" "band" "vr_ms" "n_scan" "share_lost"]);
        for i_q = 1:numel(ratio)
            q = squeeze(ratio{i_q}(:, i_b, :, i_g));
            rows.(ratio_column(i_q) + "_median") = median(q, 2, "omitnan");
            rows.(ratio_column(i_q) + "_q25") = prctile(q, 25, 2);
            rows.(ratio_column(i_q) + "_q75") = prctile(q, 75, 2);
        end
        summary = [summary; rows]; %#ok<AGROW>
    end
end
summary{:, 5:end} = round(summary{:, 5:end}, 2);
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_K_cfp_birds.csv"))

%% Plot: all ratios, groups of cuts (rows) by range band (columns); control as dashed lines -------
fig = bmm4d.plot.create_export_figure(3600, 1900, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, n_range_band, TileSpacing="tight", Padding="tight");
quantity_colors = [0.5 0.5 0.5; 0.1 0.4 0.8; 0.42 0.27 0.66];
for i_b = 1:n_band
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_q = 1:numel(ratio)
        q = squeeze(ratio{i_q}(:, i_b, :, 1));
        fill(ax, [vr_center fliplr(vr_center)], [prctile(q, 25, 2)' fliplr(prctile(q, 75, 2)')], quantity_colors(i_q, :), ...
            FaceAlpha=0.15, EdgeColor="none", HandleVisibility="off")
        plot(ax, vr_center, median(q, 2, "omitnan"), "-o", Color=quantity_colors(i_q, :), LineWidth=2, DisplayName=ratio_name(i_q))
        plot(ax, vr_center, median(squeeze(ratio{i_q}(:, i_b, :, 2)), 2, "omitnan"), "--", Color=quantity_colors(i_q, :), ...
            LineWidth=1.2, DisplayName=ratio_name(i_q) + ", " + geometry_name(2))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, FontSize=13, YScale="log"); xlim(ax, [0 10]); ylim(ax, [0.4 3]); yticks(ax, [0.5 0.7 1 1.5 2 3])
    title(ax, sprintf("%s (%d scans)", band_label(i_b), nnz(enough(1, i_b, :))), FontWeight="bold", FontSize=14)
    if band_group(i_b) == 2
        xlabel(ax, "Expected |radial velocity| of the birds (m/s)")
    end
    if band_range(i_b) == 1
        ylabel(ax, "Relative to |v_r| ≥ 4 m/s")
    end
    if i_b == 1
        legend(ax, Location="northeast", FontSize=10)
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_K_cfp_birds.png"), Resolution=150)
close(fig)

%% Summary: share filtered and output by expected |v_r|, lowest and 2.4-5.1 deg cuts ----------
% Fitted flight direction, 4-10 and 10-20 km, median over scans
summary_bands = [1 2 4 5];                           % <= 1.3 deg and 1.8-5.1 deg, 4-10 and 10-20 km
band_colors = [0.15 0.25 0.6; 0.15 0.25 0.6; 0.85 0.5 0.1; 0.85 0.5 0.1];
band_styles = ["-o" "--o" "-o" "--o"];
fig = bmm4d.plot.create_export_figure(2800, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
ax1 = nexttile(t); hold(ax1, "on"); box(ax1, "on"); grid(ax1, "on")
ax2 = nexttile(t); hold(ax2, "on"); box(ax2, "on"); grid(ax2, "on")
for k = 1:numel(summary_bands)
    i_b = summary_bands(k);
    plot(ax1, vr_center, 100 * median(squeeze(share_lost(:, i_b, :, 1)), 2, "omitnan"), band_styles(k), Color=band_colors(k, :), ...
        LineWidth=2, MarkerSize=5, DisplayName=band_label(i_b))
    plot(ax2, vr_center, median(squeeze(ratio{2}(:, i_b, :, 1)), 2, "omitnan"), band_styles(k), Color=band_colors(k, :), ...
        LineWidth=2, MarkerSize=5, DisplayName=band_label(i_b))
end
yline(ax2, 1, "k:", HandleVisibility="off")
set([ax1 ax2], FontSize=13); xlim([ax1 ax2], [0 10])
ylim(ax1, [0 80]); ylim(ax2, [0.5 1.3])
xlabel([ax1 ax2], "Expected |radial velocity| of the birds (m/s)")
ylabel(ax1, "Gates that lost ≥ 3 dB (%)"); ylabel(ax2, "Output relative to |v_r| ≥ 4 m/s")
title(ax1, "Where the filter acts", FontWeight="bold", FontSize=14)
title(ax2, "What the output keeps", FontWeight="bold", FontSize=14)
legend(ax1, Location="northeast", FontSize=12)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_K_cfp_summary.png"), Resolution=150)
close(fig)

%% Example scan: output and CFP with the zero isodop ----------------------------------------
% The strongest night of KSGF; flight direction from the VAD rings at 4-20 km of each cut
[~, i_x] = max(scans.bird_density .* (scans.station == "KSGF"));
level2 = bmm4d.radar.read_level2_moments(scans.raw_file(i_x), ["REF" "VEL" "CFP"], max_elevation_deg);
has_cfp = find(arrayfun(@(s) ~isempty(s.REF) && ~isempty(s.CFP), level2));
example_elevation = [0.5 2.4];
map_range_km = 20;
angle_edges = -180:10:180;
angle_center = angle_edges(1:end - 1) + 5;
fig = bmm4d.plot.create_export_figure(3300, 2000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, 3, TileSpacing="tight", Padding="tight");
for i_c = 1:numel(example_elevation)
    [~, k] = min(abs([level2(has_cfp).elev] - example_elevation(i_c)));
    i_e = has_cfp(k);
    cut = level2(i_e);
    vel = gd_cut_velocity(level2, i_e);
    n_gate = find(cut.REF_range <= map_range_km * 1000, 1, "last");
    range_km = cut.REF_range(1:n_gate) / 1000;
    az = cut.az(:)';
    code = double(cut.CFP_code(1:n_gate, :));
    cfp_db = max(code - 8, 0);
    cfp_db(code < 8) = nan;                          % filter not applied
    output = z_to_refl(idb(cut.REF(1:n_gate, :))) / 11;
    output(isnan(cut.REF(1:n_gate, :))) = 0;
    % Flight direction: median of the VAD rings at 4-20 km
    ring_coef = nan(0, 3);
    for r = 4:19
        rows = range_km >= r & range_km < r + 1;
        use = code(rows, :) < 8 & ~isnan(vel(rows, :));
        if nnz(use) >= min_fit_gates
            [coef, fit_rms] = gd_vad_fit(vel(rows, :), az, use);
            if fit_rms <= max_rms
                ring_coef(end + 1, :) = coef'; %#ok<SAGROW>
            end
        end
    end
    heading = atan2d(median(ring_coef(:, 2)), median(ring_coef(:, 3)));
    speed = hypot(median(ring_coef(:, 2)), median(ring_coef(:, 3)));
    [az_grid, r_grid] = meshgrid([az az(1) + 360], range_km);
    x = r_grid .* sind(az_grid);
    y = r_grid .* cosd(az_grid);
    panels = {log10(max(output, 0.1)), cfp_db};
    for i_p = 1:2
        ax = nexttile(t); hold(ax, "on")
        surf(ax, x, y, zeros(size(x)), panels{i_p}(:, [1:end 1]), EdgeColor="none")
        view(ax, 2); axis(ax, "equal", "tight"); set(ax, FontSize=12, Color=[0.85 0.92 0.98])
        % Zero isodop (dashed) and flight direction (arrow)
        for side = [-90 90]
            plot3(ax, [0 map_range_km * sind(heading + side)], [0 map_range_km * cosd(heading + side)], [1 1], "r--", LineWidth=2)
        end
        quiver3(ax, 0, 0, 1, 12 * sind(heading), 12 * cosd(heading), 0, 0, "r", LineWidth=2.5, MaxHeadSize=0.6)
        if i_p == 1
            colormap(ax, parula); clim(ax, [0 3])
            cb = colorbar(ax, Color="k"); cb.Label.String = "log_{10} birds km^{-3}";
            title(ax, sprintf("%.1f°: radar output", cut.elev), FontWeight="bold", FontSize=14)
            ylabel(ax, "y (km)")
        else
            colormap(ax, flipud(gray(256))); clim(ax, [0 30])
            cb = colorbar(ax, Color="k"); cb.Label.String = "CFP (dB); blue: filter not applied";
            title(ax, sprintf("%.1f°: power removed by the filter", cut.elev), FontWeight="bold", FontSize=14)
        end
        if i_c == 2
            xlabel(ax, "x (km)")
        end
    end
    % Share of gates that lost >= 3 dB by azimuth relative to the flight direction, 4-15 km
    rows = range_km >= 4 & range_km < 15;
    relative = mod(az - heading + 180, 360) - 180;
    i_angle = discretize(repmat(relative, nnz(rows), 1), angle_edges);
    share = accumarray(i_angle(:), reshape(code(rows, :) >= 11, [], 1), [numel(angle_center) 1]) ./ ...
        accumarray(i_angle(:), 1, [numel(angle_center) 1]);
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    xline(ax, [-90 90], "r--", LineWidth=1.5)
    plot(ax, angle_center, 100 * share, "-o", Color=[0.3 0.3 0.3], LineWidth=2, MarkerSize=4)
    set(ax, FontSize=12); xlim(ax, [-180 180]); xticks(ax, -180:90:180); ylim(ax, [0 100])
    ylabel(ax, "Gates that lost ≥ 3 dB at 4-15 km (%)")
    title(ax, "Filtered gates by azimuth", FontWeight="bold", FontSize=14)
    if i_c == 2
        xlabel(ax, "Beam azimuth minus flight direction (°)")
    end
end
title(t, [sprintf("%s, %s UTC, %.0f birds km^{-3} at 20-60 km", scans.station(i_x), string(scans.time(i_x), "d MMMM yyyy HH:mm"), ...
    scans.bird_density(i_x)); sprintf("Birds at %.0f m/s towards %.0f° (red arrow); red dashed: zero isodop", speed, mod(heading, 360))], ...
    FontSize=15, FontWeight="bold", FontName="Helvetica", Color="k")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_K_cfp_example.png"), Resolution=150)
close(fig)
