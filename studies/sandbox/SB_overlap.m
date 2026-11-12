%% Sandbox: do neighbouring radars agree where they see the same air?
% On the volumes of SB_night_maps (9 dry nights, one volume per radar,
% same UTC time), each pair of radars less than 300 km apart is compared on
% the lowest cut, gate by gate, where both beams are at the same height
% (within 100 m) and about the same range (within 10 km): the same air,
% at the same height, with the same beam size. The difference in dBZ there
% is the pair's level difference; birds or beam geometry cannot explain it.
% Gates are dropped if either is missing (below threshold), is rain
% (REF >= 20 dBZ and RHO >= 0.97) or above 35 dBZ.
% Per night, one offset per radar (dB, mean zero) is fitted to the pair
% differences by least squares. An offset that repeats every night is a
% calibration (or processing) difference; one that changes is not.
% The same is done on the second cut (~0.9 deg, about 800 m higher at the
% matched gates): a calibration offset is the same on both cuts, an offset
% from the beam height (propagation, pointing) is not.
% Near a radar's detection threshold, reported values are censored and
% compressed, so offsets are also computed on strong gates only: both
% values at least STRONG_MARGIN_DB above their radar's threshold, taken
% per sweep as c + 20 log10(range) with c the median over 20-150 km of the
% 1st percentile of reported values minus 20 log10(range).
% Outputs: figures/overlap_pairs.csv, overlap_offsets.csv, overlap.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
radar = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
radar = radar(ismember(radar.name, c.radars), :);
n_radar = height(radar);
files = readtable(fullfile(figure_dir, "night_maps_files.csv"), TextType="string");
files.raw_file(ismissing(files.raw_file)) = "";
nights = unique(files.night);
n_night = numel(nights);
max_height_diff_m = 100;
max_range_diff_m = 10000;
strong_margin_db = 5;
earth_radius_m = 6371000;
effective_radius_m = 4 / 3 * earth_radius_m;

%% Read the lowest cut (REF, RHO) ----------------------------------------------
cut_limits = [0 0.7; 0.7 1.1];                 % elevation bands of the two cuts (deg)
cut_names = ["0.5°" "0.9°"];
n_cut = size(cut_limits, 1);
lowest = cell(height(files), n_cut);
parfor k = 1:height(files)
    if files.raw_file(k) == ""
        continue
    end
    try
        sweeps = bmm4d.radar.read_level2_moments(files.raw_file(k), ["REF" "RHO"], 1.1);
        sweeps = sweeps(arrayfun(@(s) ~isempty(s.REF) && ~isempty(s.RHO), sweeps));
        one_file = cell(1, 2);
        for i_c = 1:2
            in_band = find([sweeps.elev] >= cut_limits(i_c, 1) & [sweeps.elev] < cut_limits(i_c, 2), 1);
            if isempty(in_band)
                continue
            end
            s = sweeps(in_band);
            n_common = min(size(s.REF, 1), size(s.RHO, 1));
            ref = s.REF(1:n_common, :);
            % Detection threshold: c + 20 log10(range km), from the weakest reported values
            range_km = double(s.REF_range(1:n_common)) / 1000;
            in_fit = range_km >= 20 & range_km <= 150;
            threshold_c = median(prctile(ref(in_fit, :), 1, 2) - 20 * log10(range_km(in_fit)), "omitnan");
            rain = ref >= 20 & s.RHO(1:n_common, :) >= 0.97;
            ref(rain | ref > 35) = nan;
            one_file{i_c} = struct("REF", ref, "range", double(s.REF_range(1:n_common)), "az", double(s.az(:)), "elev", s.elev, ...
                "threshold_c", threshold_c);
        end
        lowest(k, :) = one_file;
    catch err
        fprintf("%s: %s\n", files.raw_file(k), err.message)
    end
end

%% Compare each pair on matched gates --------------------------------------------
radar_km = 6371 * acos(min(1, sind(radar.lat) .* sind(radar.lat') + cosd(radar.lat) .* cosd(radar.lat') .* cosd(radar.lon - radar.lon')));
[pair_a, pair_b] = find(triu(radar_km < 300 & radar_km > 0));
pairs = table();
for i_c = 1:n_cut
for i_n = 1:n_night
    for i_p = 1:numel(pair_a)
        a = pair_a(i_p);
        b = pair_b(i_p);
        sweep_a = lowest{files.night == nights(i_n) & files.station == radar.name(a), i_c};
        sweep_b = lowest{files.night == nights(i_n) & files.station == radar.name(b), i_c};
        if isempty(sweep_a) || isempty(sweep_b)
            continue
        end
        % Gates of A between 20 and 150 km, as ground positions
        keep_range = sweep_a.range >= 20000 & sweep_a.range <= 150000;
        [range_grid, az_grid] = ndgrid(sweep_a.range(keep_range), sweep_a.az);
        ref_a = sweep_a.REF(keep_range, :);
        ground_angle = range_grid / effective_radius_m * (effective_radius_m / earth_radius_m);
        lat = asind(sind(radar.lat(a)) .* cos(ground_angle) + cosd(radar.lat(a)) .* sin(ground_angle) .* cosd(az_grid));
        lon = radar.lon(a) + atan2d(sind(az_grid) .* sin(ground_angle) .* cosd(radar.lat(a)), ...
            cos(ground_angle) - sind(radar.lat(a)) .* sind(lat));
        % The same positions seen from B
        range_b = earth_radius_m * acos(min(1, sind(lat) .* sind(radar.lat(b)) + cosd(lat) .* cosd(radar.lat(b)) .* cosd(lon - radar.lon(b))));
        near = abs(range_b - range_grid) <= max_range_diff_m & range_b <= 150000 & ~isnan(ref_a);
        if nnz(near) < 100
            continue
        end
        az_b = mod(atan2d(sind(lon(near) - radar.lon(b)) .* cosd(lat(near)), ...
            cosd(radar.lat(b)) .* sind(lat(near)) - sind(radar.lat(b)) .* cosd(lat(near)) .* cosd(lon(near) - radar.lon(b))), 360);
        height_a = radar.height(a) + range_grid(near) * sind(sweep_a.elev) + range_grid(near) .^ 2 / (2 * effective_radius_m);
        height_b = radar.height(b) + range_b(near) * sind(sweep_b.elev) + range_b(near) .^ 2 / (2 * effective_radius_m);
        range_id = interp1(sweep_b.range, 1:numel(sweep_b.range), range_b(near), "nearest", NaN);
        [az_sorted, order] = sort(sweep_b.az);
        az_id = order(interp1([az_sorted - 360; az_sorted; az_sorted + 360], repmat((1:numel(az_sorted))', 3, 1), az_b, "nearest"));
        ref_b = nan(nnz(near), 1);
        ok = ~isnan(range_id);
        ref_b(ok) = sweep_b.REF(sub2ind(size(sweep_b.REF), range_id(ok), az_id(ok)));
        a_values = ref_a(near);
        match = abs(height_a - height_b) <= max_height_diff_m & ~isnan(ref_b);
        if nnz(match) < 100
            continue
        end
        difference = a_values(match) - ref_b(match);
        range_a_near = range_grid(near);
        range_b_near = range_b(near);
        % Strong gates: both well above their radar's detection threshold
        margin_a = a_values(match) - sweep_a.threshold_c - 20 * log10(range_a_near(match) / 1000);
        margin_b = ref_b(match) - sweep_b.threshold_c - 20 * log10(range_b_near(match) / 1000);
        strong = margin_a >= strong_margin_db & margin_b >= strong_margin_db;
        median_strong = nan;
        if nnz(strong) >= 50
            median_strong = median(difference(strong));
        end
        pairs = [pairs; table(i_c, nights(i_n), radar.name(a), radar.name(b), a, b, nnz(match), median(difference), ...
            10 * log10(mean(10 .^ (a_values(match) / 10)) / mean(10 .^ (ref_b(match) / 10))), ...
            median(height_a(match)), median(range_a_near(match)), median(range_b_near(match)), sweep_a.elev, sweep_b.elev, ...
            nnz(strong), median_strong, sweep_a.threshold_c, sweep_b.threshold_c, median(margin_a), median(margin_b), ...
            VariableNames=["cut" "night" "radar_a" "radar_b" "id_a" "id_b" "n_gate" ...
            "median_db" "mean_db" "height_m" "range_a_m" "range_b_m" "elev_a" "elev_b" ...
            "n_strong" "median_db_strong" "threshold_c_a" "threshold_c_b" "margin_a_db" "margin_b_db"])]; %#ok<AGROW>
    end
end
end
writetable(pairs, fullfile(figure_dir, "overlap_pairs.csv"))
for i_c = 1:n_cut
    fprintf("%s cut: %d pair-nights, median %d gates at %.0f m\n", cut_names(i_c), nnz(pairs.cut == i_c), ...
        median(pairs.n_gate(pairs.cut == i_c)), median(pairs.height_m(pairs.cut == i_c)))
    fprintf("    median margin above threshold %.1f dB; strong gates (>= %d dB): %d pair-nights with >= 50, median share %.2f\n", ...
        median([pairs.margin_a_db(pairs.cut == i_c); pairs.margin_b_db(pairs.cut == i_c)]), strong_margin_db, ...
        nnz(pairs.cut == i_c & ~isnan(pairs.median_db_strong)), median(pairs.n_strong(pairs.cut == i_c) ./ pairs.n_gate(pairs.cut == i_c)))
end

%% One offset per radar and night --------------------------------------------------
% median_db(a, b) = offset(a) - offset(b), offsets summing to zero
offset = nan(n_radar, n_night, n_cut);
offset_strong = nan(n_radar, n_night, n_cut);
for i_c = 1:n_cut
for i_n = 1:n_night
    one = pairs(pairs.cut == i_c & pairs.night == nights(i_n), :);
    offset_strong(:, i_n, i_c) = radar_offsets(one(~isnan(one.median_db_strong), :), "median_db_strong", n_radar);
    design = zeros(height(one) + 1, n_radar);
    design(sub2ind(size(design), (1:height(one))', one.id_a)) = 1;
    design(sub2ind(size(design), (1:height(one))', one.id_b)) = -1;
    design(end, :) = 1;
    seen = any(design(1:end - 1, :) ~= 0, 1);
    fit = design(:, seen) \ [one.median_db; 0];
    offset(seen, i_n, i_c) = fit;
    fprintf("%s %s: %d pairs, sd %.2f dB, residual rms after offsets %.2f dB\n", cut_names(i_c), ...
        string(nights(i_n), "yyyy-MM-dd"), height(one), std(one.median_db), rms(one.median_db - design(1:end - 1, seen) * fit))
end
end
% Lowest cut: mean over nights, night-to-night spread, share that repeats
low = offset(:, :, 1);
offset_mean = mean(low, 2, "omitnan");
offset_sd = std(low, 0, 2, "omitnan");
residual = low - offset_mean;
fprintf("0.5° cut: spread of the mean offsets %.2f dB; night-to-night sd %.2f dB; repeating share %.2f\n", ...
    std(offset_mean, "omitnan"), sqrt(mean(residual(:) .^ 2, "omitnan")), var(offset_mean, "omitnan") / var(low(:), "omitnan"))
% Same offset on both cuts? (radar-nights seen on both)
both = ~isnan(offset(:, :, 1)) & ~isnan(offset(:, :, 2));
low_values = offset(:, :, 1);
high_values = offset(:, :, 2);
cut_fit = polyfit(low_values(both), high_values(both), 1);
fprintf("0.9° vs 0.5° offsets (%d radar-nights): r = %.2f, slope %.2f\n", nnz(both), corr(low_values(both), high_values(both)), cut_fit(1))
both_strong = ~isnan(offset_strong(:, :, 1)) & ~isnan(offset_strong(:, :, 2));
low_strong = offset_strong(:, :, 1);
high_strong = offset_strong(:, :, 2);
strong_fit = polyfit(low_strong(both_strong), high_strong(both_strong), 1);
fprintf("strong gates only: 0.9° vs 0.5° offsets (%d radar-nights): r = %.2f, slope %.2f\n", nnz(both_strong), ...
    corr(low_strong(both_strong), high_strong(both_strong)), strong_fit(1))
residual_strong = low_strong - mean(low_strong, 2, "omitnan");
fprintf("strong gates, 0.5° cut: spread of the mean offsets %.2f dB; night-to-night sd %.2f dB\n", ...
    std(mean(low_strong, 2, "omitnan"), "omitnan"), sqrt(mean(residual_strong(:) .^ 2, "omitnan")))
% Do the offsets follow the detection thresholds (noise level) of the radars?
threshold = nan(n_radar, n_night);
for i_n = 1:n_night
    one = pairs(pairs.cut == 1 & pairs.night == nights(i_n), :);
    threshold(:, i_n) = accumarray([one.id_a; one.id_b], [one.threshold_c_a; one.threshold_c_b], [n_radar 1], @mean, nan);
end
threshold_relative = threshold - mean(threshold, 1, "omitnan");
ok = ~isnan(threshold_relative) & ~isnan(low);
fprintf("detection threshold (dBZ at 1 km): spread over radars %.2f dB; corr with the 0.5° offsets r = %.2f\n", ...
    std(threshold_relative(:), "omitnan"), corr(threshold_relative(ok), low(ok)))
offsets = table(radar.name, round(offset_mean, 2), round(offset_sd, 2), sum(~isnan(low), 2), ...
    VariableNames=["radar" "mean_offset_db" "sd_offset_db" "n_nights"]);
offsets = [offsets array2table(round(low, 2), VariableNames=cellstr(string(nights, "yyyy-MM-dd")))];
writetable(offsets, fullfile(figure_dir, "overlap_offsets.csv"))
disp(sortrows(offsets(:, 1:4), "mean_offset_db"))

%% Plot results -----------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 1500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 4, TileSpacing="tight", Padding="tight");
title(t, "Radar levels from matched gates of neighbouring radars (lowest cut, 9 dry nights)", FontSize=14)
% Offsets per radar: one dot per night, mean as a bar
[~, order] = sort(offset_mean);
ax = nexttile(t, [1 2]); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
bar(ax, 1:n_radar, offset_mean(order), FaceColor=[0.75 0.75 0.75], EdgeColor="none", DisplayName="mean")
night_colors = lines(n_night);
for i_n = 1:n_night
    scatter(ax, (1:n_radar) + (i_n - (n_night + 1) / 2) * 0.07, low(order, i_n), 36, night_colors(i_n, :), "filled", ...
        DisplayName=string(nights(i_n), "yyyy-MM-dd"))
end
xticks(ax, 1:n_radar); xticklabels(ax, radar.name(order)); xtickangle(ax, 60)
ylabel(ax, "Radar offset (dB, mean zero)"); set(ax, FontSize=12)
legend(ax, Location="northwest", NumColumns=3, FontSize=10)
title(ax, "0.5° cut, offset per radar: bar = mean over nights, dots = nights", FontWeight="normal")
% Map of the mean offsets
ax = nexttile(t); hold(ax, "on"); box(ax, "on")
set(ax, Color=[0.85 0.92 0.97], FontSize=12)
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
for s = states'
    plot(ax, s.Lon, s.Lat, Color=[0.5 0.5 0.5], LineWidth=0.5)
end
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
scatter(ax, radar.lon, radar.lat, 220, offset_mean, "filled", MarkerEdgeColor="k")
text(ax, radar.lon, radar.lat - 0.35, radar.name, FontSize=8, HorizontalAlignment="center")
colormap(ax, signed_map); clim(ax, [-3 3]); cb = colorbar(ax); cb.Label.String = "Mean offset (dB)";
xlim(ax, c.figure_bounds(1:2)); ylim(ax, c.figure_bounds(3:4))
title(ax, "Mean offset, 0.5° cut", FontWeight="normal")
% Same offset on both cuts?
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on"); axis(ax, "equal")
scatter(ax, low_values(both), high_values(both), 30, [0.3 0.3 0.3], "filled", MarkerFaceAlpha=0.6)
limits = [-12 8];
plot(ax, limits, limits, "k--")
xlim(ax, limits); ylim(ax, limits); set(ax, FontSize=12)
xlabel(ax, "Offset, 0.5° cut (dB)"); ylabel(ax, "Offset, 0.9° cut (dB)")
title(ax, sprintf("r = %.2f, slope %.2f", corr(low_values(both), high_values(both)), cut_fit(1)), FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "overlap.png"), Resolution=150)
close(fig)

function offset = radar_offsets(one, column, n_radar)
% One offset per radar (dB, mean zero over the radars seen) from the pair
% differences in ONE.(COLUMN) = offset(a) - offset(b), by least squares.
offset = nan(n_radar, 1);
if isempty(one)
    return
end
design = zeros(height(one) + 1, n_radar);
design(sub2ind(size(design), (1:height(one))', one.id_a)) = 1;
design(sub2ind(size(design), (1:height(one))', one.id_b)) = -1;
design(end, :) = 1;
seen = any(design(1:end - 1, :) ~= 0, 1);
% Minimum-norm solution: with few pairs the network can split in groups
offset(seen) = lsqminnorm(design(:, seen), [one.(column); 0]);
end
