function sandbox_case_figure(S, results, file, i_t, figure_title)
%SANDBOX_CASE_FIGURE Prior mean, reconstruction and scores of sandbox runs.
%   RESULTS is a cell array of sandbox_run outputs. One row per run at
%   frame I_T: column density of the prior mean and of the reconstruction,
%   and the west-east section through KDMX of both. A last row has the
%   scores for all runs: the radar-distance curve (column / 150 km regional
%   mean, 1 km bins, mean over frames) and the share of the column below
%   1 km by distance (solution solid, prior mean dashed). The prior mean is
%   the one the run used (radar prior, or its layer mean for
%   prior_mean "profile"; zero for "zero"). Writes FILE (png).

g = S.g;
inv = S.c.inversion;
n_run = numel(results);
[ix, iy, layer] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
cell_air = bmm4d.grid.air_fraction(g, S.active_id);   % columns count the air part of the lowest cell
is_active = false(g.sz(1:2));
is_active(column_id) = true;
n_time = size(results{1}.field, 2);

% Prior mean of each run
prior_means = cell(n_run, 1);
for i_r = 1:n_run
    option = results{i_r}.option;
    radar_prior = double(S.prior_cache(sprintf("%g|%s|%g", option.lambda_profile, option.prior_interpolation, ...
        option.prior_length_scale_m)));
    switch string(option.prior_mean)
        case "profile"
            layer_mean = zeros(g.sz(3), n_time);
            for k = 1:n_time
                layer_mean(:, k) = accumarray(layer, radar_prior(:, k), [g.sz(3) 1], @mean, 0);
            end
            prior_means{i_r} = layer_mean(layer, :);
        case "zero"
            prior_means{i_r} = zeros(size(radar_prior));
        otherwise
            prior_means{i_r} = radar_prior;
    end
end

% Distance to the nearest radar
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
distance_km = inf(g.sz(1:2));
for k = 1:height(S.radar_selected)
    distance_km = min(distance_km, 6371 * acos(min(1, sind(lat_grid) .* sind(S.radar_selected.lat(k)) + ...
        cosd(lat_grid) .* cosd(S.radar_selected.lat(k)) .* cosd(lon_grid - S.radar_selected.lon(k)))));
end
fine_edges = 0:1:250;
fine_center = fine_edges(1:end - 1)' + 0.5;
fine_bin = discretize(distance_km, fine_edges);
sigma_cells = 150000 / mean([g.dx g.dy]);
active_weight = imgaussfilt(double(is_active), sigma_cells);
z_layer = squeeze(g.z);
low = z_layer(layer) < 1000;

states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
k_radar = find(S.radar_selected.name == "KDMX", 1);
[~, section_y] = min(abs(g.lat - S.radar_selected.lat(k_radar)));
fig = bmm4d.plot.create_export_figure(3600, 950 * (n_run + 1), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, n_run + 1, 4, TileSpacing="tight", Padding="tight");
title(t, sprintf("%s (maps and sections at %s UTC)", figure_title, string(S.c.times(i_t), "HH:mm")), FontSize=15)
colors = lines(n_run);
for i_r = 1:n_run
    fields = {prior_means{i_r}(:, i_t), double(results{i_r}.field(:, i_t))};
    map_titles = {"Prior mean", {"Reconstruction", sprintf("bias %.3f, trusted misfit %.3f, detail %.2f", ...
        results{i_r}.bias_200km, results{i_r}.misfit_trusted, results{i_r}.detail)}};
    section_titles = ["Prior mean, KDMX section" "Reconstruction, KDMX section"];
    for i_f = 1:2
        column = reshape(accumarray(column_id, fields{i_f} .* cell_air, [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
        column(~is_active) = nan;
        ax = nexttile(t);
        bmm4d.plot.plot_night_map(ax, g, background, column, [0 200], states, S.radar_selected, S.c.figure_bounds, ones(g.sz(1:2)));
        yline(ax, g.lat(section_y), "w--")
        title(ax, map_titles{i_f}, FontSize=12, FontWeight="normal")
        if i_f == 1
            ylabel(ax, split(results{i_r}.name, ", "), FontSize=12, FontWeight="bold", Color=colors(i_r, :), Interpreter="none")
        end
        if i_r == n_run && i_f == 1
            cb = colorbar(ax, "southoutside");
            cb.Label.String = "Column density (birds km^{-2})";
        end
    end
    for i_f = 1:2
        field = nan(g.sz);
        field(S.active_id) = fields{i_f};
        section = squeeze(field(:, section_y, :))';
        ax = nexttile(t);
        imagesc(ax, g.lon, z_layer, section, AlphaData=~isnan(section))
        hold(ax, "on")
        plot(ax, S.radar_selected.lon(k_radar), 300, "kv", MarkerFaceColor="k")
        set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12)
        clim(ax, [0 150])
        colormap(ax, parula)
        xlim(ax, S.radar_selected.lon(k_radar) + [-6 6])
        ylim(ax, [0 3000])
        ylabel(ax, "Height (m a.s.l.)")
        title(ax, section_titles(i_f), FontSize=12, FontWeight="normal")
        if i_r == n_run && i_f == 1
            cb = colorbar(ax, "southoutside");
            cb.Label.String = "Density (birds km^{-3})";
        end
    end
end

% Scores of all runs
ax_curve = nexttile(t, [1 2]);
hold(ax_curve, "on"); box(ax_curve, "on"); grid(ax_curve, "on")
ax_low = nexttile(t, [1 2]);
hold(ax_low, "on"); box(ax_low, "on"); grid(ax_low, "on")
use = is_active & ~isnan(fine_bin);
for i_r = 1:n_run
    field = mean(double(results{i_r}.field), 2);
    column = reshape(accumarray(column_id, field .* cell_air, [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    column = max(column, 0) .* is_active;
    regional = imgaussfilt(column, sigma_cells) ./ active_weight;
    curve = accumarray(fine_bin(use), column(use) ./ max(regional(use), 1), [numel(fine_center) 1], @mean, nan);
    plot(ax_curve, fine_center, curve / mean(curve, "omitnan"), Color=colors(i_r, :), LineWidth=1.5, ...
        DisplayName=results{i_r}.name)
    sources = {field, mean(prior_means{i_r}, 2)};
    styles = ["-" "--"];
    for i_s = 1:2
        low_column = accumarray(column_id, sources{i_s} .* cell_air .* low, [prod(g.sz(1:2)) 1]);
        all_column = accumarray(column_id, sources{i_s} .* cell_air, [prod(g.sz(1:2)) 1]);
        share = accumarray(fine_bin(use), low_column(use) ./ max(all_column(use), eps), [numel(fine_center) 1], @mean, nan);
        plot(ax_low, fine_center, movmean(share, 5, "omitnan"), styles(i_s), Color=colors(i_r, :), LineWidth=1.5, ...
            HandleVisibility=string(matlab.lang.OnOffSwitchState(i_s == 1)), DisplayName=results{i_r}.name)
    end
end
yline(ax_curve, 1, "k--", HandleVisibility="off")
xlim(ax_curve, [0 250]); ylim(ax_curve, [0.6 1.4]); set(ax_curve, FontSize=12)
xlabel(ax_curve, "Distance to the nearest radar (km)"); ylabel(ax_curve, "Column / 150 km regional mean")
legend(ax_curve, Location="northeast", FontSize=11, Interpreter="none")
title(ax_curve, "Radar-distance curve, 1 km bins (flat at 1 = no imprint)", FontWeight="normal")
xlim(ax_low, [0 250]); set(ax_low, FontSize=12)
xlabel(ax_low, "Distance to the nearest radar (km)"); ylabel(ax_low, "Share of the column below 1 km")
legend(ax_low, Location="northeast", FontSize=11, Interpreter="none")
title(ax_low, "Share below 1 km: reconstruction (solid), prior mean (dashed)", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, file, Resolution=150)
close(fig)
end
