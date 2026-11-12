%% Gate density H: example sweeps of the radar output and the estimator
% Maps within 40 km of the radar output, its moving mean over the
% estimator's window (GD_WINDOW_MEAN) and the estimator (no density cap,
% no clutter filter information) for scans drawn at random from the GD_J
% reference set: 4 per density class of the scan (mean output at 20-60 km
% on the lowest cut below 20, 20-100 and above 100 birds km^-3), each from a
% different radar. One figure per scan, with three elevation cuts (the
% lowest, about 1.5° and about 3.4°). One colour scale per scan, from 1/30
% to 3 times the estimator's median at 20-40 km on the lowest cut.
% Outputs: studies/gate_density/figures/GD_H_example_maps.csv and
% GD_H_example_<station>_<yyyyMMdd>.png

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
beam = bmm4d.radar.gate.beam();                                   % cross-section log-sd 1.25 (assumed, BMM4D.RADAR.GATE.BEAM)
table_s = bmm4d.radar.gate.bird_sum_table(beam, 10000, 20000);
s_floor_cm2 = 0.1;
radius_m = 1000;
map_range_m = 40000;
target_elevation = [0.5 1.5 3.4];
n_per_class = 4;
class_edges = [0 20 100 Inf];                        % birds km^-3, as GD_I
class_name = ["below 20" "20-100" "above 100"];

%% Pick the scans: random, one radar each ----------------------------------------
scans = load(fullfile(bmm4d.util.local_data_dir(), "gate_density", "post2020_scans", "scans.mat")).scans;
scans.density_class = discretize(scans.bird_density, class_edges);
rng(7)
picked = [];
for i_class = 1:numel(class_name)
    candidates = find(scans.density_class == i_class & ~ismember(scans.station, scans.station(picked)));
    candidates = candidates(randperm(numel(candidates)));
    [~, first] = unique(scans.station(candidates), "stable");
    picked = [picked; candidates(first(1:n_per_class))]; %#ok<AGROW>
end
examples = scans(picked, ["station" "time" "bird_density" "density_class"]);
examples.density_class = class_name(examples.density_class)';
examples.file = "GD_H_example_" + examples.station + "_" + string(examples.time, "yyyyMMdd") + ".png";

%% Estimate and plot each scan ------------------------------------------------------
for i_x = 1:height(examples)
    level2 = bmm4d.radar.read_level2_moments(scans.raw_file(picked(i_x)), "REF", max(target_elevation) + 0.5);
    level2 = level2(arrayfun(@(s) ~isempty(s.REF), level2));
    fig = bmm4d.plot.create_export_figure(2400, 2150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
    t = tiledlayout(fig, 3, numel(target_elevation), TileSpacing="tight", Padding="tight");
    for i_e = 1:numel(target_elevation)
        [~, i_cut] = min(abs([level2.elev] - target_elevation(i_e)));
        in_map = level2(i_cut).REF_range <= map_range_m;
        range_m = level2(i_cut).REF_range(in_map);
        dbz = level2(i_cut).REF(in_map, :);
        az_step_deg = 360 / size(dbz, 2);
        output = z_to_refl(idb(dbz)) / 11;
        output(isnan(dbz)) = 0;
        detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        [lo, hi] = bmm4d.radar.gate.bounds(dbz, isnan(dbz), range_m, beam, detection_dbz, s_floor_cm2);
        estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, az_step_deg, beam, table_s, radius_m=radius_m);
        if i_e == 1
            far = range_m >= 20000;
            limits = log10(median(estimate.density(far, :), "all", "omitnan")) + log10([1 / 30 3]);
        end
        az_edges = 0:az_step_deg:360;
        [az_grid, r_grid] = meshgrid(az_edges, range_m / 1000);
        x = r_grid .* sind(az_grid);
        y = r_grid .* cosd(az_grid);
        values = {output, gd_window_mean(output, range_m, az_step_deg, radius_m), estimate.density};
        row_label = ["Radar output" "Moving mean" "Estimator"];
        for i_p = 1:3
            ax = nexttile(t, (i_p - 1) * numel(target_elevation) + i_e); hold(ax, "on")
            plotted = log10(max(values{i_p}, 0.1));
            surf(ax, x, y, zeros(size(x)), plotted(:, [1:end 1]), EdgeColor="none")
            view(ax, 2); axis(ax, "equal", "tight"); set(ax, FontSize=12, Color="w")
            colormap(ax, parula); clim(ax, limits)
            if i_p == 1
                title(ax, sprintf("%.1f°", level2(i_cut).elev), FontWeight="bold", FontSize=14)
            elseif i_p == 3
                xlabel(ax, "x (km)")
            end
            if i_e == 1
                ylabel(ax, ["\bf" + row_label(i_p); "\rmy (km)"])
            end
        end
    end
    cb = colorbar(ax, Color="k"); cb.Layout.Tile = "east"; cb.Label.String = "log_{10} birds km^{-3}";
    title(t, sprintf("%s, %s UTC, %.0f birds km^{-3} at 20-60 km", examples.station(i_x), ...
        string(examples.time(i_x), "d MMMM yyyy HH:mm"), examples.bird_density(i_x)), FontSize=15, FontWeight="bold", FontName="Helvetica")
    bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, examples.file(i_x)), Resolution=150)
    close(fig)
    fprintf("[example maps] %d/%d %s\n", i_x, height(examples), examples.file(i_x))
end
disp(examples)
writetable(examples, fullfile(figure_dir, "GD_H_example_maps.csv"))
