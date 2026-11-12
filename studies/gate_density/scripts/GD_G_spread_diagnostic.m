%% Gate density G: where does the fitted cross-section spread come from?
% One cross-section log-sd fitted to all tiles is larger at 15-45 km than
% at 4-15 km. A property of the birds should not depend on range, so this
% script fits the log-sd separately by tile groups to see what it follows:
% - range and height of the beam centre (species, insects change with height);
% - expected number of birds per gate, lambda = density * V_sample: the
%   spread is identified from single-bird spikes where lambda is small;
%   where lambda is large, any extra gate-to-gate variation (patchiness,
%   residual clutter, background) can only be explained by a larger spread;
% - the radar's clutter filter (share of gates in the tile that lost
%   >= 3 dB, newer files only).
% It also compares, per tile, the estimate with the output mean: where
% gates hold many birds the output mean is reliable, so the log-sd that
% makes them agree is a second, independent check. The result motivates
% the fixed log-sd 1.25 of BMM4D.RADAR.GATE.BEAM (best fit at 10-30 birds per gate).
% Data: the 100 selected scans (as GD_E), cuts up to 10 deg, gates to 80 km,
% 2 km tiles. Outputs: studies/gate_density/figures/GD_G_*; the tiles are
% cached in the local data folder (delete GD_G_tiles.mat to recompute).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------------
scans_project_dir = "/Users/rafnuss/Library/CloudStorage/OneDrive-Vogelwarte/BMM-within";
sigma_grid = 0.5:0.25:2.5;
n_sigma = numel(sigma_grid);
density_grid = logspace(-1, 3.7, 60);
n_max = 10000;
s_floor_cm2 = 0.1;
tile_m = 2000;
min_known = 10;                                     % informative gates per tile
max_range_m = 80000;
max_elevation_deg = 10;
max_density = 3000;
earth_radius_m = 4 / 3 * 6371000;

%% Real sweeps of the 100 selected scans ------------------------------------------
scans = load(fullfile(scans_project_dir, "data", "derived", "scans.mat")).scans;
radar_table = readtable(fullfile(scans_project_dir, "data", "reference", "radar.csv"), TextType="string");
[~, base, ext] = fileparts(scans.filename);
scans.raw_file = fullfile(scans_project_dir, "data", "raw_scans", "nexrad", ...
    extractBefore(base, 5) + "_" + extractBetween(base, 5, 12), base + ext);
[~, i_radar] = ismember(scans.station, radar_table.name);
scans.antenna_agl = radar_table.height(i_radar) - radar_table.dem(i_radar);
n_scan = height(scans);

%% Tile log likelihood and estimate for every log-sd (cached) ------------------------
% Each worker reads its own scans and fits every log-sd (tables shared), so
% the sweeps never have to be held or sent in one piece. The tiles (about
% 3 million, 180 MB) are kept in the local data folder, not in the repository.
tiles_file = fullfile(bmm4d.util.local_data_dir(), "gate_density", "GD_G_tiles.mat");
if ~isfile(tiles_file)
    tables = cell(n_sigma, 1);
    beams = cell(n_sigma, 1);
    for i_g = 1:n_sigma
        beams{i_g} = bmm4d.radar.gate.beam(sigma_log_sd=sigma_grid(i_g));
        tables{i_g} = bmm4d.radar.gate.bird_sum_table(beams{i_g}, n_max, 20000);
    end
    if isempty(gcp("nocreate"))
        parpool("Processes", 4);
    end
    scan_tiles = cell(n_scan, 1);
    raw_files = scans.raw_file;
    antenna_agl = scans.antenna_agl;
    parfor i_s = 1:n_scan
        level2 = bmm4d.radar.read_level2_moments(raw_files(i_s), ["REF" "CFP"], max_elevation_deg);
        has_ref = arrayfun(@(s) ~isempty(s.REF), level2);
        % Split cuts repeat an elevation: keep the first (surveillance) cut; CFP may come with the other cut
        ref_cuts = find(has_ref);
        [~, keep] = unique(round([level2(ref_cuts).elev] * 2.5) / 2.5, "stable");
        rows = {};
        for i_e = ref_cuts(keep)
            n_gate = min(size(level2(i_e).REF, 1), find(level2(i_e).REF_range <= max_range_m, 1, "last"));
            range_m = level2(i_e).REF_range(1:n_gate);
            dbz = level2(i_e).REF(1:n_gate, :);
            az_step_deg = 360 / size(dbz, 2);
            % Gates that lost >= 3 dB to the clutter filter (NaN where the file has no CFP)
            cfp_lost = nan(size(dbz));
            i_cfp = find(abs([level2.elev] - level2(i_e).elev) < 0.2 & arrayfun(@(s) ~isempty(s.CFP), level2), 1);
            if ~isempty(i_cfp) && size(level2(i_cfp).CFP_code, 2) == size(dbz, 2)
                n_cfp = min(n_gate, size(level2(i_cfp).CFP_code, 1));
                cfp_lost(1:n_cfp, :) = level2(i_cfp).CFP_code(1:n_cfp, :) >= 11;
            end
            % Clutter and rain (> 3000 birds km^-3): no information
            excluded = z_to_refl(idb(dbz)) / 11 > max_density;
            output = z_to_refl(idb(dbz)) / 11;
            output(isnan(dbz)) = 0;
            output(excluded) = nan;
            detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
            ll = [];
            dens = [];
            for i_g = 1:n_sigma
                [lo, hi] = bmm4d.radar.gate.bounds(dbz, isnan(dbz), range_m, beams{i_g}, detection_dbz, s_floor_cm2);
                lo(excluded) = -Inf; hi(excluded) = Inf;
                [tiles, tile_id] = gd_estimate_tiles(lo, hi, range_m, az_step_deg, beams{i_g}, tables{i_g}, density_grid, tile_m);
                ll = [ll, tiles.log_likelihood];
                dens = [dens, tiles.density];
            end
            % Per tile: scan, elevation, antenna height, centre range, informative gates, output mean,
            % share of gates that lost >= 3 dB, then log likelihood and density for every log-sd
            n_tile = height(tiles);
            valid = ~isnan(output);
            output_mean = accumarray(tile_id(valid), output(valid), [n_tile 1]) ./ accumarray(tile_id(valid), 1, [n_tile 1]);
            lost = accumarray(tile_id(:), cfp_lost(:), [n_tile 1]) ./ accumarray(tile_id(:), 1, [n_tile 1]);
            rows{end + 1} = [repmat([i_s, level2(i_e).elev, antenna_agl(i_s)], n_tile, 1), ...
                tiles.center_km, tiles.n_known, output_mean, lost, ll, dens];
        end
        scan_tiles{i_s} = cat(1, rows{:});
        fprintf("[spread] scan %d done\n", i_s)
    end

    %% One row per tile -----------------------------------------------------------------
    all_tiles = cat(1, scan_tiles{:});
    tile = array2table(all_tiles(:, 1:7), VariableNames=["scan" "elevation" "antenna_agl" "range_km" "n_known" "output_mean" "cfp_lost"]);
    tile.height_m = tile.antenna_agl + tile.range_km * 1000 .* sind(tile.elevation) + (tile.range_km * 1000) .^ 2 / (2 * earth_radius_m);
    ll = all_tiles(:, 7 + (1:n_sigma));                               % [tile sigma]
    density = all_tiles(:, 7 + n_sigma + (1:n_sigma));
    keep = tile.n_known >= min_known;
    tile = tile(keep, :); ll = ll(keep, :); density = density(keep, :);
    fprintf("[spread] %d tiles\n", height(tile))
    % Expected birds per gate (counting volume) at the reference log-sd 1.5
    beam_reference = bmm4d.radar.gate.beam(sigma_log_sd=1.5);
    tile.lambda = density(:, sigma_grid == 1.5) .* beam_reference.volume_sample_km3(tile.range_km * 1000);
    save(tiles_file, "tile", "ll", "density", "sigma_grid", "scans", "-v7.3")
else
    load(tiles_file, "tile", "ll", "density")
end

%% Best log-sd per group -------------------------------------------------------------
% Profile over the grid, refined by a parabola; bootstrap over scans for the spread
best_of = @(profile) gd_profile_max(profile, sigma_grid);
groups = struct();
groups.range = discretize(tile.range_km, [2 5 10 15 20 30 45 60 80]);
groups.height = discretize(tile.height_m, [0 300 500 800 1200 1700 2500 4000]);
groups.lambda = discretize(tile.lambda, [0 0.3 1 3 10 30 100 300 1e5]);
groups.cfp = discretize(tile.cfp_lost, [0 0.05 0.2 0.5 1.0001]);
group_names = string(fieldnames(groups))';
summary = table();
rng(4)
for name = group_names
    g = groups.(name);
    for k = unique(g(~isnan(g)))'
        in = g == k;
        scans_in = unique(tile.scan(in));
        % Profile of each scan within the class, then resample scans
        scan_profile = zeros(n_scan, n_sigma);
        for i_g = 1:n_sigma
            scan_profile(:, i_g) = accumarray(tile.scan(in), ll(in, i_g), [n_scan 1]);
        end
        boot = nan(200, 1);
        for b = 1:200
            count = accumarray(scans_in(randi(numel(scans_in), numel(scans_in), 1)), 1, [n_scan 1])';
            boot(b) = best_of(count * scan_profile);
        end
        far = in & tile.range_km >= 20;                                % where gates hold many birds
        ratio = sum(density(far, :) .* tile.n_known(far), 1) ./ sum(tile.output_mean(far) .* tile.n_known(far), 1);
        summary = [summary; table(name, k, nnz(in), numel(scans_in), median(tile.range_km(in)), median(tile.height_m(in)), ...
            median(tile.lambda(in)), best_of(sum(ll(in, :), 1)), prctile(boot, 2.5), prctile(boot, 97.5), ...
            gd_crossing(sigma_grid, ratio, 1), ...
            VariableNames=["group" "class" "n_tile" "n_scan" "range_km" "height_m" "lambda" "best_sigma" ...
            "boot_2_5" "boot_97_5" "sigma_mean_match"])]; %#ok<AGROW>
    end
end
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_G_spread_diagnostic.csv"))

%% Range x height map of the best log-sd ---------------------------------------------------
range_edges = [2 5 10 15 20 30 45 60 80];
height_edges = [0 300 500 800 1200 1700 2500 4000];
best_map = nan(numel(height_edges) - 1, numel(range_edges) - 1);
for i_r = 1:numel(range_edges) - 1
    for i_h = 1:numel(height_edges) - 1
        in = groups.range == i_r & groups.height == i_h;
        if nnz(in) >= 200 && numel(unique(tile.scan(in))) >= 10
            best_map(i_h, i_r) = best_of(sum(ll(in, :), 1));
        end
    end
end
lambda_map = nan(size(best_map));
for i_r = 1:numel(range_edges) - 1
    for i_h = 1:numel(height_edges) - 1
        lambda_map(i_h, i_r) = median(tile.lambda(groups.range == i_r & groups.height == i_h), "omitnan");
    end
end
save(fullfile(figure_dir, "GD_G_spread_diagnostic.mat"), "best_map", "lambda_map", "range_edges", "height_edges", "summary", "sigma_grid")

%% Plot ------------------------------------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3000, 1150, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
ax = nexttile(t);
% Cells by range and height bin index; tick labels give the bin edges
imagesc(ax, best_map, AlphaData=~isnan(best_map)); set(ax, YDir="normal", Color="w", FontSize=12)
colormap(ax, parula(256)); clim(ax, [0.5 2.5]); cb = colorbar(ax); cb.Label.String = "Best s_\sigma";
xticks(ax, 0.5:1:numel(range_edges) - 0.5); xticklabels(ax, string(range_edges))
yticks(ax, 0.5:1:numel(height_edges) - 0.5); yticklabels(ax, string(height_edges))
xlabel(ax, "Range (km)"); ylabel(ax, "Beam centre height (m AGL)"); title(ax, "Best s_\sigma by range and height", FontWeight="normal")
for name = ["lambda" "cfp"]
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    rows = summary(summary.group == name, :);
    x = rows.lambda;
    if name == "cfp"
        x = (1:height(rows))';
    end
    errorbar(ax, x, rows.best_sigma, rows.best_sigma - rows.boot_2_5, rows.boot_97_5 - rows.best_sigma, "-o", ...
        Color="k", LineWidth=2, DisplayName="tile likelihood")
    plot(ax, x, rows.sigma_mean_match, "s--", Color=[0.1 0.4 0.8], LineWidth=1.5, DisplayName="estimate = output mean (tiles beyond 20 km)")
    set(ax, FontSize=12); ylim(ax, [0.5 2.5]); ylabel(ax, "Best s_\sigma")
    if name == "lambda"
        set(ax, XScale="log"); xlim(ax, [0.1 1000]); xticks(ax, [0.1 1 10 100 1000])
        xlabel(ax, "Birds per counting volume, \lambda (median of class)")
        title(ax, "By expected birds per gate", FontWeight="normal"); legend(ax, Location="south", FontSize=10)
    else
        xticks(ax, x); xticklabels(ax, ["<5%" "5-20%" "20-50%" ">50%"]); xlabel(ax, "Gates in the tile that lost \geq 3 dB to the clutter filter")
        title(ax, "By clutter filter (files with CFP)", FontWeight="normal")
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_G_spread_diagnostic.png"), Resolution=150)
close(fig)

