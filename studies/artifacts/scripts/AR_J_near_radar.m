%% Artefacts J: the reconstruction near the radars (dev)
% The full-night maps show patterns centred on the radars. Two questions:
% 1. Is the column near a radar higher or lower than further out? Column
%    against the distance to the nearest radar, each radar and frame
%    relative to its own 60-100 km ring (the azimuthal mean, so real
%    spatial structure averages out over ~100 radars), for the us_east
%    runs (us_east_fix, us_east_shared) and, with the prior, for the
%    Midwest case (ar_midwest_shared, production settings, white prior).
%    And: are the plumes stretched along the flight direction (column
%    within 30 deg of the track against 60-90 deg from it, 30-100 km)?
% 2. Do the data near a radar agree with the model? Observed over
%    predicted density of the Midwest solve by range and elevation, by
%    sweep at 5-15 km and by predicted density: the forward model holds
%    the beam, so a residual that depends on the elevation at a fixed
%    range is in the data (cleaning, scan mode, bird cross-section).
% 3. Which cuts make the dip (solves without the batch or the split cuts
%    within 20 km), and what in the raw reflectivity differs (detection,
%    the gate-density estimator, the clutter filter).
% Outputs: figures/AR_J_distance.csv, AR_J_residuals.csv, AR_J_sweeps.csv,
% AR_J_drop.csv, AR_J_raw_5_15km.csv, AR_J_near_radar.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
prod_dir = fullfile(bmm4d.util.local_data_dir(), "prod");
distance_edges = 0:10:150;      % km
n_bin = numel(distance_edges) - 1;
reference_bins = 7:10;          % 60-100 km

%% Midwest case: solve with production settings ------------------------------------
c = ar_test_case("ar_midwest_shared", [-100 -85 36 46], ...
    datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30"), "read_options", {"shape", false});
S = ar_load(c);
result = ar_solve(S, "shared function");
predicted = S.observations.F(:, S.active_id) * double(result.field);
observed = full(S.observations.values);
n_time = size(observed, 2);
elev = repmat(S.observations.row_elev_deg, 1, n_time);
range_km = repmat(S.observations.row_range_m / 1000, 1, n_time);
use = observed > 0.5 & predicted > 0.5;

%% Residuals by range and elevation -------------------------------------------------
range_edges = [5 10 15 20 30 40 60 90 150];
elev_edges = [0 0.7 1.1 1.6 2.1 3 4.5 7];
i_range = discretize(range_km, range_edges);
i_elev = discretize(elev, elev_edges);
ok = use & ~isnan(i_range) & ~isnan(i_elev);
sz = [numel(range_edges) - 1, numel(elev_edges) - 1];
ratio = accumarray([i_range(ok) i_elev(ok)], observed(ok), sz) ./ accumarray([i_range(ok) i_elev(ok)], predicted(ok), sz);
ratio(accumarray([i_range(ok) i_elev(ok)], 1, sz) < 2000) = NaN;
residuals = array2table(round(ratio, 3), ...
    RowNames=compose("%d-%d km", range_edges(1:end-1)', range_edges(2:end)'), ...
    VariableNames=compose("elev %.1f-%.1f", elev_edges(1:end-1)', elev_edges(2:end)'));
disp(residuals)
writetable(residuals, fullfile(figure_dir, "AR_J_residuals.csv"), WriteRowNames=true)

% By sweep at 5-15 km (VCP elevations: split cuts up to 1.3 deg, batch cuts above)
near = use & range_km >= 5 & range_km < 15;
sweep_elev = unique(round(elev(near), 1));
sweeps = table();
for e = sweep_elev(:)'
    sel = near & abs(elev - e) < 0.05;
    if nnz(sel) > 2000
        sweeps = [sweeps; table(e, nnz(sel), sum(observed(sel)) / sum(predicted(sel)), ...
            VariableNames=["elev_deg" "rows" "observed_over_predicted"])]; %#ok<AGROW>
    end
end
disp(sweeps)
writetable(sweeps, fullfile(figure_dir, "AR_J_sweeps.csv"))

% Split against batch cuts at 5-15 km by predicted density quartile
q = quantile(predicted(near & elev >= 1.6), 0:0.25:1);
for k = 1:4
    in_q = near & predicted >= q(k) & predicted < q(k + 1);
    fprintf("predicted %5.1f-%5.1f birds km^-3: obs/pred %.2f at <= 1.3 deg, %.2f at >= 1.8 deg\n", q(k), q(k + 1), ...
        sum(observed(in_q & elev < 1.6)) / sum(predicted(in_q & elev < 1.6)), ...
        sum(observed(in_q & elev >= 1.6)) / sum(predicted(in_q & elev >= 1.6)))
end

%% Column against the distance to the nearest radar ----------------------------------
% Midwest: prior and reconstruction (columns weighted by each cell's air share)
[dist_km, nearest] = nearest_radar(S.g, S.radar_selected);
n_xy = prod(S.g.sz(1:2));
xy = mod(S.active_id - 1, n_xy) + 1;
distance = table(distance_edges(1:end-1)' + 5, VariableNames="distance_km");
fields = {S.prior, result.field};
names = ["midwest_prior" "midwest_reconstruction"];
for i_f = 1:2
    columns = zeros(n_xy, n_time);
    for t = 1:n_time
        columns(:, t) = accumarray(xy, double(fields{i_f}(:, t)) .* S.cell_air, [n_xy 1]);
    end
    distance.(names(i_f)) = distance_profile(columns, xy, n_xy, dist_km, nearest, distance_edges, reference_bins);
end

% us_east runs
us = load(fullfile(prod_dir, "us_east_fix", "case.mat"), "g", "radar_selected", "radar_scans");
[dist_us, nearest_us, bearing_us] = nearest_radar(us.g, us.radar_selected);
n_xy_us = prod(us.g.sz(1:2));
runs = ["us_east_fix" "us_east_shared"];
along_cross = nan(height(us.radar_selected), numel(us.radar_scans(1, :)), numel(runs));
for i_run = 1:numel(runs)
    run_field = load(fullfile(prod_dir, runs(i_run), "inversion_fix.mat"), "active_field", "active_id");
    xy_us = mod(run_field.active_id - 1, n_xy_us) + 1;
    air = bmm4d.grid.air_fraction(us.g, run_field.active_id);
    columns = zeros(n_xy_us, size(run_field.active_field, 2));
    for t = 1:size(columns, 2)
        columns(:, t) = accumarray(xy_us, double(run_field.active_field(:, t)) .* air, [n_xy_us 1]);
    end
    clear run_field
    distance.(runs(i_run)) = distance_profile(columns, xy_us, n_xy_us, dist_us, nearest_us, distance_edges, reference_bins);

    % Along the track (within 30 deg of it or its reverse) against across (60-90 deg), 30-100 km
    has = accumarray(xy_us, 1, [n_xy_us 1]) > 0;
    ring = find(has & dist_us > 30 & dist_us < 100);
    for t = 1:size(columns, 2)
        for k = 1:height(us.radar_selected)
            scan = us.radar_scans{k, t};
            if isempty(scan), continue, end
            uv = mean(scan.uv(all(isfinite(scan.uv), 2), :), 1);
            if any(~isfinite(uv)) || hypot(uv(1), uv(2)) < 3, continue, end
            mine = ring(nearest_us(ring) == k);
            cos_psi = abs(cosd(bearing_us(mine) - atan2d(uv(1), uv(2))));
            along_cross(k, t, i_run) = mean(columns(mine(cos_psi > cosd(30)), t)) / mean(columns(mine(cos_psi < cosd(60)), t));
        end
    end
    v = along_cross(:, :, i_run);
    v = v(isfinite(v) & v > 0);
    fprintf("%s: along/across-track column at 30-100 km, median %.3f (%d radar-frames)\n", runs(i_run), median(v), numel(v))
end
disp(distance)
writetable(distance, fullfile(figure_dir, "AR_J_distance.csv"))

%% Which cuts make the dip: solve without the batch or the split cuts within 20 km
% Row weight 0; the dropped rows are then predicted from the others.
drops = struct("name", {"no batch cuts < 20 km", "no split cuts < 20 km"}, ...
    "rows", {range_km < 20 & elev >= 1.6, range_km < 20 & elev < 1.6});
drop_profiles = table(distance.distance_km, distance.midwest_reconstruction, VariableNames=["distance_km" "all_rows"]);
for i_d = 1:numel(drops)
    S_drop = S;
    S_drop.row_weights(drops(i_d).rows) = 0;
    result_drop = ar_solve(S_drop, drops(i_d).name);
    columns = zeros(n_xy, n_time);
    for t = 1:n_time
        columns(:, t) = accumarray(xy, double(result_drop.field(:, t)) .* S.cell_air, [n_xy 1]);
    end
    drop_profiles.(matlab.lang.makeValidName(drops(i_d).name)) = distance_profile(columns, xy, n_xy, dist_km, nearest, distance_edges, reference_bins);
    predicted_drop = S.observations.F(:, S.active_id) * double(result_drop.field);
    dropped = drops(i_d).rows & use;
    fprintf("%s: dropped rows observed / predicted %.2f\n", drops(i_d).name, sum(observed(dropped)) / sum(predicted_drop(dropped)))
    clear S_drop result_drop predicted_drop
end
disp(drop_profiles(1:6, :))
writetable(drop_profiles, fullfile(figure_dir, "AR_J_drop.csv"))

%% Raw reflectivity near the radar: detection, the density estimator, the clutter filter
% Raw Archive II cuts of 12 Midwest radars at 04-07 UTC (one volume per
% hour), against the cleaned volumes (us_east_shared cache, aspect undone).
% - echo share and range-normalised 1st percentile of dBZ per cut at 5-15 km;
% - estimated density per unit linear Z on the kept (thinned) gates,
%   pooled over volumes, per sweep: the same constant converts Z to
%   density for every sweep, so the ratio should not depend on the sweep;
% - clutter filter power removed (CFP; codes: 0 below threshold, 1 range
%   folded, 2 filter not applied, >= 8: dB = code - 8): share of gates the
%   filter removed power from and their mean Z over that of the others,
%   same range and cut. Above ~1.8 deg the main lobe is far above the
%   ground, so a lower Z there is bird echo removed by the filter, or the
%   filter decision leaving the strong bird echo unfiltered.
addpath(fullfile(project_root, "studies", "aspect", "functions"))   % as_read_volume
raw_dir = bmm4d.radar.nexrad_raw_dir();
volume_dir = fullfile(prod_dir, "volumes", "us_east_shared");
raw_radars = ["KDMX" "KDVN" "KOAX" "KEAX" "KARX" "KLSX" "KMKX" "KICT" "KTWX" "KSGF" "KILX" "KLOT"];
raw_rows = table();
for radar_name = raw_radars
    volumes = dir(fullfile(volume_dir, radar_name + "20261022_*.mat"));
    if isempty(volumes), continue, end
    stamp = datetime(extractBetween(string({volumes.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
    for target = datetime("2026-10-22 04:00"):hours(1):datetime("2026-10-22 07:00")
        [gap, i_v] = min(abs(stamp - target));
        raw_file = fullfile(raw_dir, radar_name + "_20261022", erase(volumes(i_v).name, ".mat"));
        if gap > minutes(10) || ~isfile(raw_file), continue, end
        volume = as_read_volume(fullfile(volumes(i_v).folder, volumes(i_v).name), 3000);
        cuts = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], 7);
        for k = 1:numel(volume.dens.sweeps)
            sweep = volume.dens.sweeps(k);
            % Raw cut of this sweep: same number of radials, closest elevation, longest range (surveillance cut)
            candidates = find(arrayfun(@(cut) numel(cut.az) == numel(sweep.az) && ~isempty(cut.REF) && ~isempty(cut.CFP), cuts));
            if isempty(candidates), continue, end
            [~, best] = min(abs([cuts(candidates).elev] - sweep.elev) - 1e-6 * arrayfun(@(cut) max(cut.REF_range), cuts(candidates)));
            cut = cuts(candidates(best));
            n_gate = min(numel(cut.REF_range), numel(cut.CFP_range));
            z_lin = 10 .^ (cut.REF(1:n_gate, :) / 10);
            z_lin(cut.REF_code(1:n_gate, :) == 0) = 0;
            z_lin(cut.REF_code(1:n_gate, :) == 1) = NaN;
            filter_code = double(cut.CFP_code(1:n_gate, :));
            gate_range = cut.REF_range(1:n_gate);
            % Raw gate of each kept (thinned) gate of the cleaned sweep
            [~, i_az] = min(abs(mod(sweep.az(:)' - cut.az(:) + 180, 360) - 180), [], 1);
            [~, i_rg] = min(abs(sweep.range(:) - gate_range(:)'), [], 2);
            in_band = sweep.range(:) >= 5000 & sweep.range(:) <= 15000;
            density = sweep.data(in_band, :);
            z_kept = z_lin(i_rg(in_band), i_az);
            ok = isfinite(density) & isfinite(z_kept);
            near_gates = gate_range >= 5000 & gate_range <= 15000;
            z = z_lin(near_gates, :);
            code = filter_code(near_gates, :);
            filtered = code >= 9 & isfinite(z);
            unfiltered = (code == 8 | code == 2) & isfinite(z);
            normalised_dbz = 10 * log10(z) - 20 * log10(gate_range(near_gates) / 1000);   % dBZ at 1 km
            echo = z > 0;
            raw_rows = [raw_rows; table(radar_name, target, round(sweep.elev, 1), numel(sweep.az), mean(echo(isfinite(z))), ...
                prctile(normalised_dbz(echo), 1), sum(density(ok)), sum(z_kept(ok)), ...
                nnz(filtered), nnz(unfiltered), sum(z(filtered)), sum(z(unfiltered)), ...
                VariableNames=["radar" "time" "elev" "n_az" "echo_share" "normalised_dbz_p1" "sum_density" "sum_z" ...
                "n_filtered" "n_unfiltered" "z_filtered" "z_unfiltered"])]; %#ok<AGROW>
        end
    end
end
raw_summary = groupsummary(raw_rows, "elev", {"median", "sum"}, ["echo_share" "normalised_dbz_p1" "sum_density" "sum_z" ...
    "n_filtered" "n_unfiltered" "z_filtered" "z_unfiltered"]);
raw_summary = raw_summary(raw_summary.GroupCount >= 20, :);
density_per_z = raw_summary.sum_sum_density ./ raw_summary.sum_sum_z;
raw_summary.density_per_z_vs_1p3 = density_per_z / density_per_z(raw_summary.elev == 1.3);
raw_summary.share_filtered = raw_summary.sum_n_filtered ./ (raw_summary.sum_n_filtered + raw_summary.sum_n_unfiltered);
raw_summary.z_filtered_over_unfiltered = (raw_summary.sum_z_filtered ./ raw_summary.sum_n_filtered) ./ ...
    (raw_summary.sum_z_unfiltered ./ raw_summary.sum_n_unfiltered);
raw_summary = raw_summary(:, ["elev" "GroupCount" "median_echo_share" "median_normalised_dbz_p1" "density_per_z_vs_1p3" ...
    "share_filtered" "z_filtered_over_unfiltered"]);
disp(raw_summary)
writetable(raw_summary, fullfile(figure_dir, "AR_J_raw_5_15km.csv"))

%% Plot results ---------------------------------------------------------------------
fig = bmm4d.plot.create_export_figure(2000, 800);
set(fig, "DefaultAxesFontName", "Helvetica", "DefaultTextFontName", "Helvetica")
layout = tiledlayout(fig, 1, 2, "TileSpacing", "tight", "Padding", "tight");

ax = nexttile(layout);
hold(ax, "on")
plot(ax, distance.distance_km, distance.us_east_fix, "-o", Color=[0.85 0.33 0.10], MarkerSize=3, DisplayName="us\_east, per-volume fit")
plot(ax, distance.distance_km, distance.us_east_shared, "-s", Color=[0.55 0.10 0.05], MarkerSize=3, DisplayName="us\_east, shared function")
plot(ax, distance.distance_km, distance.midwest_reconstruction, "-^", Color=[0.1 0.45 0.6], MarkerSize=3, DisplayName="Midwest, reconstruction")
plot(ax, distance.distance_km, distance.midwest_prior, "--", Color=[0.5 0.5 0.5], DisplayName="Midwest, prior")
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Distance to the nearest radar (km)")
ylabel(ax, "Column / column at 60-100 km")
legend(ax, Location="southeast", Box="off")
title(ax, "Column near the radars", FontWeight="bold")
grid(ax, "on")

ax = nexttile(layout);
hold(ax, "on")
range_mid = (range_edges(1:end-1) + range_edges(2:end)) / 2;
colors = copper(5);
show = [1 2 3 5 7];
for i_show = 1:numel(show)
    plot(ax, (elev_edges(1:end-1) + elev_edges(2:end)) / 2, ratio(show(i_show), :), "-o", Color=colors(i_show, :), MarkerSize=3, ...
        DisplayName=sprintf("%d-%d km", range_edges(show(i_show)), range_edges(show(i_show) + 1)))
end
xline(ax, 1.55, ":", "split | batch cuts", Color=[0.4 0.4 0.4], FontSize=6, HandleVisibility="off")
yline(ax, 1, ":", Color=[0.4 0.4 0.4], HandleVisibility="off")
xlabel(ax, "Elevation (deg)")
ylabel(ax, "Observed / predicted")
legend(ax, Location="southwest", Box="off")
title(ax, "Midwest residuals by range", FontWeight="bold")
grid(ax, "on")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_J_near_radar.png"), Resolution=300)
close(fig)

%% Local functions ------------------------------------------------------------------
function [dist_km, nearest, bearing_deg] = nearest_radar(g, radars)
% Distance (km, local tangent plane) and compass bearing from the nearest
% radar to each grid column, [x y] order
[lon, lat] = ndgrid(double(g.lon(:)), double(g.lat(:)));
lon = lon(:);
lat = lat(:);
dist_km = inf(numel(lon), 1);
nearest = zeros(numel(lon), 1);
for k = 1:height(radars)
    d = 6371 * hypot(deg2rad(lon - radars.lon(k)) .* cosd(radars.lat(k)), deg2rad(lat - radars.lat(k)));
    closer = d < dist_km;
    dist_km(closer) = d(closer);
    nearest(closer) = k;
end
bearing_deg = atan2d(deg2rad(lon - radars.lon(nearest)) .* cosd(radars.lat(nearest)), deg2rad(lat - radars.lat(nearest)));
end

function profile = distance_profile(columns, xy, n_xy, dist_km, nearest, edges, reference_bins)
% Median over radars and frames of the mean column per distance bin,
% each radar and frame relative to its reference bins
bin = discretize(dist_km, edges);
use = accumarray(xy, 1, [n_xy 1]) > 0 & ~isnan(bin);
n_radar = max(nearest);
n_bin = numel(edges) - 1;
relative = nan(n_radar, size(columns, 2), n_bin);
for t = 1:size(columns, 2)
    m = accumarray([nearest(use) bin(use)], columns(use, t), [n_radar n_bin], @mean, NaN);
    relative(:, t, :) = m ./ mean(m(:, reference_bins), 2, "omitnan");
end
relative = reshape(relative, [], n_bin);
complete = all(isfinite(relative) & relative > 0, 2);
profile = median(relative(complete, :), 1)';
end
