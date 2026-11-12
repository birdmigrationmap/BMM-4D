%% Sandbox: does beam propagation explain the radar level differences?
% SB_overlap found that neighbouring radars differ by 2-3 dB on matched
% gates (same height and range under the 4/3-Earth model), mostly as one
% offset per radar that changes from night to night and is a third as
% large on the 0.9 deg cut. Test, without fitting: if the real beams bend
% differently around each radar, the matched gates are not at the same
% height, and the bird profile turns that height difference into dB.
% 1. Refractivity N(h) = 77.6 P/T + 3.73e5 e/T^2 from radiosondes (Iowa
%    Environmental Mesonet archive) at 00 and 12 UTC around each night,
%    interpolated to each radar (inverse distance, 3 nearest sites) and in
%    time to the volume time.
% 2. Beam centre height by ray tracing in a layered spherical atmosphere
%    (dtheta/ds = cos(theta) (1/(R+h) + dn/dh / n)), against 4/3 Earth.
% 3. Bird profile (linear Z by 100 m height) from each radar's own volume
%    at 5-35 km, where bending hardly matters; the air at a pair's matched
%    gates gets the mean of the two radars' profiles.
% 4. Predicted pair difference = beam-weighted profile at A's true height
%    minus the same at B's (two-way Gaussian beam), against the measured
%    median difference of SB_overlap (figures/overlap_pairs.csv).
% Outputs: local_data_dir/sandbox_nights/sondes (json), figures/propagation.png,
% propagation_pairs.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
sonde_dir = fullfile(bmm4d.util.local_data_dir(), "sandbox_nights", "sondes");
if ~isfolder(sonde_dir)
    mkdir(sonde_dir);
end
c = case_settings("sandbox");
radar = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
radar = radar(ismember(radar.name, c.radars), :);
n_radar = height(radar);
files = readtable(fullfile(figure_dir, "night_maps_files.csv"), TextType="string");
files.raw_file(ismissing(files.raw_file)) = "";
pairs = readtable(fullfile(figure_dir, "overlap_pairs.csv"), TextType="string");
nights = unique(files.night);
n_night = numel(nights);
earth_radius_m = 6371000;
effective_radius_m = 4 / 3 * earth_radius_m;
sigma_beam_rad = deg2rad(0.925 / (2 * sqrt(2 * log(2)))) / sqrt(2);   % two-way, as compute_F_elev
height_grid = (0:10:6000)';                                          % m a.s.l.

%% 1. Radiosonde refractivity profiles ------------------------------------------
sites = table(["KABR" "KBIS" "KMPX" "KGRB" "KAPX" "KDTX" "KDVN" "KOAX" "KLBF" "KILX" "KILN" "KTOP" "KDDC" ...
    "KSGF" "KBNA" "KOUN" "KLZK" "KAMA"]', ...
    [45.46 46.77 44.85 44.48 44.91 42.70 41.61 41.32 41.13 40.15 39.42 39.07 37.76 37.23 36.25 35.18 34.84 35.23]', ...
    [-98.41 -100.76 -93.57 -88.13 -84.72 -83.47 -90.58 -96.37 -100.68 -89.34 -83.82 -95.63 -99.97 -93.40 -86.56 ...
    -97.44 -92.26 -101.71]', VariableNames=["name" "lat" "lon"]);
% Volume time of each night (all radars within 5 min of it); soundings at
% 00 and 12 UTC of that UTC day
volume_time = arrayfun(@(d) files.target_time(find(files.night == d, 1)), nights);
sounding_time = [dateshift(volume_time, "start", "day"), dateshift(volume_time, "start", "day") + hours(12)];
refractivity = nan(numel(height_grid), height(sites), n_night, 2);
for i_n = 1:n_night
    for i_s = 1:2
        stamp = string(sounding_time(i_n, i_s), "yyyyMMddHHmm");
        for i_site = 1:height(sites)
            file = fullfile(sonde_dir, sites.name(i_site) + "_" + stamp + ".json");
            if ~isfile(file)
                try
                    websave(file, "https://mesonet.agron.iastate.edu/json/raob.py?ts=" + stamp + "&station=" + sites.name(i_site));
                catch
                    continue
                end
            end
            refractivity(:, i_site, i_n, i_s) = sounding_refractivity(file, height_grid);
        end
    end
end
fprintf("%d of %d soundings usable\n", nnz(any(~isnan(refractivity), 1)), height(sites) * n_night * 2)

% At each radar and volume time: inverse distance (3 nearest sites with a
% sounding, power 2), then linear in time between 00 and 12 UTC
site_km = 6371 * acos(min(1, sind(radar.lat) .* sind(sites.lat') + cosd(radar.lat) .* cosd(sites.lat') .* cosd(radar.lon - sites.lon')));
radar_refractivity = nan(numel(height_grid), n_radar, n_night);
for i_n = 1:n_night
    time_weight = hours(volume_time(i_n) - sounding_time(i_n, 1)) / 12;
    for i_r = 1:n_radar
        both_times = nan(numel(height_grid), 2);
        for i_s = 1:2
            available = find(all(~isnan(refractivity(:, :, i_n, i_s)), 1));
            [d, order] = sort(site_km(i_r, available));
            use = available(order(1:min(3, numel(order))));
            w = 1 ./ max(d(1:numel(use)), 10) .^ 2;
            both_times(:, i_s) = refractivity(:, use, i_n, i_s) * (w(:) / sum(w));
        end
        radar_refractivity(:, i_r, i_n) = both_times * [1 - time_weight; time_weight];
    end
end

%% 2. Check the ray tracer on the standard atmosphere ----------------------------
standard = 315 * exp(-height_grid / 7000);             % dN/dh ~ -40 km^-1 near the ground
test_range = (0:1000:150000)';
ray = ray_height(height_grid, standard, 300, 0.5, test_range);
four_thirds = 300 + test_range * sind(0.5) + test_range .^ 2 / (2 * effective_radius_m);
fprintf("ray tracer, standard atmosphere: height - 4/3 Earth at 150 km = %.0f m\n", ray(end) - four_thirds(end))

%% 3. Bird profiles near each radar ----------------------------------------------
profile_edges = 0:100:5000;
bird_profile = nan(numel(profile_edges) - 1, n_radar, n_night);
profile_sum = cell(height(files), 1);
parfor k = 1:height(files)
    if files.raw_file(k) == ""
        continue
    end
    try
        sweeps = bmm4d.radar.read_level2_moments(files.raw_file(k), ["REF" "RHO"], 5);
    catch
        continue
    end
    [~, i_r] = ismember(files.station(k), radar.name);
    z_sum = zeros(numel(profile_edges) - 1, 1);
    count = zeros(numel(profile_edges) - 1, 1);
    sweeps = sweeps(arrayfun(@(x) ~isempty(x.REF) && ~isempty(x.RHO), sweeps));
    for i_e = 1:numel(sweeps)
        s = sweeps(i_e);
        n_common = min(size(s.REF, 1), size(s.RHO, 1));
        range_m = double(s.REF_range(1:n_common));
        in_range = range_m >= 5000 & range_m <= 35000;
        ref = s.REF(in_range, :);
        rain = ref >= 20 & s.RHO(in_range, :) >= 0.97;
        z = 10 .^ (ref / 10);
        z(isnan(ref)) = 0;                              % below threshold: no echo
        keep = ~rain & ~(ref > 35);
        gate_height = radar.height(i_r) + range_m(in_range) * sind(s.elev) + range_m(in_range) .^ 2 / (2 * effective_radius_m);
        bin = discretize(repmat(gate_height, 1, size(ref, 2)), profile_edges);
        ok = keep & ~isnan(bin);
        z_sum = z_sum + accumarray(bin(ok), z(ok), size(z_sum));
        count = count + accumarray(bin(ok), 1, size(count));
    end
    profile_sum{k} = z_sum ./ count;
end
for k = 1:height(files)
    if ~isempty(profile_sum{k})
        [~, i_r] = ismember(files.station(k), radar.name);
        bird_profile(:, i_r, nights == files.night(k)) = profile_sum{k};
    end
end
profile_center = profile_edges(1:end - 1)' + 50;

%% 4. Predicted pair differences --------------------------------------------------
pairs.predicted_db = nan(height(pairs), 1);
pairs.delta_h_a = nan(height(pairs), 1);
pairs.delta_h_b = nan(height(pairs), 1);
for i_p = 1:height(pairs)
    i_n = find(nights == pairs.night(i_p));
    a = pairs.id_a(i_p);
    b = pairs.id_b(i_p);
    % The air at the matched gates: mean of the two radars' profiles
    air = mean([bird_profile(:, a, i_n) bird_profile(:, b, i_n)], 2, "omitnan");
    if all(isnan(air)) || any(isnan(radar_refractivity(:, [a b], i_n)), "all")
        continue
    end
    air = fillmissing(air, "nearest");
    air_fine = interp1(profile_center, air, height_grid, "linear", "extrap");
    air_fine = max(air_fine, 0);
    seen = zeros(1, 2);
    sides = [a b];
    ranges = [pairs.range_a_m(i_p) pairs.range_b_m(i_p)];
    elevations = [pairs.elev_a(i_p) pairs.elev_b(i_p)];
    delta_h = zeros(1, 2);
    for i_side = 1:2
        r = sides(i_side);
        true_height = ray_height(height_grid, radar_refractivity(:, r, i_n), radar.height(r), elevations(i_side), ranges(i_side));
        model_height = radar.height(r) + ranges(i_side) * sind(elevations(i_side)) + ranges(i_side) ^ 2 / (2 * effective_radius_m);
        delta_h(i_side) = true_height - model_height;
        weight = exp(-(height_grid - true_height) .^ 2 / (2 * (ranges(i_side) * sigma_beam_rad) ^ 2));
        seen(i_side) = 10 * log10(sum(air_fine .* weight) / sum(weight));
    end
    pairs.predicted_db(i_p) = seen(1) - seen(2);
    pairs.delta_h_a(i_p) = delta_h(1);
    pairs.delta_h_b(i_p) = delta_h(2);
end
writetable(pairs, fullfile(figure_dir, "propagation_pairs.csv"))

cut_names = ["0.5°" "0.9°"];
for i_c = 1:2
    use = pairs.cut == i_c & ~isnan(pairs.predicted_db) & isfinite(pairs.predicted_db);
    m = pairs.median_db(use);
    p = pairs.predicted_db(use);
    slope = p \ m;
    fprintf("%s cut, %d pairs: measured sd %.2f dB, predicted sd %.2f dB, r = %.2f, slope %.2f, variance explained %.2f\n", ...
        cut_names(i_c), nnz(use), std(m), std(p), corr(m, p), slope, 1 - var(m - p) / var(m))
end
fprintf("beam height - 4/3 Earth at the matched gates: median %.0f m, 10-90%% %.0f to %.0f m\n", ...
    median([pairs.delta_h_a; pairs.delta_h_b], "omitnan"), prctile([pairs.delta_h_a; pairs.delta_h_b], [10 90]))

%% Plot results -----------------------------------------------------------
fig = bmm4d.plot.create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "Does beam bending (radiosonde refractivity) explain the level differences between neighbouring radars?", FontSize=14)
cut_colors = [0.1 0.3 0.6; 0.85 0.33 0.1];
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on"); axis(ax, "equal")
for i_c = 1:2
    use = pairs.cut == i_c & isfinite(pairs.predicted_db);
    scatter(ax, pairs.predicted_db(use), pairs.median_db(use), 18, cut_colors(i_c, :), "filled", MarkerFaceAlpha=0.5, ...
        DisplayName=sprintf("%s cut, r = %.2f", cut_names(i_c), corr(pairs.predicted_db(use), pairs.median_db(use))))
end
plot(ax, [-15 15], [-15 15], "k--", HandleVisibility="off")
xlim(ax, [-12 12]); ylim(ax, [-12 12]); set(ax, FontSize=12)
xlabel(ax, "Predicted from beam bending (dB)"); ylabel(ax, "Measured, matched gates (dB)")
legend(ax, Location="northwest"); title(ax, "Pair differences (A - B)", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
histogram(ax, [pairs.delta_h_a(pairs.cut == 1); pairs.delta_h_b(pairs.cut == 1)], -600:25:400, FaceColor=cut_colors(1, :))
xline(ax, 0, "k--"); set(ax, FontSize=12)
xlabel(ax, "Beam height - 4/3 Earth at the matched gates (m)"); ylabel(ax, "Radar sides of pairs, 0.5° cut")
title(ax, "How far the beams bend from the model", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
gradient_low = squeeze(mean(diff(radar_refractivity(height_grid >= 0 & height_grid <= 3000, :, :)), 1)) * 100;   % N per km over the radar layer
for i_n = 1:n_night
    agl = height_grid - mean(radar.height);
    plot(ax, mean(gradient(squeeze(radar_refractivity(:, :, i_n))')', 2) * 100, height_grid, LineWidth=1.2, ...
        DisplayName=string(nights(i_n), "yyyy-MM-dd"))
end
xline(ax, -40, "k--", "4/3 Earth", HandleVisibility="off"); xline(ax, -157, "k:", "ducting", HandleVisibility="off")
ylim(ax, [0 3000]); xlim(ax, [-250 50]); set(ax, FontSize=12)
xlabel(ax, "dN/dh (km^{-1}), mean over radars"); ylabel(ax, "Height (m a.s.l.)")
legend(ax, Location="southeast", FontSize=9); title(ax, "Refractivity gradient at volume time", FontWeight="normal")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "propagation.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function N = sounding_refractivity(file, height_grid)
% Refractivity on HEIGHT_GRID (m a.s.l.) from one IEM sounding json; NaN
% if the sounding is missing. Below and above the sounding: -40 km^-1.
N = nan(size(height_grid));
data = jsondecode(fileread(file));
if isempty(data.profiles) || isempty(data.profiles(1).profile)
    return
end
levels = data.profiles(1).profile;
if ~isstruct(levels)
    return
end
values = @(name) arrayfun(@(l) double_or_nan(l.(name)), levels);
pressure = values("pres");
height_m = values("hght");
temperature = values("tmpc") + 273.15;
dewpoint = values("dwpc");
ok = ~isnan(pressure) & ~isnan(height_m) & ~isnan(temperature) & ~isnan(dewpoint);
if nnz(ok) < 5
    return
end
[height_m, order] = unique(height_m(ok));
pressure = pressure(ok); pressure = pressure(order);
temperature = temperature(ok); temperature = temperature(order);
dewpoint = dewpoint(ok); dewpoint = dewpoint(order);
vapour = 6.112 * exp(17.67 * dewpoint ./ (dewpoint + 243.5));          % hPa
level_N = 77.6 * pressure ./ temperature + 3.73e5 * vapour ./ temperature .^ 2;
N = interp1(height_m, level_N, height_grid, "linear");
below = height_grid < height_m(1);
above = height_grid > height_m(end);
N(below) = level_N(1) - 0.04 * (height_grid(below) - height_m(1));
N(above) = level_N(end) - 0.04 * (height_grid(above) - height_m(end));
end

function value = double_or_nan(x)
if isempty(x)
    value = nan;
else
    value = double(x);
end
end

function h = ray_height(height_grid, N, h0, elev_deg, range_out)
% Beam centre height (m a.s.l.) at slant ranges RANGE_OUT (m) of a ray
% leaving h0 at ELEV_DEG through refractivity N(height_grid), by Euler
% steps of 50 m along the ray in a layered spherical atmosphere.
R = 6371000;
ds = 50;
n = 1 + N * 1e-6;
dn_dh = gradient(n, height_grid);
n_step = ceil(max(range_out) / ds);
h_path = zeros(n_step + 1, 1);
h = h0;
theta = deg2rad(elev_deg);
h_path(1) = h;
for k = 1:n_step
    i_h = min(max(round((h - height_grid(1)) / 10) + 1, 1), numel(height_grid));
    dtheta = cos(theta) * (1 / (R + h) + dn_dh(i_h) / n(i_h));
    h = h + sin(theta) * ds;
    theta = theta + dtheta * ds;
    h_path(k + 1) = h;
end
h = interp1((0:n_step)' * ds, h_path, range_out);
end
