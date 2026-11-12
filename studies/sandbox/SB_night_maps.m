%% Sandbox: lowest-sweep mosaics of the dry nights
% For each night of SB_night_screen (figures/night_screen.csv), downloads
% one NEXRAD Level II volume per sandbox radar, the one closest to 3 h
% after sunset at the domain centre (one UTC time per night, so the radars
% are compared at the same moment), and maps the reflectivity of the
% lowest cut. Each map pixel shows its nearest radar within 150 km, so
% level steps between radars appear as edges between their areas. The
% data are raw (no rain or clutter removal).
% Outputs: local_data_dir/sandbox_nights/raw (volumes),
% figures/night_maps.png and night_maps_files.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
raw_dir = fullfile(bmm4d.util.local_data_dir(), "sandbox_nights", "raw");
if ~isfolder(raw_dir)
    mkdir(raw_dir);
end
c = case_settings("sandbox");
radar = readtable(bmm4d.util.data_path("reference", "radar.csv"), TextType="string");
radar = radar(ismember(radar.name, c.radars), :);
selected = readtable(fullfile(figure_dir, "night_screen.csv"), TextType="string");
n_night = height(selected);
n_radar = height(radar);

%% Download one volume per radar and night ---------------------------------
% Target time: 3 h after sunset at the domain centre, for the night that
% starts on the local date
centre = [mean(c.figure_bounds(3:4)) mean(c.figure_bounds(1:2))];
target_time = NaT(n_night, 1);
for i_n = 1:n_night
    sunset = bmm4d.util.getSunriseSunset(centre(1), centre(2), datetime(selected.date(i_n)) + hours(24 - centre(2) / 15));
    target_time(i_n) = sunset + hours(3);
end
[i_night, i_radar] = ndgrid(1:n_night, 1:n_radar);
stations = radar.name(i_radar(:));
times = target_time(i_night(:));
raw_files = strings(numel(stations), 1);
parfor k = 1:numel(stations)
    try
        raw_files(k) = gd_download_level2(stations(k), times(k), raw_dir);
    catch err
        fprintf("%s %s: %s\n", stations(k), string(times(k)), err.message)
    end
end
files = table(selected.date(i_night(:)), stations, times, raw_files, VariableNames=["night" "station" "target_time" "raw_file"]);
writetable(files, fullfile(figure_dir, "night_maps_files.csv"))
fprintf("%d of %d volumes\n", nnz(raw_files ~= ""), numel(raw_files))

%% Read the lowest cut ---------------------------------------------------------
lowest = cell(height(files), 1);
parfor k = 1:height(files)
    if files.raw_file(k) == ""
        continue
    end
    try
        sweeps = bmm4d.radar.read_level2_moments(files.raw_file(k), "REF", 0.7);
        sweeps = sweeps(arrayfun(@(s) ~isempty(s.REF), sweeps));
        lowest{k} = struct("REF", sweeps(1).REF, "range", sweeps(1).REF_range, "az", sweeps(1).az, "elev", sweeps(1).elev);
    catch err
        fprintf("%s: %s\n", files.raw_file(k), err.message)
    end
end

%% Plot results -----------------------------------------------------------
% Mosaic on a 0.025 deg grid: nearest radar within 150 km, nearest gate
lon = c.figure_bounds(1):0.025:c.figure_bounds(2);
lat = c.figure_bounds(3):0.025:c.figure_bounds(4);
[lon_grid, lat_grid] = meshgrid(lon, lat);
distance_km = inf([size(lon_grid) n_radar]);
azimuth_deg = zeros([size(lon_grid) n_radar]);
for i_r = 1:n_radar
    distance_km(:, :, i_r) = 6371 * acos(min(1, sind(lat_grid) .* sind(radar.lat(i_r)) + ...
        cosd(lat_grid) .* cosd(radar.lat(i_r)) .* cosd(lon_grid - radar.lon(i_r))));
    azimuth_deg(:, :, i_r) = mod(atan2d(sind(lon_grid - radar.lon(i_r)) .* cosd(lat_grid), ...
        cosd(radar.lat(i_r)) .* sind(lat_grid) - sind(radar.lat(i_r)) .* cosd(lat_grid) .* cosd(lon_grid - radar.lon(i_r))), 360);
end
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 3, ceil(n_night / 3), TileSpacing="tight", Padding="tight");
title(t, "Lowest cut reflectivity, 3 h after sunset (raw data; each pixel from its nearest radar)", FontSize=14)
for i_n = 1:n_night
    mosaic = nan(size(lon_grid));
    nearest_km = inf(size(lon_grid));
    for i_r = 1:n_radar
        sweep = lowest{(i_r - 1) * n_night + i_n};
        if isempty(sweep)
            continue
        end
        d = distance_km(:, :, i_r);
        use = d < nearest_km & d <= 150;
        range_id = interp1(double(sweep.range) / 1000, 1:numel(sweep.range), d(use), "nearest", NaN);
        a = azimuth_deg(:, :, i_r);
        [~, az_id] = min(abs(mod(a(use) - double(sweep.az(:))' + 180, 360) - 180), [], 2);
        value = nan(nnz(use), 1);
        ok = ~isnan(range_id);
        value(ok) = sweep.REF(sub2ind(size(sweep.REF), range_id(ok), az_id(ok)));
        mosaic(use) = value;
        nearest_km(use) = d(use);
    end
    ax = nexttile(t); hold(ax, "on"); box(ax, "on")
    set(ax, Color=[0.85 0.92 0.97], FontSize=11)
    imagesc(ax, lon, lat, mosaic, AlphaData=~isnan(mosaic))
    for s = states'
        plot(ax, s.Lon, s.Lat, Color=[0.4 0.4 0.4], LineWidth=0.5)
    end
    plot(ax, radar.lon, radar.lat, "k.", MarkerSize=10)
    colormap(ax, parula); clim(ax, [-10 25])
    axis(ax, "xy"); xlim(ax, c.figure_bounds(1:2)); ylim(ax, c.figure_bounds(3:4))
    title(ax, sprintf("%s %s UTC", string(selected.date(i_n), "yyyy-MM-dd"), string(target_time(i_n), "HH:mm")), ...
        FontWeight="bold", FontSize=12)
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Label.String = "Reflectivity, lowest cut (dBZ)";
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "night_maps.png"), Resolution=150)
close(fig)
