%% Load paths and reference data
project_root = setup_paths();

radar = readtable(data_path("reference", "radar.csv"), TextType="string");
load(data_path("reference", "DEM", "dem.mat"), "dem", "dem_lat", "dem_lon")

dem = double(dem);
dem(dem == -9999) = -2;
dem_f = griddedInterpolant({dem_lat, dem_lon}, dem, "nearest");

%% Define the polar terrain-horizon grid
bb_range = 1000:250:500000;
bb_az = 0:1:365;
bb = nan(numel(bb_range), numel(bb_az), height(radar), "single");

[AZ, R] = meshgrid(bb_az, bb_range);
[x, y] = pol2cart(deg2rad(AZ), R);

%% Build a terrain-horizon elevation map for each radar
for i_r = 1:height(radar)
    lat_m = deg2km(distance(radar.lat(i_r), radar.lon(i_r), radar.lat(i_r) + 1, radar.lon(i_r))) * 1000;
    lon_m = deg2km(distance(radar.lat(i_r), radar.lon(i_r), radar.lat(i_r), radar.lon(i_r) + 1)) * 1000;

    lat_r = y / lat_m + radar.lat(i_r);
    lon_r = x / lon_m + radar.lon(i_r);
    dem_r = max(0, dem_f(lat_r, lon_r));

    [~, terrain_elev, ~] = slant2range(R(:), radar.height(i_r), dem_r(:));
    terrain_elev = reshape(terrain_elev, size(R));

    % Once terrain rises above a sight line, it blocks all farther gates on
    % the same azimuth until an equal or higher elevation is used.
    bb(:,:,i_r) = single(cummax(terrain_elev, 1));
end

%% Save reference products
output_dir = data_path("reference", "beam_blockage");
bb_radar = radar.name;
save(fullfile(output_dir, "beam_blockage.mat"), "bb", "bb_range", "bb_az", "bb_radar", "-v7.3")

for i_r = 1:height(radar)
    bbr = bb(:,:,i_r);
    save(fullfile(output_dir, "beam_blockage_" + radar.name(i_r) + ".mat"), "bbr", "bb_range", "bb_az", "-v7.3")
end
