function c = nr_night_case(name, source_dir, radar_bounds, max_frames)
%NR_NIGHT_CASE Regional test case cut from a cleaned multi-radar night (dev).
%   C = NR_NIGHT_CASE(NAME, SOURCE_DIR, RADAR_BOUNDS) returns the settings
%   of a case restricted to the radars inside RADAR_BOUNDS ([lon_min
%   lon_max lat_min lat_max]) and to all frames of the night whose case
%   folder is SOURCE_DIR (e.g. products/prod/night_us_YYYYMMDD of the
%   shared night set, whose volumes are already cleaned). As ar_test_case,
%   but from any night: the grid is the radar extent plus 160 km, the
%   cleaning, thinning, solve grid and domain are those of us_east.
%   <local_data_dir>/prod/<NAME> gets the source's manifest rows of these
%   radars (raw and cleaned volume files, as this machine sees them:
%   nr_shared_path) and its static clutter masks (the us_east masks when
%   the source has none), so bmm4d.build_observations(C) builds the case
%   with no download and no cleaning (it opens a 6-worker pool).
%   MAX_FRAMES (default Inf) keeps the first frames only.

if nargin < 4
    max_frames = inf;
end
us_east = case_settings("us_east");
source = load(fullfile(source_dir, "manifest.mat"), "radars", "raw_files", "volume_files", "times");
frames = 1:min(numel(source.times), max_frames);
source.times = source.times(frames);
source.raw_files = source.raw_files(:, frames);
source.volume_files = source.volume_files(:, frames);
radar = bmm4d.radar.select_radars_in_bbox(bmm4d.util.data_path("reference", "radar.csv"), radar_bounds);
radar = radar(ismember(radar.name, source.radars), :);
pad_lat = 160 / 111.32;
pad_lon = 160 / (111.32 * cosd(max(radar.lat)));
grid_bounds = [min(radar.lon) - pad_lon, max(radar.lon) + pad_lon, min(radar.lat) - pad_lat, max(radar.lat) + pad_lat];
c = bmm4d.config(name, "radars", radar.name', "times", source.times, "grid_bounds", grid_bounds, ...
    "read_options", {"shape", false}, "filter", us_east.filter, "solve_aggregate", us_east.solve_aggregate, ...
    "domain", us_east.domain);

% Manifest and static clutter masks of these radars
if ~isfolder(c.data_dir)
    mkdir(c.data_dir)
end
% bmm4d.config can drop frames (sun limit): keep the source columns of c.times
[~, radar_id] = ismember(c.radars, source.radars);
[found_time, time_id] = ismember(c.times, source.times);
assert(all(found_time), "nr_night_case:Times", "Case times missing from the source manifest.")
radars = c.radars;
times = c.times;
raw_files = nr_shared_path(source.raw_files(radar_id, time_id));
volume_files = nr_shared_path(source.volume_files(radar_id, time_id));
save(fullfile(c.data_dir, "manifest.mat"), "radars", "raw_files", "volume_files", "times")
clutter_file = fullfile(source_dir, "static_clutter.mat");
if ~isfile(clutter_file)
    clutter_file = fullfile(us_east.data_dir, "static_clutter.mat");
end
clutter = load(clutter_file);
[found, clutter_id] = ismember(c.radars, clutter.clutter_radars);
clutter_masks = cell(numel(c.radars), 1);
clutter_masks(found) = clutter.clutter_masks(clutter_id(found));
clutter_radars = c.radars;
settings = clutter.settings;
save(fullfile(c.data_dir, "static_clutter.mat"), "clutter_masks", "clutter_radars", "settings")
end
