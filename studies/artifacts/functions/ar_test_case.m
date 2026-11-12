function c = ar_test_case(name, radar_bounds, times, varargin)
%AR_TEST_CASE Small test case cut from the us_east night (dev).
%   C = AR_TEST_CASE(NAME, RADAR_BOUNDS, TIMES, ...) returns the settings of
%   the us_east night (21-22 October 2026) restricted to the radars inside
%   RADAR_BOUNDS ([lon_min lon_max lat_min lat_max]) and to the frames
%   TIMES. The grid is the radar extent plus 160 km, the cleaning, thinning
%   and solve grid are those of us_east (shared volume cache), the solve
%   domain is the us_east domain (land east of -104 deg, 50 km offshore).
%   Other name-value pairs go to bmm4d.config after these (e.g.
%   "read_options", {"shape", false}).
%   The case folder prod/<NAME> gets the us_east manifest and static
%   clutter masks of these radars and times, so get_volumes and
%   static_clutter are not run: no download, the same clutter masks as
%   us_east. Then clean_volumes and build_observations build the case.

us_east = case_settings("us_east");
radar = bmm4d.radar.select_radars_in_bbox(bmm4d.util.data_path("reference", "radar.csv"), radar_bounds);
pad_lat = 160 / 111.32;
pad_lon = 160 / (111.32 * cosd(max(radar.lat)));
grid_bounds = [min(radar.lon) - pad_lon, max(radar.lon) + pad_lon, min(radar.lat) - pad_lat, max(radar.lat) + pad_lat];
c = bmm4d.config(name, "radars", radar.name', "times", times, "grid_bounds", grid_bounds, ...
    "read_options", us_east.read_options, "filter", us_east.filter, "solve_aggregate", us_east.solve_aggregate, ...
    "domain", us_east.domain, varargin{:});

% Manifest and static clutter masks: the us_east ones, for these radars and times
if ~isfolder(c.data_dir)
    mkdir(c.data_dir)
end
manifest_file = fullfile(c.data_dir, "manifest.mat");
if ~isfile(manifest_file)
    source = load(fullfile(us_east.data_dir, "manifest.mat"), "radars", "raw_files", "times");
    [~, radar_id] = ismember(c.radars, source.radars);
    [~, time_id] = ismember(c.times, source.times);
    assert(all(radar_id > 0) && all(time_id > 0), "ar_test_case:NotInUsEast", ...
        "Every radar and time must be in the us_east manifest.")
    radars = c.radars;
    times = c.times;
    raw_files = source.raw_files(radar_id, time_id);
    save(manifest_file, "radars", "raw_files", "times")
    clutter = load(fullfile(us_east.data_dir, "static_clutter.mat"));
    [~, radar_id] = ismember(c.radars, clutter.clutter_radars);
    clutter_masks = clutter.clutter_masks(radar_id);
    clutter_radars = c.radars;
    settings = clutter.settings;
    save(fullfile(c.data_dir, "static_clutter.mat"), "clutter_masks", "clutter_radars", "settings")
end
end
