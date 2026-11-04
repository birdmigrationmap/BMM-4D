function c = case_settings(case_name)
%CASE_SETTINGS Settings of the named test cases (bmm4d.config).
%   C = CASE_SETTINGS(CASE_NAME) for
%   - "keywest": KBYX, 23 April 2018, 01:00-07:00 UTC every 15 min (fast
%     test, same radar, times and weights as the paper Key West case);
%   - "midwest": 17 radars, 16 May 2020, 01:30-12:00 UTC every 30 min (the
%     paper Midwest case; scaling test for a network);
%   - "us_east": radars east of -104 deg, one night (21-22 Oct 2026);
%   - "sandbox": Midwest part of the us_east night, 6 frames, for fast
%     option tests (studies/sandbox).
%   - "night_YYYYMMDD": the sandbox radars and grid on the night starting
%     on that local date, every 30 min from 1 h after sunset to 30 min
%     before sunrise at the domain centre (16-23 frames); the dry nights
%     of studies/sandbox (SB_night_screen).
%   Suffixes change one setting: _rho, _dr, _mistnet or _none the rain
%   removal (MistNet by default), _coarse a 4 km x 200 m solve grid,
%   _median the paper's gate cleaning (outlier filters and moving median)
%   instead of the gate estimator.
%   E.g. "keywest_rho", "sandbox_median".

name_parts = split(string(case_name), "_");
base_name = name_parts(1);
if base_name == "us"
    base_name = "us_east";
end
radar_file = bmm4d.util.data_path("reference", "radar.csv");

% Case settings ---------------------------------------------------------------
switch base_name
    case "keywest"
        % Paper Key West published weights
        settings = {"radars", "KBYX", ...
            "times", datetime("2018-04-23 01:00:00"):minutes(15):datetime("2018-04-23 07:00:00"), ...
            "closest_accepted_min", 15, ...
            "grid_bounds", [-83.2, -80.2, 22.7, 25.6], ...
            "figure_bounds", [-83.2, -80.0, 22.8, 25.7], ...
            "water_override", [-82.05 -81.98 -81.93 -81.97 -82.04; 24.02 24.07 24.02 23.97 23.98], ...
            "inversion", struct("velocity_idw_neighbors", 1, "prior", 0.003, "smoothness", 1, "advection", 3, ...
                "coeff_horver", -0.5, "shoreline_crossing_coeff", 0.3, "maxit", 300)};
    case "midwest"
        grid_bounds = [-97, -85, 38, 45];
        settings = {"radars", bmm4d.radar.select_radars_in_bbox(radar_file, grid_bounds).name', ...
            "times", datetime("2020-05-16") + (1.5:0.5:12) ./ 24, ...
            "grid_bounds", grid_bounds, ...
            "figure_bounds", [-97, -86, 38.5, 45]};
    case "us_east"
        night_date = datetime("2026-10-21");
        radar = bmm4d.radar.select_radars_in_bbox(radar_file, [-104 -66 24 50]);
        % Grid: radar extent plus 160 km
        pad_lat = 160 / 111.32;
        pad_lon = 160 / (111.32 * cosd(max(radar.lat)));
        settings = {"radars", radar.name', ...
            "times", (night_date + hours(23.5)):minutes(30):(night_date + hours(36)), ...
            "grid_bounds", [min(radar.lon) - pad_lon, max(radar.lon) + pad_lon, min(radar.lat) - pad_lat, max(radar.lat) + pad_lat], ...
            "filter", struct("data_res", [2000 2000]), ...
            "solve_aggregate", [2 2 1], ...   % 4 km x 100 m (4x fewer unknowns than 2 km x 100 m)
            "domain", struct("west_lon", -104, "offshore_km", 50)};   % US east of the High Plains, 50 km offshore
    case "sandbox"
        % Same cleaning as us_east (shared volume cache)
        radar_bounds = [-100 -85 36 46];
        settings = {"radars", bmm4d.radar.select_radars_in_bbox(radar_file, radar_bounds).name', ...
            "times", datetime("2026-10-22 00:30"):minutes(30):datetime("2026-10-22 03:00"), ...
            "grid_bounds", radar_bounds, ...
            "filter", struct("data_res", [2000 2000]), ...
            "solve_aggregate", [2 2 1], ...
            "domain", struct("west_lon", -180, "offshore_km", 50)};
    case "night"
        radar_bounds = [-100 -85 36 46];
        night_date = datetime(name_parts(2), InputFormat="yyyyMMdd");
        centre = [mean(radar_bounds(3:4)) mean(radar_bounds(1:2))];
        [sunset, sunrise] = bmm4d.util.getSunriseSunset(centre(1), centre(2), night_date + hours(24 - centre(2) / 15));
        sunset.TimeZone = "";
        sunrise.TimeZone = "";
        first = dateshift(sunset + hours(1), "start", "hour") + minutes(30) * ceil(minutes(sunset + hours(1) - ...
            dateshift(sunset + hours(1), "start", "hour")) / 30);
        settings = {"radars", bmm4d.radar.select_radars_in_bbox(radar_file, radar_bounds).name', ...
            "times", first:minutes(30):(sunrise - minutes(30)), ...
            "grid_bounds", radar_bounds, ...
            "filter", struct("data_res", [2000 2000]), ...
            "solve_aggregate", [2 2 1], ...
            "domain", struct("west_lon", -180, "offshore_km", 50)};
    otherwise
        error("case_settings:UnknownCase", "Unknown case %s.", case_name)
end

% Suffixes ----------------------------------------------------------------------
rain_part = name_parts(ismember(name_parts, ["rho" "dr" "mistnet" "none"]));
rain_method = "mistnet";
if ~isempty(rain_part)
    rain_method = rain_part(1);
end
read_id = find(string(settings(1:2:end)) == "read_options");
if isempty(read_id)
    settings = [settings {"read_options", {}}];
    read_id = numel(settings) / 2;
end
settings{2 * read_id} = [settings{2 * read_id} {"removeRain", rain_method}];
if any(name_parts == "coarse")
    settings = [settings {"solve_aggregate", [2 2 2]}];
end
if any(name_parts == "median")
    settings{2 * read_id} = [settings{2 * read_id} {"density", "median", "removeOutliers", 3}];
end
c = bmm4d.config(case_name, settings{:});
end
