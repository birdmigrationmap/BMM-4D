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
%   Suffixes change one setting: _rho, _dr, _mistnet or _none the rain
%   removal (MistNet by default), _coarse a 4 km x 200 m solve grid,
%   _nearNN drops gates within NN km, _zeroNN keeps no-echo gates only
%   within NN km. E.g. "keywest_rho", "us_east_near15".

name_parts = split(string(case_name), "_");
base_name = name_parts(1);
if base_name == "us"
    base_name = "us_east";
end
radar_file = bmm4d.util.data_path("reference", "radar.csv");

% Case settings ---------------------------------------------------------------
switch base_name
    case "keywest"
        % Paper Key West cleaning options and published weights
        settings = {"radars", "KBYX", ...
            "times", datetime("2018-04-23 01:00:00"):minutes(15):datetime("2018-04-23 07:00:00"), ...
            "closest_accepted_min", 15, ...
            "grid_bounds", [-83.2, -80.2, 22.7, 25.6], ...
            "figure_bounds", [-83.2, -80.0, 22.8, 25.7], ...
            "water_override", [-82.05 -81.98 -81.93 -81.97 -82.04; 24.02 24.07 24.02 23.97 23.98], ...
            "read_options", {"removeOutliers", 2}, ...
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
        % Same cleaning as us_east (shared volume cache), no data cuts
        radar_bounds = [-100 -85 36 46];
        settings = {"radars", bmm4d.radar.select_radars_in_bbox(radar_file, radar_bounds).name', ...
            "times", datetime("2026-10-22 00:30"):minutes(30):datetime("2026-10-22 03:00"), ...
            "grid_bounds", radar_bounds, ...
            "filter", struct("data_res", [2000 2000]), ...
            "solve_aggregate", [2 2 1], ...
            "domain", struct("west_lon", -180, "offshore_km", 50), ...
            "near_range_m", 0, ...
            "no_echo_max_range_m", Inf};
    otherwise
        error("case_settings:UnknownCase", "Unknown case %s.", case_name)
end

% Suffixes ----------------------------------------------------------------------
rain_part = name_parts(ismember(name_parts, ["rho" "dr" "mistnet" "none"]));
rain_method = "mistnet";
if ~isempty(rain_part)
    rain_method = rain_part(1);
end
read_id = find(strcmp(cellfun(@string, settings(1:2:end), UniformOutput=false), "read_options"));
if isempty(read_id)
    settings = [settings {"read_options", {"removeRain", rain_method}}];
else
    settings{2 * read_id} = [settings{2 * read_id} {"removeRain", rain_method}];
end
if any(name_parts == "coarse")
    settings = [settings {"solve_aggregate", [2 2 2]}];
end
near_part = regexp(case_name, "near(\d+)", "tokens", "once");
if ~isempty(near_part)
    settings = [settings {"near_range_m", str2double(near_part{1}) * 1000}];
end
zero_part = regexp(case_name, "zero(\d+)", "tokens", "once");
if ~isempty(zero_part)
    settings = [settings {"no_echo_max_range_m", str2double(zero_part{1}) * 1000}];
end
c = bmm4d.config(case_name, settings{:});
end
