function export_dir = export_webmap(c, field, options)
%EXPORT_WEBMAP Export an inversion for the web map viewer (webmap/).
%   BMM4D.EXPORT_WEBMAP(C) exports the inversion of the case
%   (C.data_dir/inversion.mat, active_field and active_id) to
%   C.data_dir/webmap/inversion/ and prints the link that opens it in the
%   viewer (python3 webmap/serve.py, http://127.0.0.1:8642/).
%
%   BMM4D.EXPORT_WEBMAP(C, FIELD, active_id=..., run=..., note=...) exports
%   any field instead, e.g. a study variant or a test: FIELD is
%   [n_active n_time], density in birds km^-3 on the active cells of the
%   case in C.data_dir/case.mat, as bmm4d.invert returns it.
%
%   The export folder (~10 MB and seconds to a minute for a us_east night):
%   - col_XX.png, one per frame: the column density (birds km^-2) on a grid
%     regular in longitude and in Web Mercator y (so the map places it
%     without distortion), red = log-encoded density, green = data mask;
%   - wind_XX.png, one per interval: the birds' mean wind (column-density
%     weighted over the layers, frames XX and XX+1 together), red = u,
%     green = v, used to morph between frames and to move the particles;
%   - meta.json: times, bounds, encodings, colour scale, radars, and what
%     the run is (case, run, note, time, host, commit or cluster job).
%
%   Options (name-value):
%     run          folder name and label of the run ("inversion"); an
%                  existing export of the same run is overwritten
%     note         one line shown in the viewer: what this run tests
%     active_id    active cells of FIELD (default: rebuilt from the case masks)
%     case_data    struct with g, radar_scans, radar_selected (e.g. the
%                  output of bmm4d.build_observations), to skip reading case.mat
%     out_dir      export folder (default C.data_dir/webmap/<run>)
%     meta_extra   struct of fields added to meta.json (e.g. solve time)
%     column_floor birds km^-2 of the log encoding (5)
%     wind_limit   m/s range of the wind encoding (50)
%
%   Example:
%     [field, active_id] = bmm4d.invert(c);
%     bmm4d.export_webmap(c, field, active_id=active_id, run="prior_shared", note="shared prior, 2 rounds")

arguments
    c (1,1) struct
    field (:,:) {mustBeNumeric} = []
    options.run (1,1) string = "inversion"
    options.note (1,1) string = ""
    options.active_id (:,1) double = []
    options.case_data (1,1) struct = struct()
    options.out_dir (1,1) string = ""
    options.meta_extra (1,1) struct = struct()
    options.column_floor (1,1) double = 5
    options.wind_limit (1,1) double = 50
end

%% Inversion and case: field, grid, radars, active cells
active_id = options.active_id;
if isempty(field)
    saved = load(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id");
    field = saved.active_field;
    active_id = saved.active_id;
end
case_data = options.case_data;
if isempty(fieldnames(case_data))
    case_data = load(fullfile(c.data_dir, "case.mat"), "g", "radar_scans", "radar_selected");
end
g = case_data.g;
radar_scans = case_data.radar_scans;
if isempty(active_id)
    column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
    active_mask = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
    active_id = find(active_mask(:));
end
assert(numel(active_id) == size(field, 1), "bmm4d:export_webmap:Size", ...
    "The field has %d rows but the case has %d active cells.", size(field, 1), numel(active_id))
n_time = size(field, 2);
cell_air = zeros(g.sz, "single");
cell_air(active_id) = bmm4d.grid.air_fraction(g, active_id);
active_column = false(g.sz(1:2));
[ix, iy, ~] = ind2sub(g.sz, active_id);
active_column(sub2ind(g.sz(1:2), ix, iy)) = true;

%% Target grid: regular in longitude and in Web Mercator y
% Same longitudes as the grid; rows from north to south with the
% longitude step in Mercator y (square pixels on the map).
mercator_y = @(lat) log(tan(pi / 4 + deg2rad(lat) / 2));
lon_step = mean(diff(g.lon));
y_north = mercator_y(max(g.lat));
n_rows = floor((y_north - mercator_y(min(g.lat))) / deg2rad(lon_step)) + 1;
row_y = y_north - (0:n_rows - 1)' * deg2rad(lon_step);
row_lat = rad2deg(atan(sinh(row_y)));
% [nx ny] grid map -> [n_rows nx] image (north up): interpolation along
% latitude only, normalized by the interpolated data mask
to_image = @(map) interp1(g.lat(:), map', row_lat, "linear", 0);
to_masked_image = @(map, mask) to_image(map .* mask) ./ max(to_image(double(mask)), eps);
image_mask = to_image(double(active_column)) >= 0.5;

%% Columns
export_dir = options.out_dir;
if strlength(export_dir) == 0
    export_dir = fullfile(c.data_dir, "webmap", options.run);
end
if ~isfolder(export_dir)
    mkdir(export_dir)
end
old_frames = dir(fullfile(export_dir, "*.png"));     % a rerun may have fewer frames
if ~isempty(old_frames)
    delete(fullfile(export_dir, "*.png"))
end
columns = zeros(n_rows, numel(g.lon), n_time, "single");
for i_t = 1:n_time
    frame = zeros(g.sz, "single");
    frame(active_id) = max(field(:, i_t), 0);
    column = sum(frame .* cell_air, 3) * g.dz / 1000;
    columns(:, :, i_t) = to_masked_image(double(column), active_column);
end
columns(repmat(~image_mask, 1, 1, n_time)) = 0;
column_max = max(max(columns, [], "all"), eps);
column_floor = options.column_floor;
encode_column = @(v) uint8(1 + round(254 * log1p(v / column_floor) / log1p(column_max / column_floor)));
% Outside the data the red channel repeats the nearest data pixel, so the
% browser's bilinear filter does not darken the edge (the green mask hides it)
[~, nearest_data] = bwdist(image_mask);
for i_t = 1:n_time
    red = encode_column(columns(:, :, i_t));
    red(~image_mask) = red(nearest_data(~image_mask));
    imwrite(cat(3, red, uint8(255 * image_mask), zeros(size(red), "uint8")), ...
        fullfile(export_dir, sprintf("col_%02d.png", i_t - 1)))
end

%% Birds' mean wind per interval
% Radar velocity profiles of both frames, to the grid with the radar
% weights of the advection (velocity settings of the case), averaged over
% the layers with the birds of both frames
inv = c.inversion;
wind_limit = options.wind_limit;
n_radar = size(radar_scans, 1);
interpolation_weights = bmm4d.prior.build_radar_weights(g, radar_scans, interpolation=inv.velocity_interpolation, ...
    length_scale_m=inv.velocity_length_scale_m, idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power);
encode_wind = @(w) uint8(round(255 * (min(max(w, -wind_limit), wind_limit) + wind_limit) / (2 * wind_limit)));
for i_t = 1:n_time - 1
    profiles = zeros(g.sz(3), 2, n_radar);
    for i_r = 1:n_radar
        profiles(:, :, i_r) = 0.5 * (bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t}, g) + ...
            bmm4d.prior.radar_uv_on_grid(radar_scans{i_r, i_t + 1}, g));
    end
    birds = zeros(g.sz, "single");
    birds(active_id) = (max(field(:, i_t), 0) + max(field(:, i_t + 1), 0)) .* cell_air(active_id);
    birds_column = sum(birds, 3);
    wind = cell(1, 2);
    for i_c = 1:2
        profile_matrix = fillmissing(reshape(permute(profiles(:, i_c, :), [3 1 2]), n_radar, g.sz(3)), "constant", 0);
        layer_wind = reshape(interpolation_weights * profile_matrix, g.sz);
        % Where there are no birds, the plain mean over the layers
        mean_wind = sum(birds .* layer_wind, 3) ./ max(birds_column, eps);
        layer_mean = mean(layer_wind, 3);
        mean_wind(birds_column <= 0) = layer_mean(birds_column <= 0);
        wind{i_c} = to_image(double(mean_wind));
    end
    imwrite(cat(3, encode_wind(wind{1}), encode_wind(wind{2}), zeros(n_rows, numel(g.lon), "uint8")), ...
        fullfile(export_dir, sprintf("wind_%02d.png", i_t - 1)))
end

%% Meta: map placement, encodings, and what this run is
radar_selected = case_data.radar_selected;
% A cell of structs, so that one radar is still a JSON list
radars = num2cell(struct("name", cellstr(string(radar_selected.name)), "lon", num2cell(radar_selected.lon), ...
    "lat", num2cell(radar_selected.lat)));
[~, host] = system("hostname");
meta = struct( ...
    "format", "bmm4d-webmap-1", ...
    "case", string(c.name), "run", options.run, "note", options.note, ...
    "created", string(datetime("now", TimeZone="UTC"), "yyyy-MM-dd'T'HH:mm:ss'Z'"), ...
    "host", strtrim(string(host)), "commit", git_commit(), "job", cluster_job(), ...
    "times", {cellstr(string(c.times(1:n_time), "yyyy-MM-dd'T'HH:mm:ss'Z'"))}, ...
    "width", numel(g.lon), "height", n_rows, ...
    "lon_west", g.lon(1) - lon_step / 2, "lon_east", g.lon(end) + lon_step / 2, ...
    "lat_north", rad2deg(atan(sinh(y_north + deg2rad(lon_step) / 2))), ...
    "lat_south", rad2deg(atan(sinh(row_y(end) - deg2rad(lon_step) / 2))), ...
    "column_floor", column_floor, "column_max", double(column_max), ...
    "column_limit", double(quantile(columns(columns > 0), 0.995)), ...
    "wind_limit", wind_limit, "sun_limit_deg", c.sun_limit_deg, ...
    "colormap", round(255 * parula(256)), ...
    "radars", {radars});
for name = string(fieldnames(options.meta_extra))'
    meta.(name) = options.meta_extra.(name);
end
writelines(jsonencode(meta, PrettyPrint=true), fullfile(export_dir, "meta.json"))

%% Link for the viewer
% The viewer lists exports under its roots: "local" = this machine's
% local_data_dir, "cluster" = the cluster's local_data_dir (products/ on
% the Vogelzug share, cluster_products_dir on the laptop)
fprintf("Exported %d frames (%d x %d) to %s\n", n_time, numel(g.lon), n_rows, export_dir)
data_root = string(bmm4d.util.local_data_dir());
if startsWith(export_dir, data_root + filesep)
    root = "local";
    if isunix && ~ismac
        root = "cluster";
    end
    id = root + "/" + strrep(extractAfter(export_dir, strlength(data_root) + 1), filesep, "/");
    fprintf("Web map: http://127.0.0.1:8642/?data=%s\n", id)
else
    fprintf("Web map: outside local_data_dir; serve it with webmap/serve.py --root name=<folder above it>\n")
end
end

function commit = git_commit()
% Short commit of the working tree, "+" if it has uncommitted changes ("" on a cluster snapshot)
commit = "";
repository = string(fileparts(which("setup_paths")));
[status, out] = system("git -C """ + repository + """ rev-parse --short HEAD");
if status == 0
    commit = strtrim(string(out));
    [~, dirty] = system("git -C """ + repository + """ status --porcelain --untracked-files=no");
    if strlength(strtrim(string(dirty))) > 0
        commit = commit + "+";
    end
end
end

function job = cluster_job()
% Job folder of a cluster run (hpc/job.sh runs from <job folder>/code), "" elsewhere
job = "";
if strlength(string(getenv("SLURM_JOB_ID"))) > 0
    [~, job] = fileparts(fileparts(string(fileparts(which("setup_paths")))));
end
end
