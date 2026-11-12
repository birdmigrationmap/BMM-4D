function run_log = as_night_set_bins(names, out_root, n_worker)
%AS_NIGHT_SET_BINS Overlap bins of the night set, one file per frame (dev, cluster job).
%   RUN_LOG = AS_NIGHT_SET_BINS(NAMES, OUT_ROOT, N_WORKER) bins every frame
%   of the nights NAMES of studies/nights/figures/NS_A_nights.csv:
%   - full nights ("night_us_*"): the cleaned volumes of the night's
%     manifest (no shape correction) with as_night_bins;
%   - rain days ("rain_us_*"): the raw reflectivity with as_rain_bins.
%   Files go to OUT_ROOT/<name>/frame_XX.mat (XX = the frame's index in the
%   manifest); existing frames are skipped, and done.txt marks a finished
%   night. N_WORKER parallel workers over frames (0 = serial).
%   RUN_LOG: one row per night (frames, cells, seconds).

if nargin < 3
    n_worker = 6;
end
prod_root = fullfile(bmm4d.util.local_data_dir(), "prod");
radar_table = as_radar_sites(fullfile(prod_root, "night_us_20261021", "manifest.mat"));
if n_worker > 0 && isempty(gcp("nocreate"))
    parpool("Processes", n_worker);
end
run_log = table();
for name = string(names(:))'
    night_tic = tic;
    manifest = load(fullfile(prod_root, name, "manifest.mat"));
    out_dir = fullfile(out_root, name);
    if ~isfolder(out_dir)
        mkdir(out_dir)
    end
    times = manifest.times;
    n_frame = numel(times);
    n_cells = zeros(n_frame, 1);
    is_rain = startsWith(name, "rain_us_");
    volume_files = manifest.volume_files;
    raw_files = manifest.raw_files;
    parfor (i = 1:n_frame, max(n_worker, 0))
        out_file = fullfile(out_dir, sprintf("frame_%02d.mat", i));
        if isfile(out_file)
            continue
        end
        try
            if is_rain
                bins = as_rain_bins(raw_files(:, i), radar_table, times(i), i, out_file);
            else
                bins = as_night_bins(volume_files(:, i), times(i), i, out_file);
            end
            n_cells(i) = height(bins);
        catch err
            fprintf("%s frame %d failed: %s\n", name, i, err.message)
        end
    end
    times_file = fullfile(out_dir, "times.mat");
    save(times_file, "times")
    writelines(string(datetime("now")), fullfile(out_dir, "done.txt"))
    row = table(name, n_frame, sum(n_cells), toc(night_tic), VariableNames=["name" "n_frame" "n_cells_new" "seconds"]);
    disp(row)
    run_log = [run_log; row]; %#ok<AGROW>
end
end

function radar_table = as_radar_sites(manifest_file)
% Site of each radar (lat, lon, antenna height) from one cleaned volume per radar
manifest = load(manifest_file, "radars", "volume_files");
n_radar = numel(manifest.radars);
[lat, lon, height_m] = deal(nan(n_radar, 1));
for r = 1:n_radar
    files = manifest.volume_files(r, :);
    files = files(~ismissing(files) & files ~= "");
    for f = files
        if isfile(f)
            loaded = load(f, "volume");
            lat(r) = loaded.volume.lat;
            lon(r) = loaded.volume.lon;
            height_m(r) = loaded.volume.height;
            break
        end
    end
end
radar_table = table(string(manifest.radars(:)), lat, lon, height_m, VariableNames=["radar" "lat" "lon" "height_m"]);
end
