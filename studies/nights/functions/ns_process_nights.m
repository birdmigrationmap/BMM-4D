function run_log = ns_process_nights(names, max_frames)
%NS_PROCESS_NIGHTS Download and clean nights of the night set, for a cluster job (dev).
%   RUN_LOG = NS_PROCESS_NIGHTS(NAMES, MAX_FRAMES) processes the rows NAMES
%   of studies/nights/figures/NS_A_nights.csv in turn (one MATLAB job on
%   merops-02, submitted with hpc/submit.sh):
%   - full nights ("night_us_*"; purpose core, appalachians, plains,
%     reference): the us_east radars every 30 min, raw volumes
%     (bmm4d.get_volumes, downloads into nexrad_raw_dir), cleaning as
%     production but without the shape correction (read option "shape"
%     false: one cache every study builds on; the shared aspect function
%     is applied later with as_apply_aspect, without cleaning again), and
%     the static clutter masks;
%   - rain days ("rain_us_*"): raw volumes only, daytime kept (sun limit
%     90 deg). The radar-pattern analysis reads the reflectivity itself;
%     the cleaning would remove the rain;
%   - the paper cases ("keywest", "midwest"): raw volumes only, cleaned
%     when the method of the v2 rerun is chosen.
%   MAX_FRAMES (default Inf) keeps the first frames only (pilot). The
%   cleaned volumes go to the shared cache of these settings, so a pilot's
%   volumes are reused by the full run. Raw files are kept.
%   One row per night is appended to <local_data_dir>/nights/run_log.csv.

if nargin < 2
    max_frames = inf;
end
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
nights = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string");
nights.first_frame = datetime(nights.first_frame);
nights.last_frame = datetime(nights.last_frame);
us = case_settings("us_east");
log_dir = fullfile(bmm4d.util.local_data_dir(), "nights");
if ~isfolder(log_dir)
    mkdir(log_dir)
end
run_log = table();
for name = string(names(:))'
    row = nights(nights.name == name, :);
    assert(height(row) == 1, "ns_process_nights:Name", "%s is not in NS_A_nights.csv.", name)
    times = row.first_frame:minutes(row.step_min):row.last_frame;
    times = times(1:min(numel(times), max_frames));
    full_night = startsWith(name, "night_us_");
    if full_night
        c = bmm4d.config(name, "radars", us.radars, "times", times, "grid_bounds", us.grid_bounds, ...
            "filter", us.filter, "solve_aggregate", us.solve_aggregate, "domain", us.domain, ...
            "read_options", {"shape", false});
    elseif startsWith(name, "rain_us_")
        c = bmm4d.config(name, "radars", us.radars, "times", times, "grid_bounds", us.grid_bounds, ...
            "filter", us.filter, "solve_aggregate", us.solve_aggregate, "domain", us.domain, ...
            "read_options", {"shape", false}, "sun_limit_deg", 90);
    else
        c = case_settings(name);
    end

    % Raw volumes, then the cleaning of the full nights
    stage_tic = tic;
    raw_files = bmm4d.get_volumes(c);
    download_s = toc(stage_tic);
    [clean_s, n_cleaned] = deal(0);
    if full_night
        stage_tic = tic;
        volume_files = bmm4d.clean_volumes(c);
        bmm4d.static_clutter(c);
        clean_s = toc(stage_tic);
        n_cleaned = nnz(volume_files ~= "");
    end
    entry = table(name, string(datetime("now")), numel(c.times), nnz(raw_files ~= ""), n_cleaned, download_s, ...
        clean_s, string(c.data_dir), VariableNames=["name" "finished" "n_frames" "n_raw" "n_cleaned" "download_s" ...
        "clean_s" "data_dir"]);
    run_log = [run_log; entry]; %#ok<AGROW>
    writetable(entry, fullfile(log_dir, "run_log.csv"), WriteMode="append")
    fprintf("[nights] %s: %d frames, %d raw, %d cleaned, download %.0f s, cleaning %.0f s\n", name, ...
        numel(c.times), entry.n_raw, n_cleaned, download_s, clean_s)
end
end
