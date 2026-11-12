%% Artefacts A: small test cases cut from the us_east night (dev)
% Cases of the 21-22 October 2026 night, 04:30-07:30 UTC every 30 min
% (7 frames), one per artefact seen in the us_east video:
% - ar_midwest: the sandbox radars (28, incl. KDMX: wedge NE of the radar,
%   screenshot 5);
% - ar_appalachians: radars of the southern Appalachians (column deficit
%   over the ridges and streaks downwind, screenshot 2);
% - ar_plains: radars of the High Plains up to the -104 deg domain edge
%   (KLNX, KGLD, KUDX discs, screenshot 1).
% Raw files and clutter masks are those of us_east (ar_test_case). Each
% case is built with up to four cleanings, each with its own volume cache:
% 1. us_east (round 0: elliptical shape correction, beam-blockage maps read
%    with the old, mirrored lookup) - cached, only the observations are built;
% 2. ar_midwest_noshape: no shape correction;
% 3. ar_appalachians_bb, ar_plains_bb: beam-blockage maps read at 90 - az
%    (dev override, dev/+bmm4d/+radar/+filter/beam_blockage.m);
% 4. ar_*_fix: both cleaning fixes, the blockage as 3 and the symmetric
%    aspect correction (ar_clean_volume) instead of the shape template;
%    cache us_east_fix, shared with the full night (AR_E_us_east.m);
% 5. ar_midwest_none, ar_midwest_shared: the ar_midwest_fix volumes without
%    aspect correction, and with the shared aspect function of
%    studies/aspect (no re-cleaning; AR_I_aspect.m).
% Cleanings 2-4 run as batch jobs (4 workers; about 30 min per 200
% volumes). The client's own pool (build_observations, 6 workers) must fit
% beside a job: 12 workers on this machine.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
dev_folder = fullfile(project_root, "studies", "artifacts", "dev");
volumes_root = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
boxes = struct("midwest", [-100 -85 36 46], "appalachians", [-89 -78 33.5 40], "plains", [-104 -94 37.5 47.5]);
no_shape = {"read_options", {"shape", false}};

%% 1. Cases with the us_east cleaning ---------------------------------------
for name = ["midwest" "appalachians" "plains"]
    c = ar_test_case("ar_" + name, boxes.(name), times);
    if ~isfile(fullfile(c.data_dir, "case.mat"))
        bmm4d.clean_volumes(c);          % all cached: only writes volume_files
        bmm4d.build_observations(c);
    end
end

%% 2-4. Other cleanings in the background -------------------------------------------
c_noshape = ar_test_case("ar_midwest_noshape", boxes.midwest, times, no_shape{:});
jobs = batch(@ar_clean_volumes, 2, {c_noshape, project_root}, Pool=4, AutoAddClientPath=true, ...
    CurrentFolder=project_root);
bb_cases = {ar_test_case("ar_appalachians_bb", boxes.appalachians, times)
    ar_test_case("ar_plains_bb", boxes.plains, times)};
jobs(2) = batch(@ar_clean_cases, 1, {bb_cases, project_root, "dev_folder", dev_folder, ...
    "cache_dir", fullfile(volumes_root, "9c2b4ef2_bbfix")}, Pool=4, AutoAddClientPath=true, CurrentFolder=project_root);
fix_cases = {ar_test_case("ar_midwest_fix", boxes.midwest, times, no_shape{:})
    ar_test_case("ar_appalachians_fix", boxes.appalachians, times, no_shape{:})
    ar_test_case("ar_plains_fix", boxes.plains, times, no_shape{:})};
jobs(3) = batch(@ar_clean_cases, 1, {fix_cases, project_root, "dev_folder", dev_folder, ...
    "cache_dir", fullfile(volumes_root, "us_east_fix"), "aspect", true}, Pool=4, AutoAddClientPath=true, CurrentFolder=project_root);

%% Observations, once the jobs have finished --------------------------------------------
wait(jobs(1))
bmm4d.build_observations(c_noshape);
wait(jobs(2))
cellfun(@bmm4d.build_observations, bb_cases, UniformOutput=false);
wait(jobs(3))
cellfun(@bmm4d.build_observations, fix_cases, UniformOutput=false);

%% 5. Shared aspect function (studies/aspect) on the ar_midwest_fix volumes ---------------
% ar_midwest_none: the ar_midwest_fix volumes (beam-blockage fix) with the
% per-volume symmetric fit undone (as_read_volume); ar_midwest_shared: the
% same divided by the shared aspect function of studies/aspect
% (as_apply_aspect, its parameters for this night). No re-cleaning.
addpath(fullfile(project_root, "studies", "aspect", "functions"))
aspect = struct("b2", 0, "b2_per_km", 0.092, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, "b1", -0.21, "b1_ramp_m", 1000);
source = load(fullfile(fix_cases{1}.data_dir, "manifest.mat"), "volume_files");
for name = ["none" "shared"]
    c_aspect = ar_test_case("ar_midwest_" + name, boxes.midwest, times, no_shape{:});
    cache_dir = fullfile(volumes_root, "aspect_" + name);
    if ~isfolder(cache_dir)
        mkdir(cache_dir)
    end
    volume_files = strings(size(source.volume_files));
    for k = find(source.volume_files ~= "")'
        [~, file_name] = fileparts(source.volume_files(k));
        volume_files(k) = fullfile(cache_dir, file_name + ".mat");
        volume = rmfield(as_read_volume(source.volume_files(k)), "aspect");
        if name == "shared"
            volume = as_apply_aspect(volume, aspect);
        end
        save(volume_files(k), "volume")
    end
    save(fullfile(c_aspect.data_dir, "manifest.mat"), "volume_files", "-append")
    bmm4d.build_observations(c_aspect);
end
