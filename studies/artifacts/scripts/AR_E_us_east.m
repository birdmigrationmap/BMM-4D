%% Artefacts E: the us_east night with the fixes (dev)
% The fixes chosen on the test cases (AR_B-AR_D), applied to the whole
% us_east night (26 frames, 23:30-12:00 UTC), against round 0
% (studies/sandbox/SB_us_east_video.m):
% 1. cleaning: beam-blockage maps read at compass azimuth (dev override,
%    dev/+bmm4d/+radar/+filter/beam_blockage.m) and the symmetric aspect
%    correction instead of the elliptical shape template (ar_clean_volume);
% 2. inversion: advection links and prior at height above the ground
%    (ar_steps, ar_profile_prior "agl");
% 3. inversion: correlated prior, 100 km (prior_correlation_m).
% Section 1 cleans in a background batch job, cache us_east_fix (shared with
% the ar_*_fix test cases of AR_A_cases.m, whose volumes it reuses): about
% 2 h with 8 workers. The client cannot open its own pool meanwhile (12
% workers), so section 2 waits for the job.
% Outputs: prod/us_east_fix/case.mat, inversion_fix.mat (ar_invert_fix);
% video and spots by AR_F_video.m.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
dev_folder = fullfile(project_root, "studies", "artifacts", "dev");
cache_dir = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes", "us_east_fix");
us_east = case_settings("us_east");
c = ar_test_case("us_east_fix", [-104 -66 24 50], us_east.times, "read_options", {"shape", false});

%% 1. Clean in the background -------------------------------------------------------
clean_job = batch(@ar_clean_cases, 1, {{c}, project_root, "dev_folder", dev_folder, "cache_dir", cache_dir, ...
    "aspect", true}, Pool=8, AutoAddClientPath=true, CurrentFolder=project_root);

%% 2. Observations of the night, once the job has finished --------------------------
wait(clean_job)
bmm4d.build_observations(c);

%% 3. Inversion with the fixes, in the background ---------------------------------------
% About 100 min (352 PCG iterations with the GP prior) on 8 threads; the
% worker holds up to ~36 GB (most of it compressed or swapped on this
% machine), so nothing else large should run beside it.
invert_job = batch(@ar_invert_fix, 1, {c, project_root}, Pool=0, AutoAddClientPath=true, CurrentFolder=project_root);
wait(invert_job)
solve_info = fetchOutputs(invert_job);
disp(solve_info{1})

%% 4. Shared aspect function instead of the per-volume fit (AR_I_aspect.m) ---------------
% The cleaned volumes of section 1 with the per-volume symmetric fit undone
% and divided by the shared aspect function of studies/aspect (no
% re-cleaning, ~15 min), then the observations and the same inversion
% (prod/us_east_shared/inversion_fix.mat).
addpath(fullfile(project_root, "studies", "aspect", "functions"))
aspect = struct("b2", 0, "b2_per_km", 0.092, "b2_per_ms", 0, "h0_m", 0, "v0_ms", 0, "b1", -0.21, "b1_ramp_m", 1000);
c_shared = ar_test_case("us_east_shared", [-104 -66 24 50], us_east.times, "read_options", {"shape", false});
shared_dir = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes", "us_east_shared");
if ~isfolder(shared_dir)
    mkdir(shared_dir)
end
source = load(fullfile(c.data_dir, "manifest.mat"), "volume_files");
volume_files = strings(size(source.volume_files));
for k = find(source.volume_files ~= "")'
    [~, file_name] = fileparts(source.volume_files(k));
    volume_files(k) = fullfile(shared_dir, file_name + ".mat");
    if isfile(volume_files(k))
        continue
    end
    volume = as_apply_aspect(rmfield(as_read_volume(source.volume_files(k)), "aspect"), aspect);
    save(volume_files(k), "volume")
end
save(fullfile(c_shared.data_dir, "manifest.mat"), "volume_files", "-append")
bmm4d.build_observations(c_shared);
invert_job = batch(@ar_invert_fix, 1, {c_shared, project_root}, Pool=0, AutoAddClientPath=true, CurrentFolder=project_root);
