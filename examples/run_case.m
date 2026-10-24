%% Reconstruct one case: volumes, cleaning, observations, inversion, maps
% Each stage saves its products in c.data_dir (manifest, static clutter,
% case, inversion) and logs its time and memory in stage_log.csv, so a
% stage can be rerun alone once the earlier ones are done. Cleaned volumes
% are cached and shared by every case with the same cleaning settings.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))

%% Case settings
% "keywest" (fast test), "midwest" (17 radars), "us_east" (one night,
% ~110 radars) or "sandbox"; see case_settings for the suffixes. A
% case_name already set in the workspace (batch runs) is kept.
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = case_settings(case_name);
% A new case is a call to bmm4d.config, e.g.
% c = bmm4d.config("kdox_test", radars="KDOX", grid_bounds=[-77 -74 37 40], ...
%     times=datetime("2021-09-20 01:00"):minutes(30):datetime("2021-09-20 10:00"));

%% Raw volumes
raw_files = bmm4d.get_volumes(c);

%% Cleaning
volume_files = bmm4d.clean_volumes(c);
clutter_masks = bmm4d.static_clutter(c);

%% Observation system
obs_case = bmm4d.build_observations(c);

%% Inversion
[active_field, active_id, inversion_info] = bmm4d.invert(c);

%% Maps
bmm4d.export_maps(c)
