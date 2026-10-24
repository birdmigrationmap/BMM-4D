%% Production test 2b: static clutter masks
% One mask per radar of the gates that stand out from their surroundings
% at the same place in most volumes of the night (wind farms, towers,
% persistent ground clutter). P3 removes them from the low sweeps.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "pipeline", "functions"))
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
load(fullfile(c.data_dir, "manifest.mat"), "raw_files")

%% Masks per radar, in parallel -------------------------------------------
stage_tic = tic;
pool = gcp("nocreate");
if isempty(pool) || pool.NumWorkers ~= 6
    delete(pool);
    parpool("Processes", 6);
end
settings = c.static_clutter;
n_radar = numel(c.radars);
clutter_masks = cell(n_radar, 1);
parfor i_r = 1:n_radar
    files = raw_files(i_r, raw_files(i_r, :) ~= "");
    clutter_masks{i_r} = bmm4d.radar.static_clutter_mask(files, settings);
end
delete(gcp("nocreate"))
clutter_radars = c.radars;
save(fullfile(c.data_dir, "static_clutter.mat"), "clutter_masks", "clutter_radars", "settings")
has_mask = ~cellfun(@isempty, clutter_masks);
flagged = cellfun(@(m) mean(m.flag(:)), clutter_masks(has_mask));
stage_log(c, "P2b_static_clutter", stage_tic, sprintf("median %.2f%%, max %.2f%% of low-sweep gates flagged", ...
    100 * median(flagged), 100 * max(flagged)))

