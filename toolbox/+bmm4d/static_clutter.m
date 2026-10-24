function clutter_masks = static_clutter(c)
%STATIC_CLUTTER Static clutter mask of each radar of a case.
%   CLUTTER_MASKS = BMM4D.STATIC_CLUTTER(C) returns, per radar, the gates
%   that stand out from their surroundings at the same place in most
%   volumes of the night (wind farms, towers, persistent ground clutter;
%   bmm4d.radar.static_clutter_mask with C.static_clutter). Saves
%   static_clutter.mat in C.data_dir; build_observations removes these
%   gates from the low sweeps when the file exists.

load(fullfile(c.data_dir, "manifest.mat"), "raw_files")

%% Masks per radar, in parallel
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
bmm4d.util.stage_log(c, "static_clutter", stage_tic, sprintf("median %.2f%%, max %.2f%% of low-sweep gates flagged", ...
    100 * median(flagged), 100 * max(flagged)))
end
