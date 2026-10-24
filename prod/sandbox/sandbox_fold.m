function S_fold = sandbox_fold(S, excluded)
%SANDBOX_FOLD Sandbox case with the observations of some radars withheld.
%   EXCLUDED are radar indices (rows of S.radar_scans). Their observation
%   values become missing (sparse zeros), and they are removed from the
%   radar list used for the radar-profile prior and the regional shape
%   (row_radar_id remapped, withheld rows get id 0). The advection steps
%   (wind from all radars), the grid and the radar table used for the
%   radar-distance score are kept. The prior cache is new, since priors
%   depend on the radars used. S_fold.withheld_rows marks the rows to score.

n_radar = size(S.radar_scans, 1);
keep = true(n_radar, 1);
keep(excluded) = false;
new_id = zeros(n_radar, 1, "uint16");
new_id(keep) = 1:nnz(keep);
withheld_rows = ismember(S.observations.row_radar_id, excluded);

S_fold = S;
S_fold.observations.values(withheld_rows, :) = 0;
S_fold.observations.row_radar_id = new_id(S.observations.row_radar_id);
S_fold.radar_scans = S.radar_scans(keep, :);
S_fold.prior_cache = containers.Map();
S_fold.withheld_rows = withheld_rows;
S_fold.excluded = excluded;
end
