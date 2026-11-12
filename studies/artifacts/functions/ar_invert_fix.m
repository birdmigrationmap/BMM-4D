function solve_info = ar_invert_fix(c, project_root)
%AR_INVERT_FIX Inversion of a built case with the fixes (dev, batch job).
%   SOLVE_INFO = AR_INVERT_FIX(C, PROJECT_ROOT) solves case C with the
%   production settings plus advection links and prior at height above the
%   ground (ar_steps, ar_profile_prior "agl") and the 100 km correlated
%   prior, and saves inversion_fix.mat (active_field, active_id,
%   solve_info) in C.data_dir. Made to run as a MATLAB batch job, so the
%   client stays free during the 1-2 h solve of the us_east night.

% Batch workers start with one computational thread: give the sparse
% products of PCG 8 of the 12 cores (the client keeps the others)
maxNumCompThreads(8);
setup_paths(project_root);
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
S = ar_load(c);
inv = c.inversion;
steps_agl = ar_steps(S.g, S.radar_scans, seconds(diff(c.times)), S.active_mask_3d, inv, "agl");
prior_agl = ar_profile_prior(S.g, S.radar_scans, S.radar_profiles, S.active_mask_3d, inv.prior_length_scale_m, "agl");
fixed = ar_solve(S, "fix", steps=steps_agl, prior=prior_agl, prior_correlation_m=100000);
active_field = fixed.field;
active_id = S.active_id;
solve_info = rmfield(fixed, ["field" "column"]);
save(fullfile(c.data_dir, "inversion_fix.mat"), "active_field", "active_id", "solve_info", "-v7.3")
end
