function S = nr_apply_drop(S, drop, values)
%NR_APPLY_DROP Remove observation entries and rebuild the prior and weights (dev).
%   S = NR_APPLY_DROP(S, DROP) takes a case loaded by ar_load and removes
%   the entries DROP [n_row n_time] (nr_drop_rows) from the observations
%   (values and sd set to 0, i.e. missing; a no-echo observation is eps,
%   not 0). The density prior with its radar profiles and the production
%   weights are rebuilt from the remaining entries, as ar_load builds them,
%   so the dropped data do not reach the solve through the prior.
%   S = NR_APPLY_DROP(S, DROP, VALUES) first replaces the observed values
%   (e.g. synthetic observations F x truth on the same entries).
%   S.dropped keeps DROP for scoring.

inv = S.c.inversion;
if nargin > 2
    S.observations.values = values;
end
S.observations.values(drop) = 0;
if isfield(S.observations, "sd")
    S.observations.sd(drop) = 0;
end
[S.prior, prior_info] = bmm4d.prior.build_density_prior_fast(S.g, S.radar_scans, S.observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", S.active_id);
S.radar_profiles = prior_info.radar_profiles;
S.row_weights = [];
if inv.observation_weights && isfield(S.observations, "sd")
    S.row_weights = bmm4d.solve.observation_weights(S.observations, ...
        S.observations.F(:, S.active_id) * double(S.prior), inv.observation_error_floor);
end
S.dropped = drop;
end
