function prior_next = ar_smooth_prior(S, prior, field, sigma_m)
%AR_SMOOTH_PRIOR Prior that carries the data's local level past their edge (dev).
%   PRIOR_NEXT = AR_SMOOTH_PRIOR(S, PRIOR, FIELD, SIGMA_M) adds to PRIOR
%   the smoothed deviation of a first solve FIELD from it, per layer and
%   time: N(d) = G * (s d) / max(G * s, 1/2), with d = FIELD - PRIOR, s = 1
%   in the columns seen by an observation of that time (any layer with
%   beam weight: the unseen low layers of a seen column are part of the
%   solved level too) and 0 elsewhere, and G a Gaussian of SIGMA_M. Where
%   the data reach, N(d) is the local mean deviation; past the edge of the
%   data it keeps that level and fades over about SIGMA_M instead of
%   falling back to the radar profiles at once. Both arrays are
%   [n_active n_time].

g = S.g;
n_time = size(field, 2);
sigma_cells = sigma_m / g.dx;
observed = S.observations.values ~= 0;
prior_next = prior;
for i_t = 1:n_time
    seen_cell = zeros(g.sz);
    seen_cell(S.active_id) = full(any(S.observations.F(observed(:, i_t), S.active_id), 1))';
    support = double(S.active_mask_3d & any(seen_cell, 3));
    deviation = zeros(g.sz);
    deviation(S.active_id) = double(field(:, i_t)) - double(prior(:, i_t));
    deviation = deviation .* support;
    smoothed = zeros(g.sz);
    for i_z = 1:g.sz(3)
        if ~any(support(:, :, i_z), "all")
            continue
        end
        weight = imgaussfilt(support(:, :, i_z), sigma_cells, Padding=0);
        smoothed(:, :, i_z) = imgaussfilt(deviation(:, :, i_z), sigma_cells, Padding=0) ./ max(weight, 0.5);
    end
    prior_next(:, i_t) = max(double(prior(:, i_t)) + smoothed(S.active_id), 0);
end
end
