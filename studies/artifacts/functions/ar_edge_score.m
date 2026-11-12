function score = ar_edge_score(S, column)
%AR_EDGE_SCORE How visible the end of the radar data is in the column maps (dev).
%   SCORE = AR_EDGE_SCORE(S, COLUMN) for column maps [nx ny n_time]
%   (ar_columns). Per frame, a column is "seen" when an observation of that
%   frame has beam weight in it. On every seen column next to an unseen
%   one (the edge of the data), the mean log column of the seen columns
%   within 12 km is compared with that of the unseen columns within 12 km.
%   - jump: mean |difference| over the edge columns and frames, as a
%     ratio (exp; 1 = no step where the data end);
%   - inner: the same with the edge moved 40 km into the data (the
%     ordinary contrast of the field at that scale);
%   - edge: jump / inner in log units (1 = the end of the data does not
%     stand out).

g = S.g;
radius_cells = 12000 / g.dx;
[dx, dy] = ndgrid(-ceil(radius_cells):ceil(radius_cells));
disc = double(hypot(dx, dy) <= radius_cells);
observed = S.observations.values ~= 0;
active_column = any(S.active_mask_3d, 3);
jump = [];
inner = [];
for i_t = 1:size(column, 3)
    weight = zeros(g.sz);
    weight(S.active_id) = full(sum(S.observations.F(observed(:, i_t), S.active_id), 1))';
    seen = any(weight > 0, 3) & active_column;
    log_column = log(max(column(:, :, i_t), 1));
    log_column(~active_column) = 0;
    jump = [jump; side_difference(seen, active_column & ~seen, log_column, disc)]; %#ok<AGROW>
    core = imerode(seen, strel("disk", round(40000 / g.dx)));
    inner = [inner; side_difference(core, seen & ~core, log_column, disc)]; %#ok<AGROW>
end
score = struct("jump", exp(mean(jump, "omitnan")), "inner", exp(mean(inner, "omitnan")), ...
    "edge", mean(jump, "omitnan") / mean(inner, "omitnan"));
end

function difference = side_difference(inside, outside, log_column, disc)
% |mean log column of inside cells - of outside cells| within the disc,
% on inside cells that touch an outside cell
edge = inside & imdilate(outside, ones(3));
mean_inside = conv2(log_column .* inside, disc, "same") ./ conv2(double(inside), disc, "same");
mean_outside = conv2(log_column .* outside, disc, "same") ./ conv2(double(outside), disc, "same");
difference = abs(mean_inside(edge) - mean_outside(edge));
end
