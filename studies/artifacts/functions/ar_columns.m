function column = ar_columns(S, field)
%AR_COLUMNS Column density maps of a field on the active cells (dev).
%   COLUMN = AR_COLUMNS(S, FIELD) is [nx ny n_time] in birds km^-2: the
%   sum over the column of density x air share of the cell x dz, negative
%   densities counted as zero, NaN outside active columns.

n_time = size(field, 2);
column = nan([S.g.sz(1:2) n_time]);
active_column = any(S.active_mask_3d, 3);
for i_t = 1:n_time
    frame = zeros(S.g.sz);
    frame(S.active_id) = max(double(field(:, i_t)), 0) .* S.cell_air;
    column_t = sum(frame, 3) * S.g.dz / 1000;
    column_t(~active_column) = nan;
    column(:, :, i_t) = column_t;
end
end
