function c_lim = get_case_clim(varargin)

values = [];
for i_arg = 1:nargin
    if isempty(varargin{i_arg})
        continue
    end
    values = [values; varargin{i_arg}(:)]; %#ok<AGROW>
end

values = values(isfinite(values));
c_lim = [0 max(quantile(values, 0.99), eps)];
end
