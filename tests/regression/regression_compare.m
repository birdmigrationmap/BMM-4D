function difference = regression_compare(reference, current)
%REGRESSION_COMPARE Relative difference of each output of a regression case.
%   REFERENCE and CURRENT are structs of numeric, logical or sparse arrays,
%   or cells of them. DIFFERENCE is a table with, per output, whether the
%   sizes match and the relative difference norm(current - reference) /
%   norm(reference) (Frobenius norm over all elements). NaN must be at the
%   same places; the difference is Inf when they are not or the sizes
%   differ.

names = string(fieldnames(reference));
same_size = false(numel(names), 1);
relative = inf(numel(names), 1);
for i = 1:numel(names)
    a = reference.(names(i));
    b = current.(names(i));
    if ~iscell(a)
        a = {a};
        b = {b};
    end
    same_size(i) = numel(a) == numel(b) && all(cellfun(@(x, y) isequal(size(x), size(y)), a, b));
    if ~same_size(i)
        continue
    end
    squared_difference = 0;
    squared_reference = 0;
    same_nan = true;
    for k = 1:numel(a)
        x = double(a{k});
        y = double(b{k});
        if ~issparse(x)
            known = ~isnan(x);
            same_nan = same_nan && isequal(known, ~isnan(y));
            x(~known) = 0;
            y(isnan(y)) = 0;
        end
        squared_difference = squared_difference + norm(y - x, "fro") ^ 2;
        squared_reference = squared_reference + norm(x, "fro") ^ 2;
    end
    if same_nan
        relative(i) = sqrt(squared_difference) / max(sqrt(squared_reference), eps);
    end
end
difference = table(names, same_size, relative, VariableNames=["output" "same_size" "relative_difference"]);
end
