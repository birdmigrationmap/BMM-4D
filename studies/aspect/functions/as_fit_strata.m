function T = as_fit_strata(D, strata, labels)
%AS_FIT_STRATA Overlap aspect harmonics in subsets of the pairs (dev).
%   T = AS_FIT_STRATA(D, STRATA, LABELS) fits as_overlap_fit (harmonics
%   c2 s2 c1 s1 c4 and the beam slopes) on each subset STRATA{i} (logical)
%   and returns one row per subset: label, n, n_pair, the harmonics and
%   their clustered standard errors, and the mismatch sd.

harmonics = ["c2" "s2" "c1" "s1" "c4"];
columns = [harmonics D.names(D.beam_columns)];
rows = cell(numel(strata), 1);
for i = 1:numel(strata)
    fit = as_overlap_fit(D, columns, strata{i});
    [~, id] = ismember(harmonics, fit.names);
    beta = nan(1, numel(harmonics));
    se = beta;
    beta(id > 0) = fit.beta(id(id > 0));
    se(id > 0) = fit.se(id(id > 0));
    rows{i} = [table(string(labels(i)), fit.n, fit.n_pair, fit.sigma_mismatch, ...
        VariableNames=["stratum" "n" "n_pair" "sigma_mismatch"]), ...
        array2table([beta se], VariableNames=[harmonics harmonics + "_se"])];
end
T = vertcat(rows{:});
end
