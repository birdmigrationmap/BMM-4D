function n_cleaned = ar_clean_cases(cases, project_root, varargin)
%AR_CLEAN_CASES Clean several cases in turn into one cache, for a batch job (dev).
%   N_CLEANED = AR_CLEAN_CASES(CASES, PROJECT_ROOT, ...) runs
%   ar_clean_volumes on each case of the cell array CASES in order, with
%   the same options (e.g. "dev_folder", "cache_dir", "aspect"). Cases that
%   share volumes share the cache, so small test cases listed first are
%   ready early and a larger case after them only cleans the rest.
%   N_CLEANED [n_case 1] counts the volumes cleaned for each case.

n_cleaned = zeros(numel(cases), 1);
for i_c = 1:numel(cases)
    n_cleaned(i_c) = ar_clean_volumes(cases{i_c}, project_root, varargin{:});
end
end
