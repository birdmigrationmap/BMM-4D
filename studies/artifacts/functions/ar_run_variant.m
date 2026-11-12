function solve_info = ar_run_variant(night_name, variant, holdout, range_m)
%AR_RUN_VARIANT One prior variant of the night-set test (dev, cluster job).
%   SOLVE_INFO = AR_RUN_VARIANT(NIGHT_NAME, VARIANT, HOLDOUT, RANGE_M)
%   solves the night with ar_invert_variant and the prior of VARIANT:
%   - "V0": white prior + smoothness (production);
%   - "V1": correlated prior 100 km + smoothness (fix 3);
%   - "GP3D": correlated prior with horizontal and vertical ranges
%     RANGE_M = [l_h l_v] (m), no smoothness;
%   - "B": one Matern term (PRIOR_OPTIONS.md, option B) with ranges
%     RANGE_M, no smoothness.
%   GP3D and B runs are saved as <VARIANT><l_h in km> (e.g. B20).
%   The ranges are a first guess read from the anomaly correlation of
%   three nights (range_fit.mat of ar_night_stage1, 12 Nov 2026): most of
%   the horizontal correlation is gone by 15-25 km, with a weaker tail to
%   50-70 km on some nights; vertically it holds to ~250 m and is mostly
%   gone by 400-650 m. The single-range fits themselves were unstable
%   (Matern 10-43 km) and are not used; a fit by night, time of night and
%   day of year is left for later.

if nargin < 4
    range_m = [];
end
smoothness = ar_night_config(night_name).inversion.smoothness;
label = variant;
switch variant
    case "V0"
        prior_options = {"smoothness", smoothness};
    case "V1"
        prior_options = {"smoothness", smoothness, "prior_correlation_m", 100000};
    case "GP3D"
        prior_options = {"smoothness", 0, "prior_correlation_m", range_m(1), "prior_correlation_vertical_m", range_m(2)};
        label = variant + round(range_m(1) / 1000);
    case "B"
        prior_options = {"smoothness", 0, "prior_matern_m", range_m};
        label = variant + round(range_m(1) / 1000);
    otherwise
        error("ar_run_variant:Variant", "Unknown variant %s.", variant)
end
solve_info = ar_invert_variant(night_name, label, prior_options, holdout);
end
