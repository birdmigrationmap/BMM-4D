function density = gd_old_cleaning(dbz, range_m, az_deg, elev_deg)
%GD_OLD_CLEANING Density of one sweep with the gate-by-gate cleaning of BMM-4D before the estimator.
%   DENSITY = GD_OLD_CLEANING(DBZ, RANGE_M, AZ_DEG, ELEV_DEG) applies to a
%   sweep of reflectivity DBZ [n_range n_az] (NaN: no echo) the cleaning of
%   the committed BMM-4D version (git HEAD f2a49f7..ca7928d, bmm4d.config):
%   read_rd with density "median" (MAD outliers at 3, the 5 x 5 continuity
%   filter, no-echo gates as zero density within 150 km, conversion with
%   11 cm^2 per bird, 2 km moving median), then clean_volume's cap at 3000
%   birds km^-3 and removal of echo patches smaller than 25 km^2. Rain
%   removal, beam blockage and the shape correction are left out, as for the
%   estimator in GD_I, so that only the step from gates to density differs.
%   The toolbox filters are unchanged since that version.

parms = struct("removeOutliers", 3, "thr_dens", 1e9, "noEchoRange", 150000, "smoothing", 2000, ...
    "maxRange", 150000, "plot_opts", "", "filename", "", "export_dir", "");
gate = struct("elev", elev_deg, "range", range_m(:), "az", az_deg(:)', "data", dbz, "unit", "dBZ", ...
    "no_echo", isnan(dbz), "rain_mask", false(size(dbz)));
scan = struct("gates", gate);
scan = bmm4d.radar.filter.outlier_mad(scan, parms);
scan = bmm4d.radar.filter.outlier_continuity(scan, parms);
scan = bmm4d.radar.filter.no_echo(scan, parms);
scan = bmm4d.radar.filter.to_density(scan, parms);
scan = bmm4d.radar.filter.smooth_median(scan, parms);
density = scan.gates.data;

% clean_volume of that version: density cap, then small echo patches
density(density > 3000) = nan;
gate_area_km2 = (range_m(:) * deg2rad(median(diff(az_deg))) * median(diff(range_m))) / 1e6;
patches = bwconncomp(density > 0, 8);
patch_area = cellfun(@(ids) sum(gate_area_km2(mod(ids - 1, numel(range_m)) + 1)), patches.PixelIdxList);
small = patches.PixelIdxList(patch_area < 25);
density(vertcat(small{:})) = nan;
end
