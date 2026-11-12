%% Aspect C: polarimetric sample, the flight axis and the two polarizations (dev)
% Two frames of the us_east night (03:30 and 06:30 UTC, ~200 volumes) are
% read from the raw Level-II files (as_frame_polar, ~3.5 min per frame):
% Z_H, ZDR and rhoHV at the gates kept by the cleaning (rain and clutter
% removed), Z_V = Z_H - ZDR.
% - The heading from rhoHV: rhoHV is lowest tail-on and highest head-on
%   (Stepanian et al. 2016); the minimum of its first azimuthal harmonic
%   (density-weighted, 10-60 km) against the track.
% - ZDR against the side-on axis: Z_H / Z_V sees the same birds, so its
%   azimuthal harmonic is free of real gradients; if only Z_H were
%   modulated (Rayleigh body, horizontal), ZDR would carry the amplitude.
% - Overlap fits for the density, Z_H and Z_V, and with four flight axes:
%   the track at the cell height (two radars' profiles), the rhoHV heading
%   (mean of the two volumes), a fixed 190 deg, and the volume track at
%   200-1500 m (as studies/artifacts ar_clean_volume).
% Outputs: <data>/aspect/polar/frame_XX.mat, <data>/aspect/polar_results.mat,
% figures/AS_axis_polarization.csv.

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"))
data_root = bmm4d.util.local_data_dir();
manifest = load(fullfile(data_root, "prod", "us_east", "manifest.mat"));
fix_dir = fullfile(data_root, "prod", "volumes", "us_east_fix");
polar_dir = fullfile(data_root, "aspect", "polar");
figure_dir = fullfile(project_root, "studies", "aspect", "figures");
if ~isfolder(polar_dir)
    mkdir(polar_dir)
end
frames = [9 15];   % 03:30 and 06:30 UTC

%% Read the raw polarimetric fields of the two frames
for i_frame = frames
    out_file = fullfile(polar_dir, sprintf("frame_%02d.mat", i_frame));
    if ~isfile(out_file)
        as_frame_polar(manifest, i_frame, fix_dir, out_file);
    end
end
loaded = arrayfun(@(i_frame) load(fullfile(polar_dir, sprintf("frame_%02d.mat", i_frame))), frames);
polar_volumes = [loaded.volumes];
polar_pairs = vertcat(as_overlap_pairs(loaded(1).bins), as_overlap_pairs(loaded(2).bins));

%% Heading from rhoHV and ZDR harmonics of each volume
P = as_polar_summary(polar_volumes);
ok = isfinite(P.drift);
circular_sd = @(deg) rad2deg(sqrt(-2 * log(abs(mean(exp(1i * deg2rad(deg)))))));
polar_summary = struct("n", nnz(ok), "drift_median", median(P.drift(ok)), "drift_iqr", iqr(P.drift(ok)), ...
    "drift_circular_sd", circular_sd(P.drift(ok)), "heading_circular_sd", circular_sd(P.heading_rho(ok)), ...
    "track_circular_sd", circular_sd(P.track(ok)), "zdr_c2_mean", mean(P.zdr_c2(ok)), ...
    "zdr_c2_se", std(P.zdr_c2(ok)) / sqrt(nnz(ok)), "zdr_c1_mean", mean(P.zdr_c1(ok)), ...
    "zh_c2_mean", mean(P.zh_c2(ok)), "zv_c2_mean", mean(P.zv_c2(ok)), ...
    "corr_zdr_zh_c2", corr(P.zdr_c2(ok), P.zh_c2(ok)));
disp(polar_summary)

%% Overlap: density, Z_H and Z_V, four flight axes
heading = nan(max(frames), numel(manifest.radars));
track = heading;
heading(sub2ind(size(heading), P.frame, P.radar)) = P.heading_rho;
track(sub2ind(size(track), P.frame, P.radar)) = P.track;
mean_direction = @(a, b) mod(atan2d(sind(a) + sind(b), cosd(a) + cosd(b)), 360);
id_a = sub2ind(size(heading), double(polar_pairs.frame), double(polar_pairs.radar_a));
id_b = sub2ind(size(heading), double(polar_pairs.frame), double(polar_pairs.radar_b));
axes_list = {[], mean_direction(heading(id_a), heading(id_b)), 190 * ones(height(polar_pairs), 1), ...
    mean_direction(track(id_a), track(id_b))};
axis_names = ["track at cell height" "rhoHV heading" "fixed 190 deg" "volume track 200-1500 m"];
axis_polarization = table();
for quantity = ["raw" "zh" "zv"]
    for i_axis = 1:numel(axes_list)
        D = as_overlap_design(polar_pairs, quantity, axes_list{i_axis});
        T = as_fit_strata(D, {true(size(D.y))}, quantity + " / " + axis_names(i_axis));
        axis_polarization = [axis_polarization; T]; %#ok<AGROW>
    end
end
disp(axis_polarization(:, ["stratum" "n" "c2" "c2_se" "s2" "c1" "c1_se" "c4"]))
writetable(axis_polarization, fullfile(figure_dir, "AS_axis_polarization.csv"))

%% Save
save(fullfile(data_root, "aspect", "polar_results.mat"), "P", "polar_summary", "axis_polarization")
