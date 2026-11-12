%% Aspect A: bin the us_east night for the overlap and single-radar analyses (dev)
% Every volume of the us_east night (21-22 Oct 2026, 107 radars, 26 frames
% every 30 min) is read three ways (as_frame_bins): raw (us_east_fix with
% its symmetric correction undone), ell (production cache 9c2b4ef2:
% elliptical shape correction of the paper) and fix (us_east_fix: the
% per-volume symmetric fit of studies/artifacts). Gates are moved to the
% frame time with the radar's own velocity profile and binned on
% 4 km x 4 km x 100 m cells; cells seen by two radars or more are kept
% (overlap). Azimuthal sector sums of each volume are kept for the
% single-radar harmonics. One file per frame (~1 min each); existing files
% are skipped, so the loop can be run in chunks.
% Output: <data>/aspect/bins/frame_XX.mat (bins, volumes).

%% Load packages and paths
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "aspect", "functions"))
data_root = bmm4d.util.local_data_dir();
manifest = load(fullfile(data_root, "prod", "us_east", "manifest.mat"));
fix_dir = fullfile(data_root, "prod", "volumes", "us_east_fix");
bin_dir = fullfile(data_root, "aspect", "bins");
if ~isfolder(bin_dir)
    mkdir(bin_dir)
end

%% Bin every frame
for i_frame = 1:numel(manifest.times)
    out_file = fullfile(bin_dir, sprintf("frame_%02d.mat", i_frame));
    if isfile(out_file)
        continue
    end
    as_frame_bins(manifest, i_frame, fix_dir, out_file);
    fprintf("frame %d done\n", i_frame)
end
