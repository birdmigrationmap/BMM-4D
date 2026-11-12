%% Artefacts B: the anisotropy (shape) correction and azimuthal wedges (dev, fix 1)
% Screenshot 5 of the us_east video: a wedge of low column NE of KDMX at
% 05:30-07:30 UTC, between dense areas NW and SE of the radar. KDMX's
% observations are 0.55-0.75 of the other sectors at every elevation
% (0.5-3.1 deg) and from 5 km, so it is not blockage. Re-cleaning the
% 06:30 volume without the shape correction gives the opposite pattern
% (densest NW, lowest S-SE, NE about the median): the correction divides
% each volume by an elliptical template fitted to the azimuthal density
% (different front and back axes) along the rhoHV axis, and here the
% template, fitted to a real NW-SE gradient, makes a bow-tie.
% The same 28 Midwest radars and 7 frames (04:30-07:30 UTC) cleaned three
% ways (AR_A_cases.m), solved with the production settings:
% - ar_midwest: elliptical shape correction (production);
% - ar_midwest_noshape: no correction;
% - ar_midwest_fix: symmetric aspect correction (ar_clean_volume:
%   1 + a cos 2(az - track - 90 deg), track from the volume's profile),
%   with the beam-blockage fix (flat here, little effect).
% Scores: wedges (ar_wedge_score: deepest 10 deg sector at 10-70 km
% against the sectors 30 deg either side), azimuthal harmonics of each
% radar's observations at 10-60 km (first: gradient across the disc;
% second along the side-on axis: the aspect effect), data misfit (each
% against its own observations, so not a ranking).
% Outputs: figures/AR_B_shape.png (maps around KDMX), AR_B_shape.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
box = [-100 -85 36 46];
no_shape = {"read_options", {"shape", false}};
cases = {ar_test_case("ar_midwest", box, times), ar_test_case("ar_midwest_noshape", box, times, no_shape{:}), ...
    ar_test_case("ar_midwest_fix", box, times, no_shape{:})};
names = ["shape correction (production)" "no correction" "symmetric aspect correction"];

%% Solves and scores -------------------------------------------------------------
results = cell(1, 3);
scores = table();
for i_c = 1:3
    S = ar_load(cases{i_c});
    results{i_c} = ar_solve(S, names(i_c));

    % Azimuthal harmonics of the observations at 10-60 km, per radar and frame
    harmonics = ar_obs_harmonics(S);
    first = harmonics.first;
    side_on = harmonics.side_on;

    wedge = ar_wedge_score(S, results{i_c}.column);
    kdmx = S.radar_selected.name == "KDMX";
    scores = [scores; table(names(i_c), wedge.median, wedge.share_below_08, min(wedge.per_radar(kdmx, :)), ...
        median(first, "all", "omitnan"), median(side_on, "all", "omitnan"), mean(side_on(isfinite(side_on)) > 0), ...
        results{i_c}.misfit, results{i_c}.misfit_trusted, ...
        VariableNames=["cleaning" "wedge_median" "wedge_share_below_08" "kdmx_deepest" "first_harmonic" ...
        "side_on_harmonic" "side_on_share_positive" "misfit" "misfit_30_120km"])]; %#ok<AGROW>
    if i_c == 1
        S_map = S;
    end
end
disp(scores)
writetable(scores, fullfile(figure_dir, "AR_B_shape.csv"))

%% Maps around KDMX --------------------------------------------------------------
ar_plot_maps(S_map, results, [3 5 7], fullfile(figure_dir, "AR_B_shape.png"), bounds=[-96.5 -91 39.8 43.8]);
