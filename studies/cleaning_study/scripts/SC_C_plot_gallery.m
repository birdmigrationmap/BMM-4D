%% Scan cleaning study: review gallery
% Export one review figure per selected case.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "cleaning_study", "functions"))
data_dir = fullfile(project_root, "studies", "cleaning_study", "data");
figure_dir = fullfile(project_root, "studies", "cleaning_study", "figures", "gallery");
if ~exist(figure_dir, "dir")
    mkdir(figure_dir)
end
load(fullfile(data_dir, "cases.mat"), "cases")
load(fullfile(data_dir, "screening.mat"), "screening_controls")

%% Plot every case --------------------------------------------------------
parfor i_c = 1:height(cases)
    output_file = fullfile(figure_dir, cases.case_id(i_c) + "_" + cases.group(i_c) + "_" + cases.radar(i_c) + ".png");
    try
        plot_cleaning_case(cases(i_c, :), screening_controls.method_names, ...
            screening_controls.read_options, screening_controls.z_edges, output_file);
    catch err
        warning("SC_C:PlotFailed", "%s: %s", cases.case_id(i_c), err.message)
    end
end
fprintf("[SC/C] exported %d figures to %s\n", height(cases), figure_dir)
