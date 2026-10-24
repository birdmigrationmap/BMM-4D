function project_root = setup_paths(project_root)
% Set the project root, load shared code, and normalize the working directory.

if nargin == 0 || strlength(string(project_root)) == 0
    project_root = fileparts(mfilename("fullpath"));
end

project_root = char(project_root);

while ~isfolder(fullfile(project_root, "functions")) || ~isfile(fullfile(project_root, "setup_paths.m"))
    parent_root = fileparts(project_root);
    if strcmp(parent_root, project_root)
        error("setup_paths:ProjectRootNotFound", "Could not find the BMM-4D project root.")
    end
    project_root = parent_root;
end

% Run active scripts from the project root so relative data paths resolve.
cd(project_root)

% Load shared external code.
addpath('/Users/rafnuss/Documents/GitHub/wsrlib');
wsrlib_setup();

% Load optional Gaussian-field simulator if it is available locally.
fgs_paths = [
    fullfile(project_root, "external", "FastGaussianSimulation")
    "/Users/rafnuss/Documents/Github/FastGaussianSimulation"
    "/Users/rafnuss/Documents/GitHub/FastGaussianSimulation"
];
for i_path = 1:numel(fgs_paths)
    if isfolder(fgs_paths(i_path))
        addpath(char(fgs_paths(i_path)));
        break
    end
end

% Load local project code.
addpath(genpath(fullfile(project_root, "functions")));

end
