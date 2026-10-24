function project_root = setup_paths(project_root)
% Set the project root, load shared code, and normalize the working directory.

if nargin == 0 || strlength(string(project_root)) == 0
    project_root = fileparts(mfilename("fullpath"));
end

project_root = char(project_root);

while ~isfolder(fullfile(project_root, "toolbox")) || ~isfile(fullfile(project_root, "setup_paths.m"))
    parent_root = fileparts(project_root);
    if strcmp(parent_root, project_root)
        error("setup_paths:ProjectRootNotFound", "Could not find the BMM-4D project root.")
    end
    project_root = parent_root;
end

% Run active scripts from the project root, where figure and export folders
% are relative to.
cd(project_root)

% External code: folders set in local_settings.json (not tracked; see
% local_settings.example.json), relative to the project root or absolute.
settings = struct();
settings_file = fullfile(project_root, "local_settings.json");
if isfile(settings_file)
    settings = jsondecode(fileread(settings_file));
end
external_dir = @(name, default) absolute_folder(project_root, setting_or(settings, name, default));

% wsrlib (required): NEXRAD reading (rsl2mat), VVP, MatConvNet for MistNet
wsrlib_dir = external_dir("wsrlib_dir", fullfile("external", "wsrlib"));
if ~isfolder(wsrlib_dir)
    error("setup_paths:MissingWsrlib", ...
        "wsrlib not found at %s. Clone https://github.com/darkecology/wsrlib and set ""wsrlib_dir"" in local_settings.json.", ...
        wsrlib_dir)
end
addpath(wsrlib_dir);
wsrlib_setup();

% Regularization Tools (Hansen; required by compute_vp): cgsvd, tikhonov,
% L-curve corner
regtools_dir = external_dir("regtools_dir", fullfile("external", "regtools"));
if isfolder(regtools_dir)
    addpath(genpath(regtools_dir));
end

% FastGaussianSimulation (optional): synthetic Gaussian truth fields
fgs_dir = external_dir("fast_gaussian_simulation_dir", fullfile("external", "FastGaussianSimulation"));
if isfolder(fgs_dir)
    addpath(fgs_dir);
end

% Package code: the bmm4d package (toolbox/, resources/mpackage.json),
% installed in place for this MATLAB session so edits apply at once.
toolbox_dir = fullfile(project_root, "toolbox");
if isempty(mpmlist("bmm4d"))
    mpminstall(toolbox_dir, InPlace=true, Temporary=true, PathPosition="beginning", Prompt=false, Verbosity="quiet");
end

end

function value = setting_or(settings, name, default)
value = default;
if isfield(settings, name)
    value = settings.(name);
end
end

function folder = absolute_folder(project_root, folder)
folder = string(folder);
if ~startsWith(folder, ["/" "~"]) && ~contains(folder, ":")
    folder = fullfile(project_root, folder);
end
end
