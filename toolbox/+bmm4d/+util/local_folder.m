function folder = local_folder(name, default)
%LOCAL_FOLDER Machine-specific folder from local_settings.json.
%   FOLDER = LOCAL_FOLDER(NAME, DEFAULT) returns field NAME of
%   local_settings.json in the repository root (not tracked; see
%   local_settings.example.json), or DEFAULT when the file or the field is
%   missing. Relative folders are taken from the repository root.

% The package locates its own repository from this file
% (toolbox/+bmm4d/+util/local_folder.m), so it works from any current folder.
repository_root = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
folder = string(default);
settings_file = fullfile(repository_root, "local_settings.json");
if isfile(settings_file)
    settings = jsondecode(fileread(settings_file));
    if isfield(settings, name)
        folder = string(settings.(name));
    end
end
if ~startsWith(folder, ["/" "~"]) && ~contains(folder, ":")
    folder = fullfile(repository_root, folder);
end
end
