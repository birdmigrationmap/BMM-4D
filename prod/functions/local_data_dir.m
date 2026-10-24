function data_root = local_data_dir(project_root)
%LOCAL_DATA_DIR Root folder for large generated products.
%   Large intermediate products (radar caches, case and inversion files)
%   can be rebuilt from the raw data, so they are kept outside the synced
%   project folder when local_settings.json (project root, gitignored) sets
%   "local_data_dir". Without it, PROJECT_ROOT is returned and products stay
%   in the project folders.

data_root = string(project_root);
settings_file = fullfile(project_root, "local_settings.json");
if isfile(settings_file)
    settings = jsondecode(fileread(settings_file));
    if isfield(settings, "local_data_dir")
        data_root = string(settings.local_data_dir);
    end
end
end
