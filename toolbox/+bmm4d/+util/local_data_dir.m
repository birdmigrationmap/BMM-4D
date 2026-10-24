function data_root = local_data_dir()
%LOCAL_DATA_DIR Root folder for large generated products.
%   Large intermediate products (radar caches, case and inversion files)
%   can be rebuilt from the raw data, so they are kept outside the synced
%   project folder when local_settings.json sets "local_data_dir". Without
%   it, the repository root is used and products stay in the project
%   folders.

data_root = bmm4d.util.local_folder("local_data_dir", ".");
end
