function file = data_path(varargin)
%DATA_PATH Path inside the data folder.
%   FILE = DATA_PATH("reference", "radar.csv") returns the path of that
%   file or folder in the data folder: "data_dir" in local_settings.json,
%   or data/ in the repository root. Reference layers, beam kernels and
%   older raw scans are read from there, whatever the current folder.

file = fullfile(local_folder("data_dir", "data"), varargin{:});
end
