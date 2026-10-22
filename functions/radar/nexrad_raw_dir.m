function raw_dir = nexrad_raw_dir()
%NEXRAD_RAW_DIR Folder where raw NEXRAD volumes are cached.
%   Raw volumes can always be downloaded again, so they are kept outside
%   the synced project folder when local_settings.json (project root,
%   gitignored) sets "nexrad_raw_dir". Without it, the project folder
%   data/raw_scans/nexrad is used. Volumes are stored as
%   <raw_dir>/<RADAR>_<yyyymmdd>/<file>.

raw_dir = fullfile("data", "raw_scans", "nexrad");
settings_file = "local_settings.json";
if isfile(settings_file)
    settings = jsondecode(fileread(settings_file));
    if isfield(settings, "nexrad_raw_dir")
        raw_dir = string(settings.nexrad_raw_dir);
    end
end
if ~isfolder(raw_dir)
    mkdir(raw_dir)
end
end
