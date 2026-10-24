function raw_dir = nexrad_raw_dir()
%NEXRAD_RAW_DIR Folder where raw NEXRAD volumes are cached.
%   Raw volumes can always be downloaded again, so they are kept outside
%   the synced project folder when local_settings.json sets
%   "nexrad_raw_dir". Without it, raw_scans/nexrad in the data folder
%   (data_path) is used. Volumes are stored as
%   <raw_dir>/<RADAR>_<yyyymmdd>/<file>.

raw_dir = local_folder("nexrad_raw_dir", data_path("raw_scans", "nexrad"));
if ~isfolder(raw_dir)
    mkdir(raw_dir)
end
end
