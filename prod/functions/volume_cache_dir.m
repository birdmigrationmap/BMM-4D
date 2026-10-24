function cache_dir = volume_cache_dir(c)
%VOLUME_CACHE_DIR Folder of cleaned volumes for the case's cleaning settings.
%   Cleaned volumes are shared by all cases and depend only on the cleaning
%   options, so the folder name carries a short hash of read_options and the
%   velocity-profile settings. Changing any of them starts a new cache.

settings_text = jsonencode(struct("read_options", {c.read_options}, "vp", c.vp, "filter", c.filter));
digest = java.security.MessageDigest.getInstance("MD5");
hash_bytes = typecast(int8(digest.digest(uint8(char(settings_text)))), "uint8");
tag = lower(string(reshape(dec2hex(hash_bytes(1:4))', 1, [])));
cache_dir = fullfile(fileparts(c.data_dir), "volumes", tag);
if ~isfolder(cache_dir)
    mkdir(cache_dir)
    writelines(settings_text, fullfile(cache_dir, "settings.json"))
end
end
