function mask = static_clutter_mask(raw_files, settings)
%STATIC_CLUTTER_MASK Gates of one radar that stand out at the same place all night.
%   For the lowest reflectivity sweep of each RAW_FILES volume, flags gates
%   more than SETTINGS.threshold_db above their local median
%   (local_contrast, SETTINGS.window_m). Gates flagged in more than
%   SETTINGS.min_fraction of the volumes are static clutter: wind farms,
%   towers, persistent ground clutter. Birds and rain move, so they are not
%   flagged at the same gate all night. The mask is grown by one gate in
%   range and azimuth. MASK has fields range, az (vectors) and flag
%   ([range az] logical), or is empty when no volume could be read.

mask = [];
count = [];
n_volume = 0;
for k = 1:numel(raw_files)
    try
        scan = rsl2mat_wrap(char(raw_files(k)), struct("max_elev", 1));
        sweeps = unique_elev_sweeps(scan, 'dz');
        [dz, range_m, az] = sweep2mat(sweeps(1));
    catch
        continue
    end
    keep = range_m <= settings.max_range_m;
    dz = dz(keep, :);
    if isempty(count)
        count = zeros(size(dz));
        mask_range = range_m(keep);
        mask_az = az;
    end
    if ~isequal(size(dz), size(count))
        continue   % other sweep geometry (VCP change): skip
    end
    contrast = local_contrast(10 .^ (dz / 10), mask_range, mask_az, settings.window_m, 0.1);
    count = count + (contrast > settings.threshold_db & ~isnan(dz));
    n_volume = n_volume + 1;
end
if n_volume == 0
    return
end
flag = count / n_volume > settings.min_fraction;
flag = imdilate(flag(:, [end 1:end 1]), true(3));   % wrap around north
flag = flag(:, 2:end-1);
mask = struct("range", mask_range(:), "az", mask_az(:), "flag", flag, "n_volume", n_volume);
end
