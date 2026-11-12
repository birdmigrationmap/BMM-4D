function sums = as_sector_corrected(volume, aspect)
%AS_SECTOR_CORRECTED Sector sums of a volume divided by the aspect factor (dev).
%   SUMS = AS_SECTOR_CORRECTED(VOLUME, ASPECT) divides the raw sector sums
%   of as_frame_bins ([36 sectors x 12 sweeps x 5 range bands]) by
%   as_aspect_template at the sector centre and the cell's mean height
%   above the antenna, with the volume's own velocity profile.

counts = volume.sectors.n;
height_m = volume.sectors.height ./ counts;
height_m(~isfinite(height_m)) = 0;
az = repmat((5:10:355)', 1, size(counts, 2), size(counts, 3));
template = as_aspect_template(az, height_m, volume.uv, aspect);
sums = volume.sectors.raw ./ template;
end
