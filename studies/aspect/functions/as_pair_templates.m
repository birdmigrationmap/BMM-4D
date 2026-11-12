function [log_ta, log_tb] = as_pair_templates(D, volumes, aspect)
%AS_PAIR_TEMPLATES Aspect factors of both radars of each overlap pair (dev).
%   [LOG_TA, LOG_TB] = AS_PAIR_TEMPLATES(D, VOLUMES, ASPECT) evaluates
%   log as_aspect_template for radar a and radar b of each pair of
%   as_overlap_design, as the cleaning would: each radar with its own
%   velocity profile (VOLUMES of as_frame_bins, by frame and radar) at the
%   cell's height above its own antenna. Pairs whose volume is not in
%   VOLUMES keep 0.

az = {mod(D.psi_a + D.track, 360), mod(D.psi_b + D.track, 360)};
radar = {D.radar_a, D.radar_b};
log_t = {zeros(size(D.y)), zeros(size(D.y))};
volume_key = [volumes.frame] * 1000 + [volumes.radar];
for side = 1:2
    % Rows of each volume (frame and radar) on this side of the pairs
    [key, ~, group] = unique(D.frame * 1000 + radar{side});
    rows_of = accumarray(group, (1:numel(D.y))', [], @(rows) {rows});
    [found, i_volume] = ismember(key, volume_key);
    for k = find(found)'
        vol = volumes(i_volume(k));
        rows = rows_of{k};
        log_t{side}(rows) = log(as_aspect_template(az{side}(rows), D.height_asl_m(rows) - vol.height, vol.uv, aspect));
    end
end
[log_ta, log_tb] = deal(log_t{:});
end
