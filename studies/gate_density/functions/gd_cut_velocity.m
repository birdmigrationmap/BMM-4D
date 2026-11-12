function vel = gd_cut_velocity(level2, i_e)
%GD_CUT_VELOCITY Radial velocity on the grid of cut I_E of a volume.
%   LEVEL2 is the output of READ_LEVEL2_MOMENTS with REF, VEL and CFP. Batch
%   cuts carry VEL themselves. Split cuts (the lowest elevations) do not: the
%   surveillance pass has REF and CFP, the following Doppler pass at the same
%   elevation has VEL; its nearest radial and gate are taken. Returns
%   [n_gate n_az] on the REF gates of the cut, or [] if no velocity is found.

cut = level2(i_e);
if ~isempty(cut.VEL)
    vel_cut = cut;
else
    elevation = [level2.elev];
    elevation_number = [level2.elevation_number];
    has_vel = arrayfun(@(s) ~isempty(s.VEL), level2);
    has_cfp = arrayfun(@(s) ~isempty(s.CFP), level2);
    % The Doppler pass follows its surveillance pass (SAILS repeats the lowest cut)
    candidates = find(has_vel(:) & ~has_cfp(:) & abs(elevation(:) - cut.elev) < 0.25 & elevation_number(:) > cut.elevation_number);
    if isempty(candidates)
        vel = [];
        return
    end
    vel_cut = level2(candidates(1));
end
[~, i_az] = min(abs(mod(cut.az(:) - vel_cut.az(:)' + 180, 360) - 180), [], 2);
i_gate = interp1(vel_cut.VEL_range(:), 1:numel(vel_cut.VEL_range), cut.REF_range(:), "nearest");
vel = nan(numel(cut.REF_range), numel(cut.az));
vel(~isnan(i_gate), :) = vel_cut.VEL(i_gate(~isnan(i_gate)), i_az);
end
