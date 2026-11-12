function volume = as_apply_aspect(volume, aspect)
%AS_APPLY_ASPECT Divide a cleaned volume by the shared bird aspect function (dev).
%   VOLUME = AS_APPLY_ASPECT(VOLUME, ASPECT) divides the density (and sd)
%   of every gate of VOLUME (as saved by bmm4d.radar.clean_volume, without
%   any shape or aspect correction) by as_aspect_template at the gate's
%   azimuth and beam-centre height above the antenna (4/3 earth), with the
%   volume's own velocity profile VOLUME.uv giving the track at that
%   height. The template has an azimuthal mean of one, so no rescaling.
%   ASPECT is stored in VOLUME.aspect (the factor can be undone).

for i_elev = 1:numel(volume.dens.sweeps)
    sweep = volume.dens.sweeps(i_elev);
    [~, height_m] = slant2ground(sweep.range(:), sweep.elev);
    [az, height] = ndgrid(sweep.az(:)', height_m);   % [n_az n_range]
    template = as_aspect_template(az', height', volume.uv, aspect);   % [n_range n_az]
    volume.dens.sweeps(i_elev).data = sweep.data ./ template;
    if isfield(sweep, "sd") && ~isempty(sweep.sd)
        volume.dens.sweeps(i_elev).sd = sweep.sd ./ template;
    end
end
volume.aspect = aspect;
end
