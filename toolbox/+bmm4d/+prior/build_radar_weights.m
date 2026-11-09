function interpolation_weights = build_radar_weights(g, radar_scans, options)
%BUILD_RADAR_WEIGHTS Radar-to-grid weights [nxy n_radar] of a radar profile field.
%   W = BMM4D.PRIOR.BUILD_RADAR_WEIGHTS(G, RADAR_SCANS, interpolation=...)
%   gives the horizontal weights that spread one value per radar to every
%   grid column (W * values), as used by the density prior, the velocity
%   prior, the advection steps and the web map export:
%   - "gaussian": smooth Gaussian weights of length scale length_scale_m
%     (build_gaussian_weights), no seams and no peak at each radar;
%   - "idw": inverse-distance weights of the idw_neighbors nearest radars
%     with power idw_power (build_horizontal_idw_weights), as in the paper.

arguments
    g (1,1) struct
    radar_scans cell
    options.interpolation (1,1) string {mustBeMember(options.interpolation, ["gaussian" "idw"])} = "idw"
    options.length_scale_m (1,1) double = 150000
    options.idw_neighbors (1,1) double = 3
    options.idw_power (1,1) double = 2.25
end

if options.interpolation == "gaussian"
    interpolation_weights = bmm4d.prior.build_gaussian_weights(g, radar_scans, options.length_scale_m);
else
    interpolation_weights = bmm4d.prior.build_horizontal_idw_weights(g, radar_scans, ...
        n_neighbors=options.idw_neighbors, power=options.idw_power);
end
end
