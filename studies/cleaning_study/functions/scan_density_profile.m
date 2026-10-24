function [profile, n_gate] = scan_density_profile(scan, z_edges)
%SCAN_DENSITY_PROFILE Mean bird density in height bins from a read_rd scan.
%   PROFILE(k) is the mean density (birds km^-3) of all retained gates whose
%   beam centre lies between Z_EDGES(k) and Z_EDGES(k+1) above the radar.
%   N_GATE(k) is the number of gates averaged.

density_sum = zeros(numel(z_edges) - 1, 1);
n_gate = zeros(numel(z_edges) - 1, 1);
for i_elev = 1:numel(scan.dens.sweeps)
    sweep = scan.dens.sweeps(i_elev);
    [~, gate_height] = slant2ground(sweep.range(:), sweep.elev);
    gate_height = repmat(gate_height, 1, size(sweep.data, 2));
    valid = ~isnan(sweep.data);
    bin_id = discretize(gate_height(valid), z_edges);
    density = sweep.data(valid);
    keep = ~isnan(bin_id);
    density_sum = density_sum + accumarray(bin_id(keep), density(keep), size(density_sum));
    n_gate = n_gate + accumarray(bin_id(keep), 1, size(n_gate));
end
profile = density_sum ./ n_gate;
end
