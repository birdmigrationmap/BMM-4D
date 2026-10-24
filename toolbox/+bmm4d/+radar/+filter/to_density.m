function scan = to_density(scan, ~)
%TO_DENSITY Convert reflectivity (dBZ) to bird density (birds km^-3).
%   Reflectivity factor (mm^6 m^-3) to reflectivity (cm^2 km^-3), divided
%   by a bird cross-section of 11 cm^2. Filters after this one work on
%   density.

for i_elev = 1:numel(scan.gates)
    gate = scan.gates(i_elev);
    gate.data = z_to_refl(idb(gate.data)) / 11;
    gate.unit = "birds km^-3";
    scan.gates(i_elev) = gate;
end
end
