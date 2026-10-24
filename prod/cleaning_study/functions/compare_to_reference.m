function metrics = compare_to_reference(scans, method_names, reference_id, z_edges)
%COMPARE_TO_REFERENCE Rain masks of several methods against a reference method.
%   SCANS{1} is read_rd with removeRain="none" (all echo), SCANS{reference_id}
%   the reference (e.g. MistNet), the others candidate methods, all with the
%   same other options. Shares are weighted by the "none" density:
%   - removed_<m>: signal removed by method m;
%   - over_<m>: removed by m but kept by the reference (extra removal);
%   - missed_<m>: removed by the reference but kept by m (rain left in);
%   - vid_<m>: vertically integrated density (birds km^-2), mean profile in z_edges.

reference = scans{1};
weight_total = 0;
n_method = numel(scans);
removed = zeros(1, n_method);
over = zeros(1, n_method);
missed = zeros(1, n_method);
for i_elev = 1:numel(reference.dens.sweeps)
    density = reference.dens.sweeps(i_elev).data;
    valid = ~isnan(density) & density > 0;
    weight_total = weight_total + sum(density(valid));
    ref_removed = valid & isnan(scans{reference_id}.dens.sweeps(i_elev).data);
    for i_m = 1:n_method
        m_removed = valid & isnan(scans{i_m}.dens.sweeps(i_elev).data);
        removed(i_m) = removed(i_m) + sum(density(m_removed));
        over(i_m) = over(i_m) + sum(density(m_removed & ~ref_removed));
        missed(i_m) = missed(i_m) + sum(density(ref_removed & ~m_removed));
    end
end
metrics = struct();
for i_m = 1:n_method
    profile = scan_density_profile(scans{i_m}, z_edges);
    profile(isnan(profile)) = 0;
    metrics.("removed_" + method_names(i_m)) = removed(i_m) / max(weight_total, eps);
    metrics.("over_" + method_names(i_m)) = over(i_m) / max(weight_total, eps);
    metrics.("missed_" + method_names(i_m)) = missed(i_m) / max(weight_total, eps);
    metrics.("vid_" + method_names(i_m)) = sum(profile .* diff(z_edges(:))) / 1000;
end
end
