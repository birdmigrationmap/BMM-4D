function metrics = compare_cleaning_methods(scans, method_names, z_edges)
%COMPARE_CLEANING_METHODS Summarize how cleaning methods differ on one volume.
%   SCANS is a cell of read_rd outputs for the same volume and the same
%   options except removeRain; SCANS{1} must be removeRain="none" and is the
%   reference. All outputs share the same sweep grids, so gate masks can be
%   compared directly. Removed fractions are weighted by the reference
%   density, so they measure how much of the signal each method discards.
%   The vertically integrated density (VID, birds km^-2) is computed from
%   the mean density profile in Z_EDGES (m).

n_method = numel(scans);
reference = scans{1};
n_sweep = numel(reference.dens.sweeps);

weight_total = 0;
weight_removed = zeros(1, n_method);
removed_masks = cell(n_sweep, n_method);
for i_elev = 1:n_sweep
    density_ref = reference.dens.sweeps(i_elev).data;
    valid_ref = ~isnan(density_ref);
    weight_total = weight_total + sum(density_ref(valid_ref));
    for i_m = 1:n_method
        removed = valid_ref & isnan(scans{i_m}.dens.sweeps(i_elev).data);
        removed_masks{i_elev, i_m} = removed;
        weight_removed(i_m) = weight_removed(i_m) + sum(density_ref(removed));
    end
end

% Union and intersection of the masks of every non-reference method; their
% difference is the signal on which the methods disagree.
weight_union = 0;
weight_intersection = 0;
for i_elev = 1:n_sweep
    density_ref = reference.dens.sweeps(i_elev).data;
    method_removed = cat(3, removed_masks{i_elev, 2:end});
    weight_union = weight_union + sum(density_ref(any(method_removed, 3)));
    weight_intersection = weight_intersection + sum(density_ref(all(method_removed, 3)));
end

metrics = struct();
metrics.n_gate_reference = sum(cellfun(@(d) nnz(~isnan(d)), {reference.dens.sweeps.data}));
for i_m = 1:n_method
    profile = scan_density_profile(scans{i_m}, z_edges);
    profile(isnan(profile)) = 0;
    metrics.("removed_" + method_names(i_m)) = weight_removed(i_m) / max(weight_total, eps);
    metrics.("vid_" + method_names(i_m)) = sum(profile .* diff(z_edges(:))) / 1000;
end
metrics.removed_union = weight_union / max(weight_total, eps);
metrics.disagreement = (weight_union - weight_intersection) / max(weight_total, eps);
end
