function out = regression_case_prod(raw_files, times, kernel_dir, work_dir, c)
%REGRESSION_CASE_PROD Production path (P2-P4) on a small Key West case.
%   Cleans RAW_FILES with clean_volume (MistNet), builds the radar's rows
%   on the case grid, then the Gaussian prior, the advection steps and the
%   lean solver, with the case settings C (prod_case("keywest") stored in
%   the reference file). Cleaned volumes are written to WORK_DIR. Static
%   clutter needs a whole night, so it is not applied. OUT holds the
%   compared outputs.

%% Clean the volumes (P2) ------------------------------------------------------
volume_files = fullfile(work_dir, regexprep(raw_files, ".*/", "") + ".mat");
for k = 1:numel(raw_files)
    status = clean_volume(raw_files(k), c.read_options, c.vp, c.filter, volume_files(k));
    assert(status == "ok", "clean_volume failed on %s: %s", raw_files(k), status)
end
% Sweeps differ in size: one flat list of all sweeps of all volumes
out.sweep_data = {};
out.uv = cell(numel(raw_files), 1);
for k = 1:numel(raw_files)
    loaded = load(volume_files(k), "volume");
    out.sweep_data = [out.sweep_data; {loaded.volume.dens.sweeps.data}'];
    out.uv{k} = loaded.volume.uv;
end

%% Grid and rows (P3) ----------------------------------------------------------
n_time = numel(times);
kernel_files = dir(fullfile(kernel_dir, "F_*.mat"));
kernel_ref = load(fullfile(kernel_files(1).folder, kernel_files(1).name), "x", "y", "z");
site = load(volume_files(1), "volume");
radar_meta = arrayfun(@(t) struct("lon", site.volume.lon, "lat", site.volume.lat, "datetime", t), times, ...
    UniformOutput=false);
g = create_grid({struct("x", kernel_ref.x, "y", kernel_ref.y, "z", kernel_ref.z)}, radar_meta, ...
    c.grid_bounds(1), c.grid_bounds(2), c.grid_bounds(3), c.grid_bounds(4), c.heightmax);
[lon_grid, lat_grid] = ndgrid(g.lon, g.lat);
g.mask_water(inpolygon(lon_grid, lat_grid, c.water_override(1, :), c.water_override(2, :))) = true;
g = coarsen_grid(g, c.solve_aggregate);
rows = build_radar_rows(reshape(volume_files, 1, []), n_time, g, kernel_dir, c.min_row_mass, c.solve_aggregate, ...
    [], c.static_clutter.max_elev, c.near_range_m, c.no_echo_max_range_m);
radar_scans = cell(1, n_time);
for i_t = 1:n_time
    radar_scans{i_t} = struct("station", rows.site.station, "lat", rows.site.lat, "lon", rows.site.lon, ...
        "height", rows.site.height, "datetime", times(i_t), "uv", rows.uv(:, :, i_t));
end
observations = struct("F", rows.F, "values", rows.values, "row_radar_id", ones(size(rows.F, 1), 1, "uint16"), ...
    "row_range_m", rows.row_range_m, "row_elev_deg", rows.row_elev_deg);
out.grid_size = g.sz;
out.mask_water = g.mask_water;
out.F = observations.F;
out.values = observations.values;

%% Prior, advection, solve (P4, lean path) -----------------------------------
inv = c.inversion;
active_mask_3d = build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, []);
active_id = find(active_mask_3d(:));
out.active_id = active_id;
out.prior = single(build_density_prior_fast(g, radar_scans, observations, ...
    lambda_profile=inv.lambda_profile, lambda_time_profile=inv.lambda_time_profile, ...
    idw_neighbors=inv.density_prior_idw_neighbors, idw_power=inv.density_prior_idw_power, maxit=inv.maxit, ...
    interpolation=inv.prior_interpolation, length_scale_m=inv.prior_length_scale_m, ...
    include_water=inv.prior_include_water, active_id=active_id));
steps = build_radar_advection_steps(g, radar_scans, seconds(diff(times)), ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=active_mask_3d);
field = solve_advection_lean(observations, g, active_mask_3d, steps, out.prior, ...
    prior=inv.prior, smoothness=inv.smoothness, advection=inv.advection, coeff_horver=inv.coeff_horver, ...
    tol=inv.tol, maxit=inv.maxit, store_class="single", advection_relax_exponent=inv.advection_relax_exponent, ...
    nonnegative=inv.nonnegative, projected_maxit=inv.projected_maxit);
out.field = single(field);
end
