%% Benchmark joint mosaic VP inversion at the target radar/time scale

project_root = setup_paths();
rng(1)

%% Build compact 100-radar, 20-time benchmark inputs
n_radar = 100;
n_time = 20;
g.x = (0:99) .* 1000;
g.y = (0:99) .* 1000;
g.z = [500 1500];
g.sz = [numel(g.x) numel(g.y) numel(g.z)];
g.dx = 1000;
g.dy = 1000;
g.mask_water = false(g.sz(1:2));
g.f_X = @(lon) lon .* 1000;
g.f_Y = @(lat) lat .* 1000;

radar_x = linspace(g.x(1), g.x(end), n_radar)';
radar_y = mod((0:n_radar - 1)' .* 37, g.y(end));
radar_scans = cell(n_radar, n_time);
for i_t = 1:n_time
    for i_r = 1:n_radar
        radar_scans{i_r, i_t} = struct( ...
            "lon", radar_x(i_r) ./ 1000, ...
            "lat", radar_y(i_r) ./ 1000, ...
            "datetime", datetime(2020, 1, 1) + minutes(30 * (i_t - 1)));
    end
end

info.obs_grid_size = [1 1 n_radar 1];
info.obs_block_offset = (0:n_radar - 1)';
field_ids = sub2ind(g.sz, ...
    round(linspace(1, g.sz(1), n_radar))', ...
    mod((0:n_radar - 1)', g.sz(2)) + 1, ...
    ones(n_radar, 1));
F = sparse((1:n_radar)', field_ids, 1, n_radar, prod(g.sz));
obs_matrix = sparse(rand(n_radar, n_time));
observations = prepare_radar_observations(F, obs_matrix, info);

%% Measure preparation time and retained operator sizes
benchmark_tic = tic;
[prior_field, density_prior_info] = build_radar_density_prior( ...
    g, radar_scans, observations, ...
    lambda_profile=10, ...
    maxit=20);
preparation_seconds = toc(benchmark_tic);

operator_memory = whos("prior_field", "density_prior_info");
fprintf("100-radar, 20-time density-prior preparation: %.2f s, %.1f MB retained\n", ...
    preparation_seconds, sum([operator_memory.bytes]) ./ 1024^2)

%% Measure one low-iteration full joint solve
steps = build_profile_advection_steps( ...
    g, zeros(g.sz(3), 2, n_time - 1), 1800, ...
    boundary_mode="open");
system_tic = tic;
system = prepare_scans_advection(observations, g, steps, water=false, coeff_horver=0);
system_seconds = toc(system_tic);
system_memory = whos("system");
fprintf("Prepared matrix-free system: %.2f s, %.1f MB retained\n", ...
    system_seconds, system_memory.bytes ./ 1024^2)
solve_tic = tic;
[fields] = compute_scans_advection( ...
    system, prior_field, ...
    prior=0.3, ...
    smoothness=1, ...
    advection=3, ...
    maxit=2);
solve_seconds = toc(solve_tic);
solution_memory = whos("fields");
fprintf("One low-iteration joint solve: %.2f s, %.1f MB returned solutions\n", ...
    solve_seconds, sum([solution_memory.bytes]) ./ 1024^2)
