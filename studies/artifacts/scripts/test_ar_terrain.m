function tests = test_ar_terrain
%TEST_AR_TERRAIN Checks of the terrain-following steps and prior (dev).
%   Needs the ar_appalachians case (AR_A_cases.m). Run with runtests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_appalachians", [-89 -78 33.5 40], times);
load(fullfile(c.data_dir, "case.mat"), "g", "radar_scans", "observations")
column_mask = bmm4d.grid.build_domain_mask(g, c.domain);
active_mask_3d = bmm4d.grid.build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
[has_active, lowest] = max(active_mask_3d, [], 3);
lowest(~has_active) = g.sz(3) + 1;
testCase.TestData = struct("c", c, "g", g, "radar_scans", {radar_scans}, "observations", observations, ...
    "active_mask_3d", active_mask_3d, "lowest", lowest, "dt", seconds(diff(c.times)));
end

function testAslStepsEqualToolbox(testCase)
% With heights above sea level the steps are the toolbox steps
d = testCase.TestData;
inv = d.c.inversion;
reference = bmm4d.advect.build_radar_advection_steps(d.g, d.radar_scans, d.dt, ...
    idw_neighbors=inv.velocity_idw_neighbors, idw_power=inv.velocity_idw_power, boundary_mode="open", ...
    shoreline_crossing_coeff=inv.shoreline_crossing_coeff, active_mask_3d=d.active_mask_3d);
steps = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, inv, "asl");
for k = 1:numel(steps)
    expected = sortrows([reference{k}.source_id_valid reference{k}.dest_id_valid]);
    actual = sortrows([steps{k}.source_id_valid steps{k}.dest_id_valid]);
    verifyEqual(testCase, actual, expected)
end
end

function testAglLinksKeepHeightAboveGround(testCase)
% Every AGL link joins cells at the same height above their lowest active
% cell, so a field that depends only on that height is carried exactly;
% the ASL links are not
d = testCase.TestData;
steps_agl = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "agl");
steps_asl = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "asl");
[x, y, z] = ndgrid(1:d.g.sz(1), 1:d.g.sz(2), 1:d.g.sz(3));
height = z - d.lowest(sub2ind(d.g.sz(1:2), x, y));
field = exp(-max(height, 0) / 8);                     % dense near the ground, any column
residual = @(s) field(double(s.dest_id_valid)) - field(double(s.source_id_valid));
verifyEqual(testCase, max(abs(residual(steps_agl{1}))), 0)
verifyGreaterThan(testCase, max(abs(residual(steps_asl{1}))), 0.1)
% Same number of links where the ground is flat; more links over terrain
verifyGreaterThan(testCase, numel(steps_agl{1}.source_id_valid), 0.95 * numel(steps_asl{1}.source_id_valid))
end

function testAglPriorKeepsColumnOverRidges(testCase)
% The AGL prior's column over high ground is close to the valleys'; the
% production (ASL) prior loses the low layer
d = testCase.TestData;
inv = d.c.inversion;
active_id = find(d.active_mask_3d(:));
[prior_asl, info] = bmm4d.prior.build_density_prior_fast(d.g, d.radar_scans, d.observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id);
prior_agl = ar_profile_prior(d.g, d.radar_scans, info.radar_profiles, d.active_mask_3d, inv.prior_length_scale_m, "agl");
% By sea-level height and at the production length scale it is the production prior
prior_check = ar_profile_prior(d.g, d.radar_scans, info.radar_profiles, d.active_mask_3d, inv.prior_length_scale_m, "asl");
verifyEqual(testCase, prior_check, prior_asl, AbsTol=1e-4 * max(prior_asl, [], "all"))
dem = max(d.g.dem, 0);
ridge = any(d.active_mask_3d, 3) & dem > 1000;
valley = any(d.active_mask_3d, 3) & dem < 450;
ratio = @(prior) ridge_valley(d.g, active_id, prior(:, 4), ridge, valley);
fprintf("ridge / valley column, ASL prior %.2f, AGL prior %.2f\n", ratio(prior_asl), ratio(prior_agl))
verifyGreaterThan(testCase, ratio(prior_agl), 0.8)
verifyLessThan(testCase, ratio(prior_asl), ratio(prior_agl))
end

function testSigmaLinksCarryColumnCount(testCase)
% A sigma link joins cells at the same relative height of their columns
% and carries the ratio of the columns' air thickness: a field with the
% same birds in every column (density 1 / air thickness) is carried
% exactly, with no link lost into the ground or out of the top
d = testCase.TestData;
active_id = find(d.active_mask_3d(:));
air = zeros(d.g.sz);
air(active_id) = bmm4d.grid.air_fraction(d.g, active_id);
thickness = sum(air, 3);
field = zeros(d.g.sz);
field(d.active_mask_3d) = 1 ./ thickness(mod(active_id - 1, prod(d.g.sz(1:2))) + 1);
steps_sigma = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "sigma");
steps_agl = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "agl");
s = steps_sigma{1};
residual = field(double(s.dest_id_valid)) - double(s.ratio_valid) .* field(double(s.source_id_valid));
verifyLessThan(testCase, max(abs(residual)) / max(field, [], "all"), 1e-5)
% Every source cell with an active destination column keeps a link
linked = @(s) numel(unique(s.source_id_valid));
fprintf("source cells with a link: AGL %d, sigma %d\n", linked(steps_agl{1}), linked(s))
verifyGreaterThanOrEqual(testCase, linked(s), linked(steps_agl{1}))
end

function testSigmaPriorKeepsColumnOverRidges(testCase)
% The sigma prior puts the radars' column totals over the ridges too
d = testCase.TestData;
inv = d.c.inversion;
active_id = find(d.active_mask_3d(:));
[~, info] = bmm4d.prior.build_density_prior_fast(d.g, d.radar_scans, d.observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id);
prior_sigma = ar_profile_prior(d.g, d.radar_scans, info.radar_profiles, d.active_mask_3d, inv.prior_length_scale_m, "sigma");
verifyTrue(testCase, all(isfinite(prior_sigma), "all"))
dem = max(d.g.dem, 0);
ridge = any(d.active_mask_3d, 3) & dem > 1000;
valley = any(d.active_mask_3d, 3) & dem < 450;
ratio = ridge_valley(d.g, active_id, prior_sigma(:, 4), ridge, valley);
fprintf("ridge / valley column, sigma prior %.2f\n", ratio)
verifyGreaterThan(testCase, ratio, 0.85)
end

function testSmoothInfiniteLengthIsAsl(testCase)
% Height above a ground smoothed over an infinite length is height above
% sea level: the "smooth" steps and prior equal the "asl" ones
d = testCase.TestData;
inv = d.c.inversion;
steps_asl = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, inv, "asl");
steps_inf = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, inv, "smooth", smooth_length_m=inf);
for k = 1:numel(steps_asl)
    verifyEqual(testCase, sortrows([steps_inf{k}.source_id_valid steps_inf{k}.dest_id_valid]), ...
        sortrows([steps_asl{k}.source_id_valid steps_asl{k}.dest_id_valid]))
    verifyTrue(testCase, all(steps_inf{k}.coeff_valid == 1))
end
active_id = find(d.active_mask_3d(:));
[~, info] = bmm4d.prior.build_density_prior_fast(d.g, d.radar_scans, d.observations, ...
    "lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, "maxit", inv.maxit, ...
    "interpolation", inv.prior_interpolation, "length_scale_m", inv.prior_length_scale_m, ...
    "include_water", inv.prior_include_water, "active_id", active_id);
prior_asl = ar_profile_prior(d.g, d.radar_scans, info.radar_profiles, d.active_mask_3d, inv.prior_length_scale_m, "asl");
prior_inf = ar_profile_prior(d.g, d.radar_scans, info.radar_profiles, d.active_mask_3d, inv.prior_length_scale_m, ...
    "smooth", smooth_length_m=inf);
verifyEqual(testCase, prior_inf, prior_asl, AbsTol=1e-5 * max(prior_asl, [], "all"))
end

function testDetourLinksStayAtTheirLayer(testCase)
% Detour links join a source to an air cell of its own layer (ASL links
% keep the layer), with a coefficient of at most 1
d = testCase.TestData;
steps = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "asl", detour=true);
s = steps{1};
detour = s.detour_valid;
verifyGreaterThan(testCase, nnz(detour), 0)
[~, ~, source_layer] = ind2sub(d.g.sz, double(s.source_id_valid(detour)));
[~, ~, dest_layer] = ind2sub(d.g.sz, double(s.dest_id_valid(detour)));
verifyEqual(testCase, dest_layer, source_layer)
verifyTrue(testCase, all(d.active_mask_3d(s.dest_id_valid(detour))))
verifyLessThanOrEqual(testCase, max(double(s.coeff_valid(detour))), 1)
end

function testHybridLiftsOnlyNearTheGround(testCase)
% Hybrid links: a source at least H above its ground keeps its height
% above sea level (within the one-layer overlap of its mapped interval);
% below H it is lifted by at most the ground rise
d = testCase.TestData;
H = 1000;
steps = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "hybrid", fade_height_m=H);
s = steps{1};
z = squeeze(d.g.z);
ground = max(d.g.dem, 0);
nxy = prod(d.g.sz(1:2));
[~, ~, source_layer] = ind2sub(d.g.sz, double(s.source_id_valid));
[~, ~, dest_layer] = ind2sub(d.g.sz, double(s.dest_id_valid));
source_ground = ground(mod(double(s.source_id_valid) - 1, nxy) + 1);
rise = ground(mod(double(s.dest_id_valid) - 1, nxy) + 1) - source_ground;
source_height = z(source_layer) - source_ground;
aloft = source_height >= H + d.g.dz;
verifyLessThanOrEqual(testCase, max(abs(dest_layer(aloft) - source_layer(aloft))), 1)
near = source_height < H & rise > 0;
verifyLessThanOrEqual(testCase, max(z(dest_layer(near)) - z(source_layer(near)) - rise(near)), d.g.dz)
end

function testClimbLinksGoToLowestAirCell(testCase)
% With climb, a link into the ground goes to the lowest active cell of its
% destination column, and fewer source cells are left without a link
d = testCase.TestData;
steps = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "smooth", smooth_length_m=50e3, climb=true);
plain = ar_steps(d.g, d.radar_scans, d.dt, d.active_mask_3d, d.c.inversion, "smooth", smooth_length_m=50e3);
s = steps{1};
verifyGreaterThan(testCase, nnz(s.climb_valid), 0)
[dx, dy, dz] = ind2sub(d.g.sz, double(s.dest_id_valid(s.climb_valid)));
verifyEqual(testCase, dz, d.lowest(sub2ind(d.g.sz(1:2), dx, dy)))
verifyGreaterThan(testCase, numel(unique(s.source_id_valid)), numel(unique(plain{1}.source_id_valid)))
end

function ratio = ridge_valley(g, active_id, values, ridge, valley)
frame = zeros(g.sz);
frame(active_id) = double(values) .* bmm4d.grid.air_fraction(g, active_id);
column = sum(frame, 3);
ratio = mean(column(ridge)) / mean(column(valley));
end
