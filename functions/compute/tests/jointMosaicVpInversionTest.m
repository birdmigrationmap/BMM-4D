classdef jointMosaicVpInversionTest < matlab.unittest.TestCase
    %JOINTMOSAICVPINVERSIONTEST Tests compact observations and matrix-free inversion.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "functions"), IncludingSubfolders=true));
        end
    end

    methods (TestMethodSetup)
        function resetRandomSeed(testCase)
            original_rng = rng;
            testCase.addTeardown(@() rng(original_rng));
            rng(2, "twister")
        end
    end

    methods (Test)
        function testPrepareRadarObservationsCompactsRows(testCase)
            F = sparse([1 2 3 4], [1 2 1 2], 1, 4, 2);
            obs_matrix = sparse([1 3 1], [1 1 2], [2 3 4], 4, 2);
            info.obs_grid_size = [1 1 4 1];

            observations = prepare_radar_observations(F, obs_matrix, info);

            testCase.verifyEqual(observations.row_ids, uint32([1; 3]))
            testCase.verifyEqual(observations.F, F([1 3], :))
            testCase.verifyEqual(full(observations.values), [2 4; 3 0], AbsTol=0)
            testCase.verifyEqual(observations.obs_grid_size, info.obs_grid_size)
        end

        function testScalarPriorExpandsToFullSequence(testCase)
            g = jointMosaicVpInversionTest.smallGrid();
            info.obs_grid_size = [0 1 1 1];
            observations = prepare_radar_observations( ...
                sparse(0, prod(g.sz)), sparse(0, 2), info);
            steps = build_profile_advection_steps(g, {zeros(g.sz(3), 2)}, 1000);
            system = prepare_scans_advection(observations, g, steps, smoothness=false);

            fields = compute_scans_advection( ...
                system, 4, advection=0, prior=1, smoothness=0, maxit=20, tol=1e-12);

            testCase.verifyEqual(fields, 4 .* ones([g.sz 2]), AbsTol=1e-10)
        end

        function testMatrixFreeSolverMatchesExplicitObjective(testCase)
            [g, observations, steps, prior_field, A_explicit, b_explicit] = ...
                jointMosaicVpInversionTest.explicitReferenceSystem();
            system = prepare_scans_advection(observations, g, steps, water=false, coeff_horver=0);

            fields = compute_scans_advection( ...
                system, prior_field, advection=0.9, prior=0.4, smoothness=0.3, maxit=500, tol=1e-12);
            x_explicit = lsqr(A_explicit, b_explicit, 1e-12, 500);

            testCase.verifyLessThanOrEqual( ...
                norm(A_explicit * fields(:) - b_explicit), ...
                norm(A_explicit * x_explicit - b_explicit) + 1e-7)
            testCase.verifyEqual(fields(:), x_explicit, AbsTol=1e-6)
        end

        function testIrregularObservationMasksRun(testCase)
            g = jointMosaicVpInversionTest.smallGrid();
            F = sparse(rand(5, prod(g.sz)));
            obs_matrix = sparse([1 3 2 5], [1 1 2 2], [2 4 3 5], 5, 2);
            info.obs_grid_size = [1 1 5 1];
            observations = prepare_radar_observations(F, obs_matrix, info);
            steps = build_profile_advection_steps(g, {zeros(g.sz(3), 2)}, 1000);
            system = prepare_scans_advection(observations, g, steps, water=false);

            [fields, inversion_info] = compute_scans_advection( ...
                system, 1, prior=0.3, smoothness=0.2, advection=0.5, maxit=100, tol=1e-10);

            testCase.verifySize(fields, [g.sz 2])
            testCase.verifyTrue(isfinite(inversion_info.objective.total))
        end

        function testProjectedFistaEnforcesNonnegativeSolution(testCase)
            g.x = [0 1];
            g.y = 0;
            g.z = 0;
            g.sz = [2 1 1];
            g.dx = 1;
            g.dy = 1;
            g.dz = 1;
            F = speye(prod(g.sz));
            obs_matrix = sparse([-1; 2]);
            info.obs_grid_size = [1 1 2 1];
            observations = prepare_radar_observations(F, obs_matrix, info);
            system = prepare_scans_advection(observations, g, cell(0, 1), smoothness=false);

            [fields, inversion_info] = compute_scans_advection( ...
                system, 0, solver="projected_fista", prior=0, smoothness=0, advection=0, ...
                projected_step=1, projected_maxit=20, projected_tol=1e-12);

            testCase.verifyGreaterThanOrEqual(min(fields(:)), 0)
            testCase.verifyEqual(fields(:), [0; 2], AbsTol=1e-10)
            testCase.verifyEqual(inversion_info.solver.name, "matrix_free_projected_fista")
        end

        function testTemporalDensityPriorReducesProfileJumps(testCase)
            [g, radar_scans, observations] = jointMosaicVpInversionTest.temporalPriorInputs();

            [~, unsmoothed] = build_radar_density_prior( ...
                g, radar_scans, observations, lambda_profile=0, lambda_time_profile=0, maxit=100, tol=1e-12);
            [~, smoothed] = build_radar_density_prior( ...
                g, radar_scans, observations, lambda_profile=0, lambda_time_profile=10, maxit=100, tol=1e-12);
            unsmoothed_jump = sum(abs(diff(squeeze(unsmoothed.radar_profiles(:, 1, :)), 1, 2)), "all");
            smoothed_jump = sum(abs(diff(squeeze(smoothed.radar_profiles(:, 1, :)), 1, 2)), "all");

            testCase.verifyLessThan(smoothed_jump, unsmoothed_jump)
        end

        function testSplitLandWaterDensityPriorFitsSeparateProfiles(testCase)
            g.sz = [2 1 2];
            g.x = [0 1000];
            g.y = 0;
            g.z = [500 1500];
            g.mask_water = logical([0; 1]);
            g.f_X = @(lon) lon .* 1000;
            g.f_Y = @(lat) lat .* 1000;
            radar_scans = {struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1))};
            F = speye(prod(g.sz));
            obs_matrix = sparse([1; 2; 10; 20]);
            info.obs_grid_size = [1 1 1 4];
            observations = prepare_radar_observations(F, obs_matrix, info);

            [density_prior, density_prior_info] = build_radar_density_prior( ...
                g, radar_scans, observations, ...
                lambda_profile=0, lambda_time_profile=0, ...
                split_land_water_profiles=true, maxit=50, tol=1e-12);

            testCase.verifyEqual(squeeze(density_prior), [1 10; 2 20], AbsTol=1e-10)
            testCase.verifyEqual(squeeze(density_prior_info.radar_profiles(:, 1, 1, :)), ...
                [1 2; 10 20], AbsTol=1e-10)
            testCase.verifyEqual(density_prior_info.setup.profile_labels, ["land" "water"])
        end

        function testWaterConfigurationChangesLaplacian(testCase)
            g = jointMosaicVpInversionTest.smallGrid();
            g.mask_water = logical([0 1; 0 1]);
            info.obs_grid_size = [0 1 1 1];
            observations = prepare_radar_observations(sparse(0, prod(g.sz)), sparse(0, 1), info);

            land_water = prepare_scans_advection(observations, g, cell(0, 1), water=true);
            continuous = prepare_scans_advection(observations, g, cell(0, 1), water=false);

            testCase.verifyGreaterThan(norm(land_water.laplacian - continuous.laplacian, "fro"), 0)
        end

        function testShorelineCrossingCoeffUsesWaterMask(testCase)
            g = jointMosaicVpInversionTest.smallGrid();
            g.sz = [3 1 1];
            g.x = [0 1000 2000];
            g.y = 0;
            g.z = 500;
            g.mask_water = logical([0; 1; 1]);

            continuous = build_profile_advection_steps( ...
                g, [1 0], 1000, boundary_mode="open", shoreline_crossing_coeff=1);
            disconnected = build_profile_advection_steps( ...
                g, [1 0], 1000, boundary_mode="open", shoreline_crossing_coeff=0);
            downweighted = build_profile_advection_steps( ...
                g, [1 0], 1000, boundary_mode="open", shoreline_crossing_coeff=0.25);

            testCase.verifyEqual(double(continuous{1}.source_id_valid), [1; 2])
            testCase.verifyEqual(double(disconnected{1}.source_id_valid), 2)
            testCase.verifyEqual(downweighted{1}.coeff_valid, [0.25; 1], AbsTol=0)
        end

        function testAdvectionStepsUseCompactStorage(testCase)
            g = jointMosaicVpInversionTest.smallGrid();

            steps = build_profile_advection_steps(g, {[1 0; 0 0]}, 1000, boundary_mode="open");

            testCase.verifyClass(steps{1}.source_id_valid, "uint32")
            testCase.verifyClass(steps{1}.dest_id_valid, "uint32")
            testCase.verifyEqual(steps{1}.coeff_valid, 1, AbsTol=0)
        end

        function testRadarAdvectionBuilderMatchesFullVelocityPath(testCase)
            [g, radar_scans] = jointMosaicVpInversionTest.velocityInputs();
            velocity = build_radar_velocity_prior(g, radar_scans, idw_neighbors=2, idw_power=2);
            expected = build_profile_advection_steps(g, velocity, 1000, boundary_mode="open");

            actual = build_radar_advection_steps( ...
                g, radar_scans, 1000, idw_neighbors=2, idw_power=2, boundary_mode="open");

            testCase.verifyEqual(actual, expected)
        end

        function testRadarAdvectionBuilderUsesVaryingIntervals(testCase)
            [g, radar_scans] = jointMosaicVpInversionTest.varyingIntervalInputs();

            steps = build_radar_advection_steps(g, radar_scans, [600 1200]);

            testCase.verifyEqual([steps{1}.dt_seconds; steps{2}.dt_seconds], [600; 1200], AbsTol=0)
        end

        function testAdvectionStepsDropNonfiniteDestinations(testCase)
            g = jointMosaicVpInversionTest.smallGrid();

            steps = build_profile_advection_steps( ...
                g, cat(4, nan(g.sz), zeros(g.sz)), 1000, boundary_mode="open");

            testCase.verifyEmpty(steps{1}.source_id_valid)
            testCase.verifyEmpty(steps{1}.dest_id_valid)
        end

        function testRadarAdvectionBuilderDropsNaNVelocityCells(testCase)
            [g, radar_scans] = jointMosaicVpInversionTest.varyingIntervalInputs();
            radar_scans{1, 2}.uv(1, :) = [nan nan];

            steps = build_radar_advection_steps(g, radar_scans, [600 1200], ...
                idw_neighbors=1, idw_power=1, boundary_mode="open");

            testCase.verifyEqual(double(steps{1}.source_id_valid), (5:8)')
            testCase.verifyEqual(double(steps{1}.dest_id_valid), (5:8)')
            testCase.verifyNotEmpty(steps{2}.source_id_valid)
        end

        function testInvalidStepCountErrors(testCase)
            g = jointMosaicVpInversionTest.smallGrid();
            info.obs_grid_size = [0 1 1 1];
            observations = prepare_radar_observations(sparse(0, prod(g.sz)), sparse(0, 2), info);

            testCase.verifyError( ...
                @() prepare_scans_advection(observations, g, cell(0, 1)), ...
                "prepare_scans_advection:InvalidStepCount")
        end

        function testPipelineSaveContracts(testCase)
            simple_case_script = fileread(fullfile(testCase.ProjectRoot, ...
                "06_simple_mosaic", "scripts", "S06_A_generate_case.m"));
            realistic_case_script = fileread(fullfile(testCase.ProjectRoot, ...
                "07_regional_synthetic", "scripts", "S07_A_generate_case.m"));
            real_case_script = fileread(fullfile(testCase.ProjectRoot, ...
                "11_midwest_mosaic", "scripts", "S11_A_prepare_case.m"));
            simple_script = fileread(fullfile(testCase.ProjectRoot, ...
                "06_simple_mosaic", "scripts", "S06_B_run_inversion.m"));
            realistic_script = fileread(fullfile(testCase.ProjectRoot, ...
                "07_regional_synthetic", "scripts", "S07_C_run_inversion.m"));
            real_script = fileread(fullfile(testCase.ProjectRoot, ...
                "11_midwest_mosaic", "scripts", "S11_B_run_inversion.m"));

            testCase.verifyNotEmpty(regexp(simple_case_script, ...
                '"observations", "truth_sequences", "g", "radar_scans"', "once"))
            testCase.verifyNotEmpty(regexp(realistic_case_script, ...
                '"observations", "truth_sequences", "g", "radar_scans"', "once"))
            testCase.verifyNotEmpty(regexp(real_case_script, ...
                '"observations", "g", "radar_scans", "radar_selected", "case_controls"', "once"))
            testCase.verifyNotEmpty(regexp(simple_script, ...
                '"mosaic_fields", "inversion_info", "inversion_controls", "source_case_file"', "once"))
            testCase.verifyNotEmpty(regexp(realistic_script, ...
                '"mosaic_fields", "inversion_info", "inversion_controls", "source_case_file"', "once"))
            testCase.verifyNotEmpty(regexp(real_script, ...
                '"mosaic_field", "inversion_info", "inversion_controls", "source_case_file"', "once"))
            testCase.verifyFalse(contains(simple_script, '"density_prior_info"'))
            testCase.verifyFalse(contains(realistic_script, '"density_prior_info"'))
            testCase.verifyFalse(contains(real_script, '"prior_field", "uv_field_interval"'))
            testCase.verifyFalse(contains(simple_case_script, '"mosaics"'))
            testCase.verifyFalse(contains(realistic_case_script, '"mosaics"'))
        end
    end

    methods (Static, Access=private)
        function g = smallGrid()
            g.sz = [2 2 2];
            g.x = [0 1000];
            g.y = [0 1000];
            g.z = [500 1500];
            g.dx = 1000;
            g.dy = 1000;
            g.mask_water = false(g.sz(1:2));
        end

        function [g, observations, steps, prior_field, A_explicit, b_explicit] = explicitReferenceSystem()
            g = jointMosaicVpInversionTest.smallGrid();
            nxyz = prod(g.sz);
            n_obs = 6;
            F = sparse(rand(n_obs, nxyz));
            obs_matrix = sparse(rand(n_obs, 2));
            info.obs_grid_size = [1 1 n_obs 1];
            observations = prepare_radar_observations(F, obs_matrix, info);
            prior_field = reshape(1:(2 * nxyz), [g.sz 2]) ./ 10;
            steps = build_profile_advection_steps(g, {[1 0; 0 0]}, 1000, boundary_mode="open");
            laplacian = LaplacianMatrix(g, water=false, coeff_horver=0);
            data_scale = 1 / sqrt(n_obs);
            data_norm = max(norm(data_scale .* F, "fro"), 1);
            prior_scale = 0.4 .* data_norm ./ sqrt(nxyz);
            smoothness_scale = 0.3 .* data_norm ./ norm(laplacian, "fro");
            local_A = [data_scale .* F; prior_scale .* speye(nxyz); smoothness_scale .* laplacian];
            local_b_1 = [data_scale .* full(obs_matrix(:, 1)); prior_scale .* reshape(prior_field(:, :, :, 1), [], 1); zeros(nxyz, 1)];
            local_b_2 = [data_scale .* full(obs_matrix(:, 2)); prior_scale .* reshape(prior_field(:, :, :, 2), [], 1); zeros(nxyz, 1)];
            step = steps{1};
            n_step = numel(step.source_id_valid);
            coeff = repmat(step.coeff_valid, n_step, 1);
            advection_operator = sparse( ...
                [(1:n_step)'; (1:n_step)'], ...
                [double(step.source_id_valid); nxyz + double(step.dest_id_valid)], ...
                [-coeff; coeff], n_step, 2 * nxyz);
            advection_scale = 0.9 .* data_norm ./ norm(advection_operator, "fro");
            A_explicit = [ ...
                local_A sparse(size(local_A, 1), nxyz); ...
                sparse(size(local_A, 1), nxyz) local_A; ...
                advection_scale .* advection_operator];
            b_explicit = [local_b_1; local_b_2; zeros(n_step, 1)];
        end

        function [g, radar_scans, observations] = temporalPriorInputs()
            g.sz = [1 1 2];
            g.x = 0;
            g.y = 0;
            g.z = [500 1500];
            g.mask_water = false(1, 1);
            g.f_X = @(lon) lon .* 1000;
            g.f_Y = @(lat) lat .* 1000;
            radar_scans = { ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1)) ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 10, 0)) ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 20, 0))};
            F = speye(2);
            obs_matrix = sparse([0 10 0; 0 10 0]);
            info.obs_grid_size = [1 1 1 2];
            observations = prepare_radar_observations(F, obs_matrix, info);
        end

        function [g, radar_scans] = velocityInputs()
            g = jointMosaicVpInversionTest.smallGrid();
            g.f_X = @(lon) lon .* 1000;
            g.f_Y = @(lat) lat .* 1000;
            radar_scans = cell(2, 3);
            radar_scans{1, 1} = struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1), "uv", [1 2; 3 4]);
            radar_scans{2, 1} = struct("lon", 1, "lat", 1, "datetime", datetime(2020, 1, 1), "uv", [2 3; 4 5]);
            radar_scans{1, 2} = struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 10, 0), "uv", [2 3; 4 5]);
            radar_scans{2, 2} = struct("lon", 1, "lat", 1, "datetime", datetime(2020, 1, 1, 0, 10, 0), "uv", [3 4; 5 6]);
            radar_scans{1, 3} = struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 20, 0), "uv", [3 4; 5 6]);
            radar_scans{2, 3} = struct("lon", 1, "lat", 1, "datetime", datetime(2020, 1, 1, 0, 20, 0), "uv", [4 5; 6 7]);
        end

        function [g, radar_scans] = varyingIntervalInputs()
            g = jointMosaicVpInversionTest.smallGrid();
            g.f_X = @(lon) lon .* 1000;
            g.f_Y = @(lat) lat .* 1000;
            radar_scans = { ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1), "uv", zeros(2, 2)) ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 10, 0), "uv", zeros(2, 2)) ...
                struct("lon", 0, "lat", 0, "datetime", datetime(2020, 1, 1, 0, 30, 0), "uv", zeros(2, 2))};
        end
    end
end
