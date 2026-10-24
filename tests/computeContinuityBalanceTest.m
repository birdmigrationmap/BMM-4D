classdef computeContinuityBalanceTest < matlab.unittest.TestCase
    %COMPUTECONTINUITYBALANCETEST Tests continuity balance estimates.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "toolbox"), IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testZeroVelocityUnchangedDensityGivesZeroResidual(testCase)
            g = computeContinuityBalanceTest.smallGrid([2 2 1], 2);
            density_field = 4 .* ones([g.sz 2]);
            velocity_field_interval = zeros([g.sz 2 1]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.signed_source_count_birds, zeros(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.signed_boundary_count_birds, zeros(2, 2, 1), AbsTol=1e-12)
        end

        function testUniformDensityIncreaseGivesTakeoff(testCase)
            g = computeContinuityBalanceTest.smallGrid([2 2 1], 2);
            density_field = cat(4, ones(g.sz), 2 .* ones(g.sz));
            velocity_field_interval = zeros([g.sz 2 1]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.signed_source_count_birds, ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.takeoff_count_birds, ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.landing_count_birds, zeros(2, 2, 1), AbsTol=1e-12)
        end

        function testUniformDensityDecreaseGivesLanding(testCase)
            g = computeContinuityBalanceTest.smallGrid([2 2 1], 2);
            density_field = cat(4, 2 .* ones(g.sz), ones(g.sz));
            velocity_field_interval = zeros([g.sz 2 1]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.signed_source_count_birds, -ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.takeoff_count_birds, zeros(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.landing_count_birds, ones(2, 2, 1), AbsTol=1e-12)
        end

        function testConstantEastwardVelocityOnlyCrossesBoundary(testCase)
            g = computeContinuityBalanceTest.smallGrid([3 1 1], 2);
            density_field = ones([g.sz 2]);
            velocity_field_interval = zeros([g.sz 2 1]);
            velocity_field_interval(:, :, :, 1, 1) = 1000 / 3600;

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.signed_source_count_birds, zeros(3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.signed_boundary_count_birds, reshape([-1; 0; 1], 3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.entering_count_birds, reshape([1; 0; 0], 3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(balance.leaving_count_birds, reshape([0; 0; 1], 3, 1, 1), AbsTol=1e-12)
        end

        function testOutputSizesMatchHorizontalIntervals(testCase)
            g = computeContinuityBalanceTest.smallGrid([2 3 2], 3);
            density_field = ones([g.sz 3]);
            velocity_field_interval = zeros([g.sz 2 2]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifySize(balance.signed_source_count_birds, [2 3 2])
            testCase.verifySize(balance.takeoff_count_birds, [2 3 2])
            testCase.verifySize(balance.landing_count_birds, [2 3 2])
            testCase.verifySize(balance.entering_count_birds, [2 3 2])
            testCase.verifySize(balance.leaving_count_birds, [2 3 2])
            testCase.verifySize(balance.column_density_birds_per_km2, [2 3 3])
            testCase.verifySize(balance.mass_birds, [2 3 3])
            testCase.verifySize(balance.net_outflow_count_birds, [2 3 2])
        end

        function testMassAndColumnDensityUseVolumeDensity(testCase)
            g = computeContinuityBalanceTest.smallGrid([1 1 2], 2);
            density_field = 2 .* ones([g.sz 2]);
            velocity_field_interval = zeros([g.sz 2 1]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.mass_birds, reshape([4 4], 1, 1, 2), AbsTol=1e-12)
            testCase.verifyEqual(balance.column_density_birds_per_km2, reshape([4 4], 1, 1, 2), AbsTol=1e-12)
        end

        function testUnitsDocumentVolumeDensityContract(testCase)
            g = computeContinuityBalanceTest.smallGrid([1 1 1], 2);
            density_field = ones([g.sz 2]);
            velocity_field_interval = zeros([g.sz 2 1]);

            balance = compute_continuity_balance(g, density_field, velocity_field_interval);

            testCase.verifyEqual(balance.units.density_field, "birds/km^3")
            testCase.verifyEqual( ...
                balance.units.signed_source_count_birds, ...
                "birds per horizontal grid cell per scan interval")
            testCase.verifyEqual(balance.units.signed_source_rate_birds_per_hour, ...
                "birds per horizontal grid cell per hour")
            testCase.verifyEqual(balance.units.column_density_birds_per_km2, "birds/km^2")
        end
    end

    methods (Static, Access=private)
        function g = smallGrid(grid_sz, n_time)
            g.sz = grid_sz;
            g.x = (0:(grid_sz(1) - 1))' .* 1000;
            g.y = (0:(grid_sz(2) - 1)) .* 1000;
            g.z = reshape((0:(grid_sz(3) - 1)) .* 1000 + 500, 1, 1, []);
            g.dx = 1000;
            g.dy = 1000;
            g.dz = 1000;
            g.lon = g.x ./ 1000;
            g.lat = g.y ./ 1000;
            g.time = datetime(2020, 1, 1) + hours(0:(n_time - 1));
            g.mask_water = false(grid_sz(1:2));
        end
    end
end
