classdef sinksourceTest < matlab.unittest.TestCase
    %SINKSOURCETEST Tests continuity-based source and sink estimates.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "functions"), IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testZeroVelocityUnchangedDensityGivesZeroResidual(testCase)
            grid = sinksourceTest.smallGrid([2 2 1], 2);
            density_field = 4 .* ones([grid.sz 2]);
            velocity_field_interval = zeros([grid.sz 2 1]);

            flux = sinksource(grid, density_field, velocity_field_interval);

            testCase.verifyEqual(flux.source_birds, zeros(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.boundary_flux_birds, zeros(2, 2, 1), AbsTol=1e-12)
        end

        function testUniformDensityIncreaseGivesTakeoff(testCase)
            grid = sinksourceTest.smallGrid([2 2 1], 2);
            density_field = cat(4, ones(grid.sz), 2 .* ones(grid.sz));
            velocity_field_interval = zeros([grid.sz 2 1]);

            flux = sinksource(grid, density_field, velocity_field_interval);

            testCase.verifyEqual(flux.source_birds, ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.takingoff_birds, ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.landing_birds, zeros(2, 2, 1), AbsTol=1e-12)
        end

        function testUniformDensityDecreaseGivesLanding(testCase)
            grid = sinksourceTest.smallGrid([2 2 1], 2);
            density_field = cat(4, 2 .* ones(grid.sz), ones(grid.sz));
            velocity_field_interval = zeros([grid.sz 2 1]);

            flux = sinksource(grid, density_field, velocity_field_interval);

            testCase.verifyEqual(flux.source_birds, -ones(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.takingoff_birds, zeros(2, 2, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.landing_birds, -ones(2, 2, 1), AbsTol=1e-12)
        end

        function testConstantEastwardVelocityOnlyCrossesBoundary(testCase)
            grid = sinksourceTest.smallGrid([3 1 1], 2);
            density_field = ones([grid.sz 2]);
            velocity_field_interval = zeros([grid.sz 2 1]);
            velocity_field_interval(:, :, :, 1, 1) = 1000 / 3600;

            flux = sinksource(grid, density_field, velocity_field_interval);

            testCase.verifyEqual(flux.source_birds, zeros(3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.boundary_flux_birds, reshape([-1; 0; 1], 3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.entering_birds, reshape([-1; 0; 0], 3, 1, 1), AbsTol=1e-12)
            testCase.verifyEqual(flux.leaving_birds, reshape([0; 0; 1], 3, 1, 1), AbsTol=1e-12)
        end

        function testOutputSizesMatchHorizontalIntervals(testCase)
            grid = sinksourceTest.smallGrid([2 3 2], 3);
            density_field = ones([grid.sz 3]);
            velocity_field_interval = zeros([grid.sz 2 2]);

            flux = sinksource(grid, density_field, velocity_field_interval);

            testCase.verifySize(flux.source_birds, [2 3 2])
            testCase.verifySize(flux.takingoff_birds, [2 3 2])
            testCase.verifySize(flux.landing_birds, [2 3 2])
            testCase.verifySize(flux.entering_birds, [2 3 2])
            testCase.verifySize(flux.leaving_birds, [2 3 2])
            testCase.verifySize(flux.column_density_birds_per_km2, [2 3 3])
            testCase.verifySize(flux.mass_birds, [2 3 3])
            testCase.verifySize(flux.net_outflow_birds, [2 3 2])
        end
    end

    methods (Static, Access=private)
        function grid = smallGrid(grid_sz, n_time)
            grid.sz = grid_sz;
            grid.x = (0:(grid_sz(1) - 1))' .* 1000;
            grid.y = (0:(grid_sz(2) - 1)) .* 1000;
            grid.z = reshape((0:(grid_sz(3) - 1)) .* 1000 + 500, 1, 1, []);
            grid.dx = 1000;
            grid.dy = 1000;
            grid.dz = 1000;
            grid.lon = grid.x ./ 1000;
            grid.lat = grid.y ./ 1000;
            grid.time = datetime(2020, 1, 1) + hours(0:(n_time - 1));
            grid.mask_water = false(grid_sz(1:2));
        end
    end
end
