classdef advectionConventionTest < matlab.unittest.TestCase
    %ADVECTIONCONVENTIONTEST Direction, boundary and frame conventions of advection.
    % Uses non-square grids so x/y swaps cannot pass by symmetry.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "toolbox"), IncludingSubfolders=true));
            % Synthetic tools are not part of the package
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "synthetic")));
        end
    end

    methods (Test)
        function testEastwardProfileMovesBlobAlongX(testCase)
            g = advectionConventionTest.grid([5 4 2]);
            reference = zeros(g.sz);
            reference(2, 3, 1) = 1;
            reference(2, 3, 2) = 2;
            % One cell per step in layer 1, still air in layer 2.
            steps = build_profile_advection_steps(g, {[1 0; 0 0]}, 1000);

            fields = advect_density_sequence(reference, steps, [0 1]);

            expected = zeros(g.sz);
            expected(3, 3, 1) = 1;
            expected(2, 3, 2) = 2;
            testCase.verifyEqual(fields(:, :, :, 2), expected)
        end

        function testNorthwardProfileMovesBlobAlongY(testCase)
            g = advectionConventionTest.grid([5 4 1]);
            reference = zeros(g.sz);
            reference(2, 2, 1) = 1;
            steps = build_profile_advection_steps(g, {[0 1]}, 1000); % 2 cells, dy = 500 m

            fields = advect_density_sequence(reference, steps, [0 1]);

            expected = zeros(g.sz);
            expected(2, 4, 1) = 1;
            testCase.verifyEqual(fields(:, :, :, 2), expected)
        end

        function testBackwardStepUndoesForwardStep(testCase)
            g = advectionConventionTest.grid([6 5 2]);
            reference = zeros(g.sz);
            reference(3, 3, :) = [4 5];
            steps = build_profile_advection_steps(g, {[1 -1; -1 1], [1 -1; -1 1]}, 1000);

            fields = advect_density_sequence(reference, steps, [-1 0 1]);

            % dx = 1000 m and dy = 500 m: layer 1 moves (+1, -2) cells per
            % step, layer 2 moves (-1, +2). t = -1 holds the upstream position.
            testCase.verifyEqual(fields(2, 5, 1, 1), 4)
            testCase.verifyEqual(fields(4, 1, 2, 1), 5)
            testCase.verifyEqual(fields(4, 1, 1, 3), 4)
            testCase.verifyEqual(fields(2, 5, 2, 3), 5)
            testCase.verifyEqual(squeeze(sum(fields, [1 2 3])), [9; 9; 9])
        end

        function testFullFieldVelocityMatchesProfileVelocity(testCase)
            g = advectionConventionTest.grid([5 4 3]);
            uv_profile = [1 0; 0 -1; 2 1];
            uv_field = zeros([g.sz 2]);
            for i_z = 1:g.sz(3)
                uv_field(:, :, i_z, 1) = uv_profile(i_z, 1);
                uv_field(:, :, i_z, 2) = uv_profile(i_z, 2);
            end

            profile_steps = build_profile_advection_steps(g, {uv_profile}, 1000);
            field_steps = build_profile_advection_steps(g, uv_field, 1000);

            testCase.verifyEqual(field_steps{1}.source_id_valid, profile_steps{1}.source_id_valid)
            testCase.verifyEqual(field_steps{1}.dest_id_valid, profile_steps{1}.dest_id_valid)
        end

        function testOpenBoundaryDropsOutgoingCells(testCase)
            g = advectionConventionTest.grid([4 3 1]);
            steps = build_profile_advection_steps(g, {[1 0]}, 1000, boundary_mode="open");
            [source_x, ~] = ind2sub(g.sz(1:2), double(steps{1}.source_id_valid));
            [dest_x, ~] = ind2sub(g.sz(1:2), double(steps{1}.dest_id_valid));

            testCase.verifyEqual(numel(steps{1}.source_id_valid), (g.sz(1) - 1) * g.sz(2))
            testCase.verifyFalse(any(source_x == g.sz(1)))
            testCase.verifyEqual(dest_x, source_x + 1)
        end

        function testZeroGradientSequenceIsMirrorSymmetric(testCase)
            % Eastward advection of f must equal the mirror image of
            % westward advection of the mirrored f. With zero_gradient two
            % sources share the edge destination; advect_density_sequence
            % assigns with "last index wins", which makes east and west
            % flow behave differently at the outflow edge.
            g = advectionConventionTest.grid([3 1 1]);
            field = reshape([1; 2; 3], g.sz);
            east = build_profile_advection_steps(g, {[1 0]}, 1000, boundary_mode="zero_gradient");
            west = build_profile_advection_steps(g, {[-1 0]}, 1000, boundary_mode="zero_gradient");

            east_fields = advect_density_sequence(field, east, [0 1]);
            west_fields = advect_density_sequence(flip(field, 1), west, [0 1]);

            testCase.verifyEqual(east_fields(:, :, :, 2), flip(west_fields(:, :, :, 2), 1), ...
                "zero_gradient outflow edge depends on flow direction (duplicate destinations overwrite).")
        end

        function testSampleFrameMatchesGridAdvection(testCase)
            % A feature at x0 in the reference scan is carried downstream by
            % advect_density_sequence; advect_sample_coordinates must map
            % its later position back to x0 (and earlier positions too).
            g = advectionConventionTest.grid([8 3 1]);
            dt_seconds = [1000 2000 1000];
            uv = repmat(1 + 0.5i, 1, 4);
            reference = zeros(g.sz);
            reference(3, 1, 1) = 1;
            u_profiles = arrayfun(@(i) [real(uv(1, i)) 0], 1:3, UniformOutput=false);
            steps = build_profile_advection_steps(g, u_profiles(2:3), dt_seconds(2:3));
            fields = advect_density_sequence(reference, [build_profile_advection_steps(g, u_profiles(1), dt_seconds(1)); steps], [-1 0 1 2]);
            x_feature = zeros(1, 4);
            for i_t = 1:4
                [x_id, ~] = ind2sub(g.sz(1:2), find(fields(:, :, 1, i_t)));
                x_feature(i_t) = g.x(x_id);
            end

            [x_adv, y_adv] = advect_sample_coordinates( ...
                x_feature, zeros(1, 4), g.z(1) .* ones(1, 4), [0; 1000], ...
                [uv; uv], dt_seconds ./ 86400, 2);

            testCase.verifyEqual(x_adv, g.x(3) .* ones(1, 4), AbsTol=1e-9)
            testCase.verifyEqual(y_adv, [0.5 0 -1 -1.5] .* 1000, AbsTol=1e-9)
        end
    end

    methods (Static, Access=private)
        function g = grid(sz)
            g.sz = sz;
            g.dx = 1000;
            g.dy = 500;
            g.dz = 100;
            g.x = (0:sz(1) - 1) .* g.dx;
            g.y = (0:sz(2) - 1)' .* g.dy;
            g.z = (0:sz(3) - 1)' .* g.dz + g.dz / 2;
            g.mask_water = false(sz(1:2));
        end
    end
end
