classdef laplacianAndGridTest < matlab.unittest.TestCase
    %LAPLACIANANDGRIDTEST Properties of the smoothness operator and grid helpers.
    % Grids are deliberately non-square (nx ~= ny ~= nz) so transposes and
    % y-fastest ordering mistakes show up as failures.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "toolbox"), IncludingSubfolders=true));
            % compute_vp needs Regularization Tools (external, local_settings.json)
            regtools_dir = local_folder("regtools_dir", fullfile("external", "regtools"));
            testCase.assumeTrue(isfolder(regtools_dir), "Regularization Tools not found.")
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(regtools_dir, IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testLaplacianSymmetricWithZeroRowSums(testCase)
            g = laplacianAndGridTest.grid([4 3 5]);
            L = LaplacianMatrix(g, water=false);

            testCase.verifySize(L, [prod(g.sz) prod(g.sz)])
            testCase.verifyEqual(full(L - L'), zeros(prod(g.sz)), AbsTol=1e-14)
            testCase.verifyEqual(full(sum(L, 2)), zeros(prod(g.sz), 1), AbsTol=1e-14)
            testCase.verifyEqual(full(L * ones(prod(g.sz), 1)), zeros(prod(g.sz), 1), AbsTol=1e-14)
            % Positive semi-definite with a one-dimensional null space.
            eigenvalues = eig(full(L));
            testCase.verifyGreaterThan(min(eigenvalues), -1e-12)
            testCase.verifyEqual(nnz(abs(eigenvalues) < 1e-10), 1)
        end

        function testLaplacianUsesXFastestOrdering(testCase)
            g = laplacianAndGridTest.grid([5 3 4]);
            L = LaplacianMatrix(g, water=false);
            [x_id, y_id, z_id] = ndgrid(1:g.sz(1), 1:g.sz(2), 1:g.sz(3));

            % A field quadratic in x only has a constant second difference
            % along x and must be invariant along y and z.
            response_x = reshape(L * x_id(:) .^ 2, g.sz);
            testCase.verifyEqual(response_x(2:end-1, :, :), -ones(g.sz - [2 0 0]) / 3, AbsTol=1e-12)
            testCase.verifyEqual(response_x - response_x(:, 1, 1), zeros(g.sz), AbsTol=1e-12)

            response_y = reshape(L * y_id(:) .^ 2, g.sz);
            testCase.verifyEqual(response_y(:, 2:end-1, :), -ones(g.sz - [0 2 0]) / 3, AbsTol=1e-12)
            testCase.verifyEqual(response_y - response_y(1, :, 1), zeros(g.sz), AbsTol=1e-12)

            response_z = reshape(L * z_id(:) .^ 2, g.sz);
            testCase.verifyEqual(response_z(:, :, 2:end-1), -ones(g.sz - [0 0 2]) / 3, AbsTol=1e-12)
            testCase.verifyEqual(response_z - response_z(1, 1, :), zeros(g.sz), AbsTol=1e-12)
        end

        function testCoeffHorverLimits(testCase)
            g = laplacianAndGridTest.grid([4 3 5]);
            [x_id, ~, z_id] = ndgrid(1:g.sz(1), 1:g.sz(2), 1:g.sz(3));

            % -1 keeps only vertical smoothing, 0.5 keeps only horizontal.
            L_vertical = LaplacianMatrix(g, water=false, coeff_horver=-1);
            L_horizontal = LaplacianMatrix(g, water=false, coeff_horver=0.5);
            testCase.verifyEqual(L_vertical * x_id(:), zeros(prod(g.sz), 1), AbsTol=1e-14)
            testCase.verifyEqual(L_horizontal * z_id(:) .^ 2, zeros(prod(g.sz), 1), AbsTol=1e-14)
            testCase.verifyGreaterThan(norm(L_vertical * z_id(:) .^ 2), 0)
            testCase.verifyGreaterThan(norm(L_horizontal * x_id(:) .^ 2), 0)
        end

        function testAllLandMaskMatchesPlainLaplacian(testCase)
            g = laplacianAndGridTest.grid([4 3 2]);
            g.mask_water = false(g.sz(1:2));

            testCase.verifyEqual(full(LaplacianMatrix(g, water=true)), ...
                full(LaplacianMatrix(g, water=false)), AbsTol=1e-14)
        end

        function testWaterMaskOnlyCutsShorelineLinks(testCase)
            % Only land-water pairs may be cut. With a straight coastline
            % (water along y = 1), along-shore land-land and water-water
            % links must be kept.
            g = laplacianAndGridTest.grid([4 3 1]);
            g.mask_water = false(g.sz(1:2));
            g.mask_water(:, 1) = true;

            L_plain = LaplacianMatrix(g, water=false);
            L_water = LaplacianMatrix(g, water=true);
            [row_id, col_id] = find(L_plain);
            is_offdiag = row_id ~= col_id;
            is_cross = g.mask_water(row_id) ~= g.mask_water(col_id);
            keep_id = sub2ind(size(L_plain), row_id(is_offdiag & ~is_cross), col_id(is_offdiag & ~is_cross));
            cut_id = sub2ind(size(L_plain), row_id(is_offdiag & is_cross), col_id(is_offdiag & is_cross));

            testCase.verifyEqual(full(L_water(cut_id)), zeros(numel(cut_id), 1))
            testCase.verifyEqual(full(sum(L_water, 2)), zeros(prod(g.sz), 1), AbsTol=1e-14)
            testCase.verifyEqual(full(L_water(keep_id)), full(L_plain(keep_id)), ...
                "Along-shore links between two land (or two water) cells were removed.", AbsTol=1e-14)
        end

        function testShorelineTransitionsUseXYOrientation(testCase)
            g = laplacianAndGridTest.grid([4 3 1]);
            g.mask_water = false(g.sz(1:2));
            g.mask_water(4, :) = true; % eastern column is water

            testCase.verifyEqual(get_water_mask(g), g.mask_water)
            source_id = sub2ind(g.sz(1:2), [3; 1], [2; 2]);
            dest_id = sub2ind(g.sz(1:2), [4; 1], [2; 3]);
            testCase.verifyEqual(find_land_water_transitions(g, source_id, dest_id), [true; false])
        end

        function testIdwWeightsAreConvexAndExactAtRadars(testCase)
            g = laplacianAndGridTest.grid([6 4 1]);
            radar_lonlat = [0 0; 5 3; 2 1; 2 1];
            radar_scans = cell(size(radar_lonlat, 1), 1);
            for i_r = 1:size(radar_lonlat, 1)
                radar_scans{i_r} = struct("lon", radar_lonlat(i_r, 1), "lat", radar_lonlat(i_r, 2));
            end

            [W, radar_xy] = build_horizontal_idw_weights(g, radar_scans, n_neighbors=4, power=2);

            testCase.verifySize(W, [prod(g.sz(1:2)) 4])
            testCase.verifyEqual(radar_xy, radar_lonlat .* 1000)
            testCase.verifyEqual(full(sum(W, 2)), ones(prod(g.sz(1:2)), 1), AbsTol=1e-12)
            testCase.verifyGreaterThanOrEqual(full(min(W, [], "all")), 0)
            % Cells holding a radar take that radar exactly; co-located radars split equally.
            testCase.verifyEqual(full(W(sub2ind(g.sz(1:2), 1, 1), :)), [1 0 0 0], AbsTol=1e-12)
            testCase.verifyEqual(full(W(sub2ind(g.sz(1:2), 6, 4), :)), [0 1 0 0], AbsTol=1e-12)
            testCase.verifyEqual(full(W(sub2ind(g.sz(1:2), 3, 2), :)), [0 0 0.5 0.5], AbsTol=1e-12)
            % A uniform per-radar value interpolates to a uniform field.
            testCase.verifyEqual(full(W * [7; 7; 7; 7]), 7 .* ones(prod(g.sz(1:2)), 1), AbsTol=1e-12)
        end

        function testComputeVpRecoversHorizontallyUniformProfile(testCase)
            rng(3)
            g = laplacianAndGridTest.grid([4 3 5]);
            vp_true = [5; 9; 12; 8; 3];
            field = repmat(reshape(vp_true, 1, 1, []), g.sz(1), g.sz(2), 1);
            n_obs = 40;
            F = sprand(n_obs, prod(g.sz), 0.3);
            obs_vector = sparse(F * field(:));

            vp = compute_vp(obs_vector, g, F, lambda=1e-8);

            testCase.verifyEqual(vp, vp_true, RelTol=1e-6)
        end
    end

    methods (Static, Access=private)
        function g = grid(sz)
            g.sz = sz;
            g.dx = 1000;
            g.dy = 1000;
            g.dz = 100;
            g.x = (0:sz(1) - 1) .* g.dx;
            g.y = (0:sz(2) - 1)' .* g.dy;
            g.z = (0:sz(3) - 1)' .* g.dz + g.dz / 2;
            g.f_X = @(lon) lon .* 1000;
            g.f_Y = @(lat) lat .* 1000;
        end
    end
end
