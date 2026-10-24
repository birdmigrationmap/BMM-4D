classdef forwardOperatorTest < matlab.unittest.TestCase
    %FORWARDOPERATORTEST Index conventions of compute_F_elev and construct_F.
    % compute_F_elev runs on a tiny non-square kernel grid inside a
    % temporary folder (it saves F_*.mat into the current folder), and the
    % saved ids are decoded with the same rule as load_F_elev.

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
        function testFastKernelPointsAlongAzimuth(testCase)
            kernel = forwardOperatorTest.computeKernel(testCase, "fast");
            forwardOperatorTest.verifyAzimuthDirections(testCase, kernel)
        end

        function testMemoryEfficientKernelPointsAlongAzimuth(testCase)
            kernel = forwardOperatorTest.computeKernel(testCase, "memory_efficient");
            forwardOperatorTest.verifyAzimuthDirections(testCase, kernel)
        end

        function testKernelModesAgree(testCase)
            fast = forwardOperatorTest.computeKernel(testCase, "fast");
            memory_efficient = forwardOperatorTest.computeKernel(testCase, "memory_efficient");

            fast_rows = sortrows([double(fast.Felevidc) double(fast.Felevfc)]);
            memory_rows = sortrows([double(memory_efficient.Felevidc) double(memory_efficient.Felevfc)]);
            testCase.verifyEqual(memory_rows(:, 1:3), fast_rows(:, 1:3), ...
                "memory_efficient mode stores different kernel cells than fast mode on a non-square grid.")
        end

        function testConstructFPlacesKernelAroundRadar(testCase)
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([0 0]);

            [F, obs_matrix, info] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);

            testCase.verifyEqual(info.obs_grid_size, [10 10 2 2])
            testCase.verifySize(F, [prod(info.obs_grid_size) prod(g.sz)])
            testCase.verifySize(obs_matrix, [prod(info.obs_grid_size) 1])
            % Gate (range 1, azimuth 1) of each radar sees the cell east of
            % the radar, gate (2, 1) the cell north, both in the lowest layer.
            radar_xy_id = [5 4; 2 3];
            for i_r = 1:2
                east_row = sub2ind(info.obs_grid_size, 1, 1, i_r, 1);
                north_row = sub2ind(info.obs_grid_size, 2, 1, i_r, 1);
                testCase.verifyEqual(find(F(east_row, :)), ...
                    sub2ind(g.sz, radar_xy_id(i_r, 1) + 1, radar_xy_id(i_r, 2), 1))
                testCase.verifyEqual(find(F(north_row, :)), ...
                    sub2ind(g.sz, radar_xy_id(i_r, 1), radar_xy_id(i_r, 2) + 1, 1))
            end
            testCase.verifyEqual(info.obs_block_offset, [0 200; 100 300])
            % The upper elevation uses its own row block and kernel layer.
            upper_row = sub2ind(info.obs_grid_size, 1, 1, 2, 2);
            testCase.verifyEqual(find(F(upper_row, :)), sub2ind(g.sz, 3, 3, 2))
        end

        function testConstructFWithSingleElevationKernel(testCase)
            % max(vertcat(sz_raz)) collapses to a scalar when only one
            % kernel is used, so obs_grid_size loses a dimension.
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([0 0]);
            elevation_kernels = elevation_kernels(1);
            for i_r = 1:numel(radar_scans)
                radar_scans{i_r}.kernel_idx = 1;
                radar_scans{i_r}.dens.sweeps = radar_scans{i_r}.dens.sweeps(1);
            end

            [~, ~, info] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);

            testCase.verifyEqual(info.obs_grid_size, [10 10 2 1])
        end

        function testObservationRowsRoundTripThroughFormatData(testCase)
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([0 0]);
            [F, obs_matrix, info] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);
            observations = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info);

            [obs_values, obs_id, azimuth_deg, elevation_deg, range_m] = ...
                bmm4d.forward.format_data(obs_matrix(:, 1), elevation_kernels, info);

            testCase.verifyEqual(obs_values, full(obs_matrix(obs_id, 1)))
            [i_range, i_az, i_r, i_f] = ind2sub(info.obs_grid_size, obs_id);
            expected = arrayfun(@(r, a, k, f) radar_scans{k}.dens.sweeps(f).data(r, a), i_range, i_az, i_r, i_f);
            testCase.verifyEqual(obs_values, expected + eps, AbsTol=1e-12)
            testCase.verifyEqual(range_m, reshape(elevation_kernels{1}.range(i_range), [], 1))
            testCase.verifyEqual(azimuth_deg, reshape(elevation_kernels{1}.az(i_az), [], 1))
            kernel_elevation = [0.5 1.5];
            testCase.verifyEqual(elevation_deg(:), reshape(kernel_elevation(i_f), [], 1))
            % Compacted rows keep the full-grid row ids of rows with data and kernel support.
            testCase.verifyEqual(observations.F, F(double(observations.row_ids), :))
            testCase.verifyEqual(sort(double(observations.row_ids)), ...
                find(any(obs_matrix, 2) & any(F, 2)))
        end

        function testCompactOperatorMatchesFullOperatorRows(testCase)
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([0 0]);
            % Leave part of the sweep unobserved so row compaction matters.
            for i_r = 1:numel(radar_scans)
                radar_scans{i_r}.dens.sweeps(2).data(1, 1) = nan;
            end

            [F_full, ~, info_full] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);
            [F_compact, ~, info_compact] = construct_F_compact(radar_scans, elevation_kernels, g);

            full_rows = find(any(F_full, 2));
            testCase.verifyEqual(size(F_compact, 1), numel(full_rows))
            testCase.verifyEqual(sortrows(full(F_compact)), sortrows(full(F_full(full_rows, :))))
            [~, ~, radar_id_full, kernel_id_full] = ind2sub(info_full.obs_grid_size, full_rows);
            testCase.verifyEqual(sortrows([double(info_compact.row_radar_id) double(info_compact.row_kernel_id)]), ...
                sortrows([radar_id_full kernel_id_full]))
        end

        function testRadarHeightShiftsKernelByRoundedLayers(testCase)
            % Kernel layers are relative to the antenna, so a radar at 3 m
            % must shift them by round(height / dz) = 0 layers, not one.
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([3 0]);

            [F, ~, info] = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);

            east_row = sub2ind(info.obs_grid_size, 1, 1, 1, 1);
            [~, ~, z_id] = ind2sub(g.sz, find(F(east_row, :)));
            testCase.verifyEqual(z_id, 1, ...
                "Radar at 3 m was placed in the wrong vertical layer.")
        end

        function testAggregateSumsFineCellsIntoBlocks(testCase)
            % F on a coarse grid of [2 2 2] blocks must equal the fine F
            % times the block membership matrix (row sums preserved).
            [g, elevation_kernels, radar_scans] = forwardOperatorTest.constructInputs([0 0]);
            aggregate = [2 2 2];
            F_fine = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g);
            F_coarse = bmm4d.forward.construct_F(radar_scans, elevation_kernels, g, aggregate=aggregate);

            coarse_size = ceil(g.sz ./ aggregate);
            [ix, iy, iz] = ind2sub(g.sz, (1:prod(g.sz))');
            block_id = sub2ind(coarse_size, ceil(ix / 2), ceil(iy / 2), ceil(iz / 2));
            membership = sparse((1:prod(g.sz))', block_id, 1, prod(g.sz), prod(coarse_size));
            testCase.verifyEqual(full(F_coarse), full(F_fine * membership), AbsTol=1e-12)
        end
    end

    methods (Static, Access=private)
        function kernel = computeKernel(testCase, mode)
            working_folder = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            x = -3000:500:3000;
            y = (-2000:500:2000)';
            z = reshape(50:100:250, 1, 1, []);
            elevation_kernel = struct( ...
                "name", "test", ...
                "elev", 1, ...
                "range", 1000:500:2500, ...
                "az", [0 90 180 270]);

            evalc("bmm4d.forward.compute_F_elev(x, y, z, 1, 1, 5, elevation_kernel, 0, mode)");

            files = dir(fullfile(working_folder.Folder, "F_test_*.mat"));
            testCase.assertNumElements(files, 1)
            kernel = load(fullfile(files.folder, files.name));
            delete(fullfile(files.folder, files.name))
            % Same decoding as load_F_elev.
            sz_coarse_grid = [numel(kernel.x) numel(kernel.y) numel(kernel.z)];
            [kernel.gy, kernel.gx, kernel.gz] = ind2sub(sz_coarse_grid([2 1 3]), double(kernel.Felevidc(:, 3)));
        end

        function verifyAzimuthDirections(testCase, kernel)
            testCase.assertNotEmpty(kernel.Felevidc)
            x_m = reshape(kernel.x(kernel.gx(:)), [], 1);
            y_m = reshape(kernel.y(kernel.gy(:)), [], 1);
            az_id = double(kernel.Felevidc(:, 2));
            weight = double(kernel.Felevfc);
            % Weighted mean position per azimuth: 0 = north (+y), 90 = east (+x).
            expected_sign = [0 1; 1 0; 0 -1; -1 0];
            for i_az = 1:4
                in_az = az_id == i_az;
                testCase.assertNotEmpty(find(in_az, 1))
                mean_xy = [sum(weight(in_az) .* x_m(in_az)) sum(weight(in_az) .* y_m(in_az))] ./ sum(weight(in_az));
                testCase.verifyEqual(sign(round(mean_xy)), expected_sign(i_az, :), ...
                    sprintf("Kernel mass for azimuth %d is on the wrong side.", kernel.elevation_kernel.az(i_az)))
            end
            % Heights should follow the beam: lowest layers only at 1 deg and <= 2.5 km.
            testCase.verifyLessThanOrEqual(max(kernel.gz), 2)
        end

        function [g, elevation_kernels, radar_scans] = constructInputs(radar_height)
            g.sz = [10 6 3];
            g.x = (0:9) .* 500;
            g.y = (0:5)' .* 500;
            g.z = (50:100:250)';
            g.dz = 100;
            g.lon = g.x ./ 1000;
            g.lat = g.y ./ 1000;

            % Kernel coarse grid 5 x 5 x 3 centred on the radar (3, 3).
            kernel.elev = 0.5;
            kernel.range = (1:10) .* 250;
            kernel.az = (0:9) .* 36;
            kernel.sz_raz = [10 10];
            kernel.sz_coarse_grid = [5 5 3];
            kernel.x = (-2:2) .* 500;
            kernel.y = (-2:2)' .* 500;
            kernel.z = reshape(50:100:250, 1, 1, []);
            kernel.F_raz_id = sub2ind(kernel.sz_raz, [1; 2], [1; 1]);
            kernel.F_gr_x = [4; 3];
            kernel.F_gr_y = [3; 4];
            kernel.F_gr_z = [1; 1];
            kernel.Fval = single([1; 1]);
            % Second elevation with the same geometry, one layer higher.
            kernel_upper = kernel;
            kernel_upper.elev = 1.5;
            kernel_upper.F_gr_z = [2; 2];
            elevation_kernels = {kernel kernel_upper};

            radar_lonlat = [2 1.5; 0.5 1];
            radar_scans = cell(2, 1);
            for i_r = 1:2
                sweep = struct("elev", {0.5 1.5}, "range", kernel.range, "az", kernel.az, ...
                    "data", {reshape(1:100, 10, 10) + 1000 * i_r, reshape(1:100, 10, 10) + 2000 * i_r});
                radar_scans{i_r} = struct( ...
                    "lon", radar_lonlat(i_r, 1), "lat", radar_lonlat(i_r, 2), ...
                    "height", radar_height(i_r), "station", "R" + i_r, ...
                    "datetime", datetime(2020, 1, 1), "kernel_idx", [1; 2], ...
                    "dens", struct("sweeps", {sweep}));
            end
        end
    end
end
