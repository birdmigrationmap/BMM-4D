classdef terrainFollowingTest < matlab.unittest.TestCase
    %TERRAINFOLLOWINGTEST Advection links over terrain (ground_smoothing_m).
    % Synthetic grid: ground rising eastward from 0 to 800 m with a 1500 m
    % ridge two cells wide across it, one radar, uniform wind of one cell
    % eastward per step in every layer.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            if isempty(mpmlist("bmm4d"))
                mpminstall(fullfile(testCase.ProjectRoot, "toolbox"), InPlace=true, Temporary=true, Prompt=false, Verbosity="quiet");
            end
        end
    end

    methods (Test)
        function testFlatGroundIsSeaLevel(testCase)
            % Without terrain the smoothed ground is flat: the links are the
            % sea-level links, all with coefficient 1
            [g, scans, mask] = terrainFollowingTest.synthetic_case(false);
            asl = bmm4d.advect.build_radar_advection_steps(g, scans, 1000, active_mask_3d=mask);
            smooth = bmm4d.advect.build_radar_advection_steps(g, scans, 1000, active_mask_3d=mask, ground_smoothing_m=50e3);
            testCase.verifyEqual(sortrows([smooth{1}.source_id_valid smooth{1}.dest_id_valid]), ...
                sortrows([asl{1}.source_id_valid asl{1}.dest_id_valid]))
            testCase.verifyEqual(smooth{1}.coeff_valid, ones(numel(smooth{1}.source_id_valid), 1))
        end

        function testLinksKeepHeightAboveSmoothedGround(testCase)
            % Away from the ground, the share-weighted height of a source's
            % destinations is its own height plus the rise of the smoothed
            % ground: a fraction of a layer is split, not rounded away
            [g, scans, mask] = terrainFollowingTest.synthetic_case(true);
            steps = bmm4d.advect.build_radar_advection_steps(g, scans, 1000, active_mask_3d=mask, ground_smoothing_m=50e3);
            s = steps{1};
            ground_smooth = bmm4d.grid.smooth_ground(g, 50e3);
            nxy = prod(g.sz(1:2));
            [~, ~, source_layer] = ind2sub(g.sz, double(s.source_id_valid));
            [~, ~, dest_layer] = ind2sub(g.sz, double(s.dest_id_valid));
            source_column = mod(double(s.source_id_valid) - 1, nxy) + 1;
            dest_column = mod(double(s.dest_id_valid) - 1, nxy) + 1;
            z = g.z(:);
            % Sources whose links all reach full air cells and cover the
            % whole source interval (shares add up to 1)
            full_air = z(dest_layer) - g.dz / 2 >= max(g.dem(dest_column), 0) & ~s.climb_valid;
            source_bad = accumarray(double(s.source_id_valid), double(~full_air), [prod(g.sz) 1]) > 0;
            source_share = accumarray(double(s.source_id_valid), double(s.coeff_valid) .^ 2, [prod(g.sz) 1]);
            use = ~source_bad(s.source_id_valid) & abs(source_share(s.source_id_valid) - 1) < 1e-9;
            share = double(s.coeff_valid(use)) .^ 2;
            [source, ~, group] = unique(double(s.source_id_valid(use)));
            mean_height = accumarray(group, share .* z(dest_layer(use))) ./ accumarray(group, share);
            first_link = accumarray(group, find(use), [], @min);
            expected = z(source_layer(first_link)) + ground_smooth(dest_column(first_link)) - ground_smooth(source_column(first_link));
            testCase.verifyGreaterThan(numel(source), 1000)
            testCase.verifyLessThan(max(abs(mean_height - expected)), 1e-6)
            testCase.verifyGreaterThan(max(abs(ground_smooth(dest_column) - ground_smooth(source_column))), 1)
        end

        function testLinksIntoTheRidgeClimb(testCase)
            % Links whose destination is under the ridge go to its lowest
            % air cell; the sea-level links drop them
            [g, scans, mask] = terrainFollowingTest.synthetic_case(true);
            asl = bmm4d.advect.build_radar_advection_steps(g, scans, 1000, active_mask_3d=mask);
            steps = bmm4d.advect.build_radar_advection_steps(g, scans, 1000, active_mask_3d=mask, ground_smoothing_m=50e3);
            s = steps{1};
            testCase.verifyGreaterThan(nnz(s.climb_valid), 0)
            [dx, dy, dz] = ind2sub(g.sz, double(s.dest_id_valid(s.climb_valid)));
            [~, lowest] = max(mask, [], 3);
            testCase.verifyEqual(dz, lowest(sub2ind(g.sz(1:2), dx, dy)))
            testCase.verifyGreaterThan(numel(unique(s.source_id_valid)), numel(unique(asl{1}.source_id_valid)))
        end
    end

    methods (Static)
        function [g, scans, mask] = synthetic_case(terrain)
            g.sz = [60 20 25];
            g.dx = 4000;
            g.dy = 4000;
            g.dz = 100;
            g.x = (0:g.sz(1) - 1) .* g.dx;
            g.y = (0:g.sz(2) - 1)' .* g.dy;
            g.z = reshape((0:g.sz(3) - 1) .* g.dz + g.dz / 2, 1, 1, []);
            g.mask_water = false(g.sz(1:2));
            g.f_X = @(lon) (lon + 90) * 90000;
            g.f_Y = @(lat) (lat - 35) * 111000;
            g.dem = zeros(g.sz(1:2));
            if terrain
                g.dem = repmat(linspace(0, 800, g.sz(1))', 1, g.sz(2));
                g.dem(30:31, :) = 1500;                    % ridge across the wind
            end
            scan = struct("lon", -90 + g.x(30) / 90000, "lat", 35 + g.y(10) / 111000, "height", 0, ...
                "uv", repmat([4 0], g.sz(3), 1));          % 4 m/s eastward: one cell per 1000 s
            scans = {scan scan};
            ground = max(g.dem, 0);
            mask = g.z + g.dz / 2 > ground & g.z - g.dz / 2 < ground + 2000;
        end
    end
end
