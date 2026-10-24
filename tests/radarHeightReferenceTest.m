classdef radarHeightReferenceTest < matlab.unittest.TestCase
    %RADARHEIGHTREFERENCETEST Height conventions between radar and grid.
    % Velocity profiles are binned above the antenna while grid layers are
    % above sea level, and observation rows must keep their beam inside the
    % grid (kernel rows sum to 1).

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            % The bmm4d package, installed in place for the session (as setup_paths)
            if isempty(mpmlist("bmm4d"))
                mpminstall(fullfile(testCase.ProjectRoot, "toolbox"), InPlace=true, Temporary=true, Prompt=false, Verbosity="quiet");
            end
        end
    end

    methods (Test)
        function testGroundReferencedProfileIsUnchanged(testCase)
            g = struct("z", reshape(50:100:450, 1, 1, []), "dz", 100);
            scan = struct("height", 0, "uv", [(1:5)' (11:15)']);
            testCase.verifyEqual(bmm4d.prior.radar_uv_on_grid(scan, g), scan.uv, AbsTol=1e-12)
        end

        function testProfileIsShiftedByAntennaHeight(testCase)
            % Antenna at 300 m: the grid layer centred at 350 m ASL is 50 m
            % above the antenna (first bin); layers below take the first bin.
            g = struct("z", reshape(50:100:650, 1, 1, []), "dz", 100);
            scan = struct("height", 300, "uv", [(1:5)' (11:15)']);
            uv_grid = bmm4d.prior.radar_uv_on_grid(scan, g);
            testCase.verifyEqual(uv_grid(:, 1), [1 1 1 1 2 3 4]', AbsTol=1e-12)
            testCase.verifyEqual(uv_grid(:, 2), [11 11 11 11 12 13 14]', AbsTol=1e-12)
        end

        function testPartialBeamRowsAreDropped(testCase)
            % Row 1 lies fully in the grid, row 2 keeps only 60% of its beam.
            F = sparse([1 1 2], [1 2 3], [0.5 0.5 0.6], 2, 3);
            obs_matrix = sparse([1 2], [1 1], [10 20], 2, 1);
            info = struct("obs_grid_size", [2 1 1 1]);
            observations = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info, min_row_mass=0.75);
            testCase.verifyEqual(double(observations.row_ids), 1)
            observations = bmm4d.forward.prepare_radar_observations(F, obs_matrix, info);
            testCase.verifyEqual(double(observations.row_ids), [1; 2])
        end
    end
end
