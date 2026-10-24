classdef mistnetPolarRadarTest < matlab.unittest.TestCase
    %MISTNETPOLARRADARTEST Tests Cartesian-to-polar MistNet remapping.

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            project_root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(project_root));
            setup_paths(project_root);
            fixture_dir = fullfile(fileparts(mfilename("fullpath")), "fixtures");
            addpath(fixture_dir, "-begin")
            testCase.addTeardown(@() rmpath(fixture_dir));
        end
    end

    methods (Test)
        function testPredictionsFollowYXCoordinates(testCase)
            radar = mistnetPolarRadarTest.testRadar();

            preds = bmm4d.radar.mistnet_polar_radar(radar);

            testCase.verifyEqual(preds{1}, [32 23 12 21])
        end

        function testProbabilitiesFollowYXCoordinates(testCase)
            radar = mistnetPolarRadarTest.testRadar();

            [~, ~, probs] = bmm4d.radar.mistnet_polar_radar(radar);

            testCase.verifyEqual(probs{1}(:,:,1), [32 23 12 21], AbsTol=0)
            testCase.verifyEqual(probs{1}(:,:,2), [132 123 112 121], AbsTol=0)
        end

        function testOutsideMistnetGridIsUnclassified(testCase)
            radar = mistnetPolarRadarTest.testRadar(2000);

            preds = bmm4d.radar.mistnet_polar_radar(radar);

            testCase.verifyTrue(all(isnan(preds{1}), "all"))
        end
    end

    methods (Static, Access=private)
        function radar = testRadar(range)
            if nargin == 0
                range = 1000;
            end
            sweep.data = zeros(1, 4);
            sweep.range_bin1 = range;
            sweep.gate_size = 250;
            sweep.nbins = 1;
            sweep.azim_v = [0 90 180 270];
            sweep.elev = 0;
            radar.dz.sweeps = sweep;
        end
    end
end
