classdef configTest < matlab.unittest.TestCase
    % bmm4d.config: defaults, overrides and derived fields.

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
            setup_paths();
            testCase.assumeTrue(isfile(bmm4d.util.data_path("reference", "radar.csv")), "radar.csv not found.")
        end
    end

    methods (Test)
        function testStructSettingReplacesOnlyGivenFields(testCase)
            c = configTest.nightCase("inversion", struct("prior", 1));

            testCase.verifyEqual(c.inversion.prior, 1)
            testCase.verifyEqual(c.inversion.smoothness, 0.3)
        end

        function testReadOptionReplacesOnlyGivenName(testCase)
            c = configTest.nightCase("read_options", {"removeRain", "rho"});

            options = c.read_options;
            testCase.verifyEqual(options{find(strcmp(cellfun(@string, options(1:2:end)), "removeRain")) * 2}, "rho")
            testCase.verifyEqual(numel(options), numel(bmm4d.config("x", radars="KBYX", times=datetime(2018, 4, 23, 3, 0, 0), ...
                grid_bounds=[-83 -80 23 26]).read_options))
        end

        function testDaytimeFramesAreDropped(testCase)
            % KBYX (Key West): 03:00 UTC is night, 18:00 UTC is day
            c = bmm4d.config("x", radars="KBYX", times=datetime(2018, 4, 23, [3 18], 0, 0), grid_bounds=[-83 -80 23 26]);

            testCase.verifyEqual(c.times, datetime(2018, 4, 23, 3, 0, 0))
            testCase.verifyEqual(c.figure_bounds, c.grid_bounds)
        end

        function testUnknownAndMissingSettingsError(testCase)
            testCase.verifyError(@() configTest.nightCase("not_a_setting", 1), "bmm4d:config:UnknownSetting")
            testCase.verifyError(@() bmm4d.config("x", radars="KBYX"), "bmm4d:config:MissingSetting")
        end
    end

    methods (Static, Access=private)
        function c = nightCase(varargin)
            c = bmm4d.config("x", "radars", "KBYX", "times", datetime(2018, 4, 23, 3, 0, 0), ...
                "grid_bounds", [-83 -80 23 26], varargin{:});
        end
    end
end
