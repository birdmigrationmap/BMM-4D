classdef radarFilterOrderTest < matlab.unittest.TestCase
    % read_rd applies the cleaning filters given in the filters option, in
    % that order. Uses one KBYX volume of the regression case; skipped
    % without it.

    properties
        RawFile
    end

    methods (TestClassSetup)
        function findVolume(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(fileparts(mfilename("fullpath"))))));
            setup_paths();
            reference_file = bmm4d.util.data_path("regression", "paper_v1_reference.mat");
            testCase.assumeTrue(isfile(reference_file), "No regression reference.")
            loaded = load(reference_file, "reference");
            testCase.RawFile = loaded.reference.raw_files(1);
            testCase.assumeTrue(isfile(testCase.RawFile), "Raw KBYX volume not found.")
        end
    end

    methods (Test)
        function testExplicitDefaultOrderMatchesOptions(testCase)
            options = {"removeRain", "rho", "parallel", false};
            implicit = bmm4d.radar.read_rd(testCase.RawFile, options{:});
            explicit = bmm4d.radar.read_rd(testCase.RawFile, options{:}, filters=implicit{1}.filter_log.filter);

            testCase.verifyEqual(radarFilterOrderTest.densities(explicit), radarFilterOrderTest.densities(implicit))
        end

        function testLeavingOutAFilterMatchesItsOption(testCase)
            options = {"removeRain", "rho", "parallel", false};
            default_run = bmm4d.radar.read_rd(testCase.RawFile, options{:});
            filters = default_run{1}.filter_log.filter;
            without_smoothing = bmm4d.radar.read_rd(testCase.RawFile, options{:}, filters=filters(filters ~= "smooth_median"));
            smoothing_off = bmm4d.radar.read_rd(testCase.RawFile, options{:}, smoothing=0);

            testCase.verifyEqual(radarFilterOrderTest.densities(without_smoothing), radarFilterOrderTest.densities(smoothing_off))
            testCase.verifyNotEqual(radarFilterOrderTest.densities(without_smoothing), radarFilterOrderTest.densities(default_run))
        end
    end

    methods (Static, Access=private)
        function values = densities(scans)
            values = arrayfun(@(sweep) sweep.data(:), scans{1}.dens.sweeps, UniformOutput=false);
            values = vertcat(values{:});
            values(isnan(values)) = -1;   % verifyEqual treats NaN as unequal
        end
    end
end
