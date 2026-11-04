classdef beamBlockageTest < matlab.unittest.TestCase
    % bmm4d.radar.filter.beam_blockage against known terrain around KMRX
    % (Knoxville): the Great Smoky Mountains to the south (Clingmans Dome,
    % 2025 m, 68 km at compass bearing 187 deg) and the open Tennessee
    % Valley to the south-west. The maps are indexed by math angle, so a
    % lookup at the compass azimuth itself removes the valley and keeps the
    % gates behind the mountains.

    methods (TestClassSetup)
        function addPaths(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(fileparts(mfilename("fullpath"))))));
            setup_paths();
            testCase.assumeTrue(isfile(bmm4d.util.data_path("reference", "beam_blockage", "beam_blockage.mat")), ...
                "No beam-blockage maps.")
        end
    end

    methods (Test)
        function testGatesBehindTheSmokiesAreRemoved(testCase)
            % Clingmans Dome rises about 1.1 deg above KMRX: a 0.9 deg sweep
            % (lower half-power edge 0.43 deg) is blocked beyond it
            blocked = beamBlockageTest.blockedShare("KMRX", 0.9, (75:0.25:100)' * 1000, 180:190);
            testCase.verifyGreaterThan(blocked, 0.95)
        end

        function testTennesseeValleyIsKept(testCase)
            % Down the valley toward Chattanooga the lowest sweep is clear
            blocked = beamBlockageTest.blockedShare("KMRX", 0.5, (75:0.25:100)' * 1000, 230:260);
            testCase.verifyLessThan(blocked, 0.1)
        end
    end

    methods (Static)
        function share = blockedShare(station, elevation_deg, range_m, azimuth_deg)
            % Share of gates (range x azimuth) the filter removes from one sweep
            gate = struct("range", range_m, "az", azimuth_deg, "elev", elevation_deg, "vert_half_bw", 0.475, ...
                "data", ones(numel(range_m), numel(azimuth_deg)), "no_echo", false(numel(range_m), numel(azimuth_deg)));
            scan = bmm4d.radar.filter.beam_blockage(struct("station", station, "gates", gate));
            share = mean(isnan(scan.gates.data), "all");
        end
    end
end
