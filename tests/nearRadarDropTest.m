classdef nearRadarDropTest < matlab.unittest.TestCase
    % Near-radar cut drop: bmm4d.forward.drop_near_radar_cuts on synthetic
    % volumes, the config default, and bmm4d.radar.read_vcp on one archive
    % file (skipped when the radar share is not mounted).

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
            setup_paths();
        end
    end

    methods (Test)
        function testBatchGatesDroppedOnlyNearRadar(testCase)
            scans = {nearRadarDropTest.volume()};
            [dropped, n_dropped] = bmm4d.forward.drop_near_radar_cuts(scans, "batch", 20000);

            batch = dropped{1}.dens.sweeps(2);
            near = batch.range < 20000;
            testCase.verifyTrue(all(isnan(batch.data(near, :)), "all"))
            testCase.verifyTrue(all(isnan(batch.sd(near, :)), "all"))
            testCase.verifyEqual(batch.data(~near, :), scans{1}.dens.sweeps(2).data(~near, :))
            % Zero-density gates count as observations too
            testCase.verifyEqual(n_dropped, nnz(near) * numel(batch.az))
        end

        function testSplitAndDopplerCutsKept(testCase)
            scans = {nearRadarDropTest.volume()};
            dropped = bmm4d.forward.drop_near_radar_cuts(scans, "batch", 20000);

            testCase.verifyEqual(dropped{1}.dens.sweeps([1 3]), scans{1}.dens.sweeps([1 3]))
        end

        function testEmptyRuleKeepsEverything(testCase)
            scans = {nearRadarDropTest.volume()};
            [dropped, n_dropped] = bmm4d.forward.drop_near_radar_cuts(scans, strings(1, 0), 20000);

            testCase.verifyEqual(dropped, scans)
            testCase.verifyEqual(n_dropped, 0)
        end

        function testSplitCutsCannotBeDropped(testCase)
            testCase.verifyError(@() bmm4d.forward.drop_near_radar_cuts({nearRadarDropTest.volume()}, "split", 20000), ...
                "bmm4d:drop_near_radar_cuts:Split")
        end

        function testVolumeWithoutWaveformFails(testCase)
            volume = nearRadarDropTest.volume();
            volume.dens.sweeps = rmfield(volume.dens.sweeps, "waveform");

            testCase.verifyError(@() bmm4d.forward.drop_near_radar_cuts({volume}, "batch", 20000), ...
                "bmm4d:drop_near_radar_cuts:NoWaveform")
        end

        function testConfigDefault(testCase)
            c = bmm4d.config("x", radars="KDVN", times=datetime(2021, 5, 13, 3, 0, 0), grid_bounds=[-92 -89 40 43]);

            testCase.verifyFalse(any(c.near_radar.waveforms == "split"))
            testCase.verifyGreaterThan(c.near_radar.max_range_m, 0)
        end

        function testReadVcpOfArchiveFile(testCase)
            % KDVN 13 May 2021 01:56 UTC, VCP 35: split cuts 0.5-1.3 deg
            raw_file = "/Volumes/UNIT_Vogelzug/weather_radar_us/nexrad_level2/KDVN_20210513/KDVN20210513_015652_V06";
            if ~isfile(raw_file)
                raw_file = "/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/weather_radar_us/nexrad_level2/KDVN_20210513/KDVN20210513_015652_V06";
            end
            testCase.assumeTrue(isfile(raw_file), "Radar share not mounted.")

            vcp = bmm4d.radar.read_vcp(raw_file);

            testCase.verifyEqual(vcp.number, 35)
            testCase.verifyTrue(all(vcp.cuts.split(vcp.cuts.elev < 1.4)))
            testCase.verifyFalse(any(vcp.cuts.split(vcp.cuts.elev > 1.4)))
            testCase.verifyTrue(any(vcp.cuts.waveform_name == "batch"))
        end
    end

    methods (Static)
        function volume = volume()
            % One volume: a split, a batch and a Doppler sweep, gates every km
            range_m = (1000:1000:60000)';
            az = 0.5:1:359.5;
            sweep = @(elev, waveform) struct("elev", elev, "az", az, "range", range_m, ...
                "data", repmat(mod(range_m / 1000, 3), 1, numel(az)), "sd", ones(numel(range_m), numel(az)), ...
                "waveform", waveform);
            sweeps = [sweep(0.5, "split") sweep(2.4, "batch") sweep(6.4, "doppler")];
            volume = struct("station", "KTST", "dens", struct("sweeps", sweeps));
        end
    end
end
