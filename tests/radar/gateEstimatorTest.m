classdef gateEstimatorTest < matlab.unittest.TestCase
    % bmm4d.radar.gate: the gate density estimator on sweeps whose answer is
    % known without the simulator of studies/gate_density.

    properties
        Beam
        Table
    end

    methods (TestClassSetup)
        function buildTable(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(fileparts(mfilename("fullpath"))))));
            setup_paths();
            testCase.Beam = bmm4d.radar.gate.beam();
            testCase.Table = bmm4d.radar.gate.bird_sum_table(testCase.Beam, 2000, 20000);
        end
    end

    methods (Test)
        function testFarUniformGatesGiveTheirDensity(testCase)
            % Far out a gate holds thousands of birds, so every gate reads the mean
            range_m = (60000:250:80000)';
            density = 200;
            dbz = gateEstimatorTest.dbzOf(density) * ones(numel(range_m), 720);
            estimate = gateEstimatorTest.estimate(testCase, dbz, false(size(dbz)), range_m, -30);

            testCase.verifyEqual(median(estimate.density(:)), density, RelTol=0.05)
        end

        function testEmptySweepReadsAlmostNoBirds(testCase)
            % No echo anywhere, with a detection limit of one bird on the beam axis at 5 km
            range_m = (5000:250:20000)';
            dbz = nan(numel(range_m), 720);
            estimate = gateEstimatorTest.estimate(testCase, dbz, true(size(dbz)), range_m, -20);

            testCase.verifyLessThan(max(estimate.density(:)), 1)
        end

        function testMissingGatesCarryNoInformation(testCase)
            % Rain-removed gates (NaN, not no-echo) are unknown; a window of them has no estimate.
            % Far out, where identical gates are what the model expects (thousands of birds per gate)
            range_m = (60000:250:80000)';
            full = gateEstimatorTest.dbzOf(50) * ones(numel(range_m), 720);
            dbz = full;
            dbz(:, 1:100) = nan;
            [lo, hi] = bmm4d.radar.gate.bounds(dbz, false(size(dbz)), range_m, testCase.Beam, -30 * ones(size(range_m)), 0.1);
            estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, 0.5, testCase.Beam, testCase.Table);
            reference = gateEstimatorTest.estimate(testCase, full, false(size(full)), range_m, -30);

            testCase.verifyTrue(all(isinf(lo(:, 1:100)) & isinf(hi(:, 1:100)), "all"))
            testCase.verifyTrue(all(isnan(estimate.density(:, 30:70)), "all"))
            % Away from the gap, the same estimate as without it
            testCase.verifyEqual(estimate.density(:, 200:600), reference.density(:, 200:600), RelTol=1e-10)
        end
    end

    methods (Static, Access=private)
        function dbz = dbzOf(density)
            % Reflectivity (dBZ) of a gate whose birds give DENSITY (birds km^-3), as to_density
            dbz = 10 * log10(density * 11 / z_to_refl(1));
        end

        function estimate = estimate(testCase, dbz, no_echo, range_m, detection_dbz)
            [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, testCase.Beam, detection_dbz * ones(size(range_m)), 0.1);
            estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, 0.5, testCase.Beam, testCase.Table);
        end
    end
end
