classdef (TestTags = {'Regression'}) paperRegressionTest < matlab.unittest.TestCase
    % Today's code against the paper-v1 reference on a small Key West case.
    % The reference (data/regression/paper_v1_reference.mat, local) is built
    % by tests/regression/make_regression_reference.m. Restructuring steps
    % must keep both paths unchanged; a deliberate method change fails here
    % and is recorded in MIGRATION.md with the size of the difference.
    % Needs the KBYX raw volumes and the 150 km kernels; skipped without them.
    % Runtime about 2 min (MistNet and two solves).

    properties
        ProjectRoot
        Reference
        KernelDir
    end

    methods (TestClassSetup)
        function loadReference(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
            testCase.ProjectRoot = setup_paths();
            addpath(fullfile(testCase.ProjectRoot, "pipeline", "functions"), ...
                fullfile(testCase.ProjectRoot, "tests", "regression"))
            reference_file = bmm4d.util.data_path("regression", "paper_v1_reference.mat");
            testCase.assumeTrue(isfile(reference_file), "No regression reference; run make_regression_reference.")
            loaded = load(reference_file, "reference");
            testCase.Reference = loaded.reference;
            testCase.assumeTrue(all(isfile(testCase.Reference.raw_files)), "Raw KBYX volumes not found.")
            testCase.KernelDir = bmm4d.util.data_path("derived", "F_elev", "150km");
        end
    end

    methods (Test)
        function testPaperPathUnchanged(testCase)
            current = regression_case_paper(testCase.Reference.raw_files, testCase.KernelDir, ...
                testCase.Reference.paper_settings);
            paperRegressionTest.verifySame(testCase, testCase.Reference.paper, current, 1e-6)
        end

        function testProductionPathUnchanged(testCase)
            work_dir = tempname;
            mkdir(work_dir)
            testCase.addTeardown(@() rmdir(work_dir, "s"))
            current = regression_case_prod(testCase.Reference.raw_files, testCase.Reference.times, ...
                testCase.KernelDir, work_dir, testCase.Reference.prod_settings);
            paperRegressionTest.verifySame(testCase, testCase.Reference.prod, current, 1e-5)
        end
    end

    methods (Static, Access=private)
        function verifySame(testCase, reference, current, tolerance)
            difference = regression_compare(reference, current);
            for i = 1:height(difference)
                testCase.verifyTrue(difference.same_size(i), difference.output(i) + ": size changed")
                testCase.verifyLessThanOrEqual(difference.relative_difference(i), tolerance, ...
                    sprintf("%s: relative difference %.3g", difference.output(i), difference.relative_difference(i)))
            end
        end
    end
end
