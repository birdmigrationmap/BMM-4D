classdef radarFileSelectionTest < matlab.unittest.TestCase
    %RADARFILESELECTIONTEST Local file matching in load_rd (no downloads).
    % Fake empty volume files are created in a temporary working folder,
    % which load_rd searches through data/raw_scans/nexrad.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectFunctions(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "functions"), IncludingSubfolders=true));
        end
    end

    methods (TestMethodSetup)
        function useTemporaryFolder(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            radarFileSelectionTest.touch("KXXX_20200101", "KXXX20200101_020508_V06");
            radarFileSelectionTest.touch("KXXX_20200101", "KXXX20200101_235800_V06");
            radarFileSelectionTest.touch("KXXX_20200102", "KXXX20200102_001000_V06");
        end
    end

    methods (Test)
        function testClosestFileWithinSameDay(testCase)
            files = load_rd("KXXX", datetime(2020, 1, 1, 2, 3, 0), closest_accepted=10);

            testCase.verifyTrue(endsWith(files, "KXXX20200101_020508_V06"))
        end

        function testRejectsFilesBeyondAcceptedOffset(testCase)
            files = load_rd("KXXX", datetime(2020, 1, 1, 3, 0, 0), closest_accepted=10);

            testCase.verifyEqual(files, "")
        end

        function testClosestFileAcrossMidnight(testCase)
            % Requested 00:01 UTC; the 23:58 file of the previous day is
            % 3 min away, the 00:10 file of the same day is 9 min away.
            % load_rd only searches the folder of the requested day.
            files = load_rd("KXXX", datetime(2020, 1, 2, 0, 1, 0), closest_accepted=10);

            testCase.verifyTrue(endsWith(files, "KXXX20200101_235800_V06"), ...
                "load_rd ignores closer files stored in the previous day's folder.")
        end
    end

    methods (Static, Access=private)
        function touch(folder_name, file_name)
            folder = fullfile("data", "raw_scans", "nexrad", folder_name);
            if ~isfolder(folder)
                mkdir(folder)
            end
            fid = fopen(fullfile(folder, file_name), "w");
            fclose(fid);
        end
    end
end
