function plan = buildfile
%BUILDFILE Build tasks of BMM-4D (run with buildtool).
%   buildtool              check and test (default)
%   buildtool check        Code Analyzer on all MATLAB files; fails on errors
%   buildtool test         unit tests (without the Regression tag)
%   buildtool regression   paper-v1 regression test (needs local data, ~1 min)
%   buildtool package      toolbox file release/BMM-4D_<version>.mltbx

plan = buildplan(localfunctions);
plan.DefaultTasks = ["check" "test"];
end

function checkTask(~)
% Code Analyzer on the package, examples, tests and studies
project_root = fileparts(mfilename("fullpath"));
folders = fullfile(project_root, ["toolbox" "examples" "tests" "studies"]);
issues = codeIssues(folders);
errors = issues.Issues(issues.Issues.Severity == "error", :);
fprintf("Code Analyzer: %d errors, %d warnings, %d info\n", height(errors), ...
    nnz(issues.Issues.Severity == "warning"), nnz(issues.Issues.Severity == "info"))
if height(errors) > 0
    disp(errors(:, ["FullFilename" "LineStart" "Description"]))
    error("buildfile:CodeErrors", "Code Analyzer found %d errors.", height(errors))
end
end

function testTask(~)
% Unit tests; the regression test runs in its own task
project_root = setup_paths();
suite = testsuite(fullfile(project_root, "tests"), IncludeSubfolders=true);
suite = suite(~arrayfun(@(test) any(test.Tags == "Regression"), suite));
results = run(suite);
disp(table(results))
assertSuccess(results);
end

function regressionTask(~)
% Today's code against the paper-v1 reference (tests/paperRegressionTest)
project_root = setup_paths();
results = runtests(fullfile(project_root, "tests", "paperRegressionTest.m"));
disp(table(results))
assertSuccess(results);
end

function packageTask(~)
% Toolbox file of toolbox/ with the version and metadata of mpackage.json
project_root = fileparts(mfilename("fullpath"));
toolbox_dir = fullfile(project_root, "toolbox");
manifest = jsondecode(fileread(fullfile(toolbox_dir, "resources", "mpackage.json")));
options = matlab.addons.toolbox.ToolboxOptions(toolbox_dir, manifest.id);
options.ToolboxName = manifest.displayName;
options.ToolboxVersion = manifest.version;
options.Summary = manifest.summary;
options.AuthorName = manifest.provider.name;
options.AuthorCompany = manifest.provider.organization;
options.MinimumMatlabRelease = "R2024b";
options.ToolboxMatlabPath = toolbox_dir;
release_dir = fullfile(project_root, "release");
if ~isfolder(release_dir)
    mkdir(release_dir)
end
options.OutputFile = fullfile(release_dir, "BMM-4D_" + manifest.version + ".mltbx");
matlab.addons.toolbox.packageToolbox(options)
fprintf("Packaged %s\n", options.OutputFile)
end
