function run_script_with_log(project_root, script_relative_path, log_file)
%RUN_SCRIPT_WITH_LOG Run a workflow script from the project root with diary logging.

cd(project_root)
addpath(project_root)

log_dir = fileparts(log_file);
if ~isempty(log_dir) && ~isfolder(log_dir)
    mkdir(log_dir)
end

diary(log_file)
cleanup = onCleanup(@() diary("off"));

fprintf("[runner] started %s\n", string(datetime("now", TimeZone="Europe/Paris", Format="yyyy-MM-dd HH:mm:ss z")));
fprintf("[runner] project_root: %s\n", project_root);
fprintf("[runner] script: %s\n", script_relative_path);

run(fullfile(project_root, script_relative_path))

fprintf("[runner] finished %s\n", string(datetime("now", TimeZone="Europe/Paris", Format="yyyy-MM-dd HH:mm:ss z")));
clear cleanup
diary("off")
end
