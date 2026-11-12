function ar_run_script(project_root, script_files, n_threads)
%AR_RUN_SCRIPT Run study scripts in a MATLAB batch job (dev).
%   AR_RUN_SCRIPT(PROJECT_ROOT, SCRIPT_FILES, N_THREADS) runs the scripts
%   SCRIPT_FILES (string array) one after the other from PROJECT_ROOT with
%   N_THREADS computational threads (batch workers start with one), so a
%   long study runs beside the client.
maxNumCompThreads(n_threads);
for script_file = string(script_files(:))'
    cd(project_root)
    addpath(project_root)
    run(script_file)
end
end
