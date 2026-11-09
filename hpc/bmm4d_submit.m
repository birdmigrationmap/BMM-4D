function job = bmm4d_submit(name, command, options)
%BMM4D_SUBMIT Run a MATLAB command on the HPC cluster with a snapshot of this repository.
%   JOB = BMM4D_SUBMIT(NAME, COMMAND) calls hpc/submit.sh: the code as it is
%   now on this laptop (with uncommitted changes) is copied to the cluster,
%   frozen, and COMMAND runs there after setup_paths (project_root is set).
%   JOB holds the Slurm id and the job folder on the cluster and on this
%   laptop (/Volumes/UNIT_Vogelzug), with provenance.txt, the snapshot and
%   the Slurm log.
%
%   Options: CPUs (8), Memory ("32G"), Time ("1-00:00:00"), Partition ("normal").
%
%   Example:
%     job = bmm4d_submit("sandbox", "case_name = ""sandbox""; run(fullfile(project_root, ""examples"", ""run_case.m""))");
%     type(job.log)   % once it runs
%
%   Rules and locations: hpc/README.md

arguments
    name (1,1) string
    command (1,1) string
    options.CPUs (1,1) double = 8
    options.Memory (1,1) string = "32G"
    options.Time (1,1) string = "1-00:00:00"
    options.Partition (1,1) string = "normal"
end

script = fullfile(fileparts(mfilename("fullpath")), "submit.sh");
call = sprintf("'%s' -c %d -m %s -t %s -p %s %s %s", script, options.CPUs, options.Memory, ...
    options.Time, options.Partition, name, shell_quote(command));
[status, out] = system(call);
fprintf("%s", out)
if status ~= 0
    error("bmm4d_submit:Failed", "Submission failed (see output above).")
end

job.slurm_id = string(regexp(out, "Slurm job (\d+)", "tokens", "once"));
job.cluster_dir = string(regexp(out, "cluster: (\S+)", "tokens", "once"));
job.laptop_dir = string(regexp(out, "laptop:  (\S+)", "tokens", "once"));
job.log = fullfile(job.laptop_dir, "slurm-" + job.slurm_id + ".out");
end

function quoted = shell_quote(text)
% Single-quote for the shell, keeping single quotes in the MATLAB command
quoted = "'" + replace(text, "'", "'\''") + "'";
end
