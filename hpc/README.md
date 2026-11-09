# Running BMM-4D on the Vogelwarte HPC cluster

A job runs a MATLAB command on node `merops-02` of the cluster `merula`, from a frozen snapshot of
this repository. Submit from the laptop, in the repository:

```bash
hpc/submit.sh sandbox 'case_name = "sandbox"; run(fullfile(project_root, "examples", "run_case.m"))'
hpc/submit.sh -c 16 -m 64G -t 2-00:00:00 -p normal night 'my_script'    # CPUs, memory, time, partition
```

or from MATLAB on the laptop:

```matlab
addpath(fullfile(project_root, "hpc"))
job = bmm4d_submit("sandbox", "case_name = ""sandbox""; run(fullfile(project_root, ""examples"", ""run_case.m""))");
type(job.log)
```

The command runs after `project_root = setup_paths();`, from the snapshot's root. Defaults: 8 CPUs,
32 GB, 1 day, partition `normal`.

## Rules

1. **Edit code on the laptop only.** Nothing is cloned, pulled or edited on the cluster. Each job
   gets its own read-only snapshot of the working tree at submission (tracked and untracked files
   not ignored by git, without images and other media), so later edits never reach a queued or
   running job.
2. **Every job is traceable**: its folder holds `provenance.txt` (commit, branch, uncommitted files,
   options, node, versions of the external code, start, end, exit code), `uncommitted.diff`,
   `command.m`, `code/` and the Slurm log `slurm-<id>.out`.
3. **Settings are per machine**: jobs use the cluster's `~/bmm4d/local_settings.json`, never the
   laptop's.
4. **Data stays in place** on the shared Vogelzug storage, which the laptop also mounts.
5. **External code is fixed**: wsrlib, regtools and FastGaussianSimulation are installed once in
   `~/bmm4d/deps` (wsrlib compiled for Linux by `~/bmm4d/deps/compile_wsrlib.sh`). Change them on
   purpose, then recompile.

## Where things are

| What | Cluster (`/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/…`) | Laptop (`/Volumes/UNIT_Vogelzug/…`) |
| --- | --- | --- |
| Job folders | `BMM-4D/jobs/<date>_<name>/` | same |
| Raw NEXRAD volumes | `weather_radar_us/nexrad_level2` (`nexrad_raw_dir`) | same (`nexrad_shared_dir`, read only) |
| Data folder (`data_dir`) | `BMM-4D/data` (copy of the laptop's `data/`) | same |
| Products, figures | `BMM-4D/products`, `BMM-4D/figures` | same |

On the cluster only: settings `~/bmm4d/local_settings.json`, external code `~/bmm4d/deps/`, MATLAB
`~/matlab/` (see `~/matlab/README`).

## Pitfalls

- **Jobs must run on `merops-02`**, the only node the MATLAB licence is activated for (`job.sh`
  pins it).
- **Parallel workers**: `clean_volumes`, `build_observations` and `static_clutter` start 6 workers
  (fixed in the code), hence 8 CPUs by default; more CPUs only help multithreaded MATLAB steps.
- **Reference data are a copy**: after regenerating `data/reference` or `data/derived` on the
  laptop, copy them again (`rsync -a data/reference data/derived merula:/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/BMM-4D/data/`).
- **New system calls**: the MATLAB image only has the programs listed in
  `~/apptainer/matlab-build-r2026b.def`. Code that calls another program (`system(...)`) fails on
  the cluster until it is added there.
- **Figures on Linux**: MATLAB has no MPEG-4 video writer there, so `export_maps` writes
  `column_density.avi` instead of `.mp4`; check `VideoWriter.getProfiles()` before using MPEG-4
  elsewhere. Helvetica is not installed, so cluster figures use other fonts.
- **Absolute laptop paths** in saved files don't exist on the cluster: `radarFilterOrderTest` is
  skipped there because the regression reference stores the laptop path of its raw volume.
- **Deleting a job folder**: `code/` is read-only, so `chmod -R u+w <job>` first.
  `provenance.txt` and `uncommitted.diff` are enough to rebuild a snapshot from git.
