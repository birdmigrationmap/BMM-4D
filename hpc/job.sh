#!/bin/bash
# Slurm job run from a snapshot made by hpc/submit.sh; do not sbatch it directly.
# Runs <snapshot>/command.m in MATLAB from <snapshot>/code, after setup_paths, on merops-02
# (the node the MATLAB licence is activated for).
#SBATCH --nodelist=merops-02

set -euo pipefail
SNAP=${1:?snapshot folder}
cd "$SNAP/code"

# Versions of the external code (installed once on the cluster, not part of the snapshot)
{
  echo "node:       $(hostname), slurm job ${SLURM_JOB_ID:-}, ${SLURM_CPUS_PER_TASK:-?} CPUs"
  for d in wsrlib regtools FastGaussianSimulation; do
    echo "$d: $(git -C "$HOME/bmm4d/deps/$d" rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout')"
  done
  echo "started:    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >> "$SNAP/provenance.txt"

# Keep MATLAB within the CPUs Slurm gave the job, then run the command after setup_paths
# (eval of the text, not run(), which would change the current folder away from project_root).
# The root goes on the path so that scripts run with run() (which cd's to their folder) still
# find setup_paths, as AGENTS.md prescribes.
status=0
"$HOME/matlab/matlab" -batch "maxNumCompThreads(str2double(getenv('SLURM_CPUS_PER_TASK'))); \
addpath(pwd); project_root = setup_paths(); eval(fileread('$SNAP/command.m'))" || status=$?

echo "finished:   $(date -u +%Y-%m-%dT%H:%M:%SZ), exit $status" >> "$SNAP/provenance.txt"
exit $status
