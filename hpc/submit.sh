#!/bin/bash
# Run a MATLAB command on the Vogelwarte HPC (merula) with a frozen snapshot of this repository.
#
# Usage (from anywhere in the repository, on the laptop):
#   hpc/submit.sh [-c CPUS] [-m MEM] [-t TIME] [-p PARTITION] NAME 'MATLAB command'
#   hpc/submit.sh sandbox 'case_name = "sandbox"; run(fullfile(project_root, "examples", "run_case.m"))'
# Rules and locations: hpc/README.md
#
# What it does:
#   1. Copies the code as it is now on the laptop (tracked and untracked, not ignored, files;
#      media files left out) to <JOBS>/<yyyymmdd_HHMMSS>_<NAME>/code on the cluster.
#   2. Saves provenance next to it: git commit, branch, status, diff of uncommitted changes,
#      the command (command.m) and the submit options.
#   3. Adds the cluster's own local_settings.json, makes the snapshot read-only, and submits
#      hpc/job.sh from the snapshot with sbatch. The job only ever runs that snapshot.
#
# The code is edited on the laptop only. The cluster never pulls or edits it; every job has its
# own snapshot, so later edits never change a running or queued job.

set -euo pipefail

HOST=merula
JOBS=/mnt/ceph/DOM_Forschung/UNIT_Vogelzug/BMM-4D/jobs   # laptop: /Volumes/UNIT_Vogelzug/BMM-4D/jobs
CPUS=8; MEM=32G; TIME=1-00:00:00; PARTITION=normal

usage() { sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts "c:m:t:p:h" opt; do
  case $opt in
    c) CPUS=$OPTARG ;; m) MEM=$OPTARG ;; t) TIME=$OPTARG ;; p) PARTITION=$OPTARG ;; *) usage ;;
  esac
done
shift $((OPTIND - 1))
[[ $# -eq 2 ]] || usage
NAME=$1; COMMAND=$2
[[ "$NAME" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "NAME may only contain letters, digits, _ and -"; exit 1; }

ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$ROOT"
JOB_ID=$(date +%Y%m%d_%H%M%S)_$NAME
REMOTE=$JOBS/$JOB_ID
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

# Files of the snapshot: what git would see (tracked + untracked, minus .gitignore), without media
# and without tracked files deleted in the working tree
git ls-files -co --exclude-standard | grep -vE '\.(png|gif|jpe?g|pdf|mp4|html|fig|mltbx)$' \
  | while IFS= read -r f; do [[ -f "$f" ]] && printf '%s\n' "$f"; done > "$STAGE/files.txt" || true

# Provenance
{
  echo "job:        $JOB_ID"
  echo "submitted:  $(date -u +%Y-%m-%dT%H:%M:%SZ) by $(whoami)@$(hostname -s)"
  echo "commit:     $(git rev-parse HEAD)"
  echo "branch:     $(git rev-parse --abbrev-ref HEAD)"
  echo "options:    -c $CPUS -m $MEM -t $TIME -p $PARTITION"
  echo "files:      $(wc -l < "$STAGE/files.txt")"
  echo "uncommitted changes (git status --short):"
  git status --short | sed 's/^/  /'
} > "$STAGE/provenance.txt"
git diff HEAD > "$STAGE/uncommitted.diff"
printf '%s\n' "$COMMAND" > "$STAGE/command.m"

# Copy the snapshot and provenance, add the cluster settings, freeze the code
ssh -o BatchMode=yes $HOST "mkdir -p '$REMOTE/code'"
rsync -a --files-from="$STAGE/files.txt" "$ROOT/" "$HOST:$REMOTE/code/"
rsync -a "$STAGE/provenance.txt" "$STAGE/uncommitted.diff" "$STAGE/command.m" "$HOST:$REMOTE/"
ssh -o BatchMode=yes $HOST "cp ~/bmm4d/local_settings.json '$REMOTE/code/' && chmod -R a-w '$REMOTE/code'"

# Submit from the snapshot
SLURM_ID=$(ssh -o BatchMode=yes $HOST "cd '$REMOTE' && sbatch --parsable --job-name='$NAME' \
  --partition=$PARTITION --cpus-per-task=$CPUS --mem=$MEM --time=$TIME \
  --output='$REMOTE/slurm-%j.out' code/hpc/job.sh '$REMOTE'")
echo "slurm:      $SLURM_ID" | ssh -o BatchMode=yes $HOST "cat >> '$REMOTE/provenance.txt'"

echo "Submitted $JOB_ID as Slurm job $SLURM_ID"
echo "  cluster: $REMOTE"
echo "  laptop:  /Volumes/UNIT_Vogelzug/BMM-4D/jobs/$JOB_ID"
echo "  status:  ssh $HOST squeue -j $SLURM_ID    log: slurm-$SLURM_ID.out"
