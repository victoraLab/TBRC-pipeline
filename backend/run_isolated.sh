#!/bin/bash
# Run a pipeline in an isolated working directory so runs can go in parallel.
#
# snakemake locks .snakemake/ in its working directory, and every run has always
# shared the pipeline root. Two concurrent runs therefore collide: the second
# either refuses to start, or -- worse, and this is what happened during
# development -- both proceed and one run's failure cleanup deletes the other's
# intermediate files, producing a "No such file or directory" that points nowhere
# near the real cause.
#
# The fix is a per-run directory containing symlinks to everything read-only
# (workflow, config, barcodes, igblast databases) and its own .snakemake, logs
# and results. Nothing is copied, so the cost is a few symlinks per run.
#
# goServer.bash, goUpload.bash and `tbrc run` all launch through this script,
# named <user>-<run>, so a second submission no longer waits for the first.
# Two runs that would write the same results/<user>/<sample>/ still collide;
# give them different run names.
#
# Usage: run_isolated.sh <run-name> <DATA_TABLE> <INPUTPATH> [snakemake arguments...]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAME="$1"; TABLE="$2"; INPUT="$3"
WORK="${ROOT}/.parallel/${NAME}"

mkdir -p "${WORK}"
# Read-only, shared: symlink rather than copy.
for item in workflow config igblast snakemake.sh; do
    [[ -e "${WORK}/${item}" ]] || ln -s "${ROOT}/${item}" "${WORK}/${item}"
done
# data/ is shared for barcodes but each run reads its own subtree under runs/,
# so the whole tree is linked and nothing is duplicated.
[[ -e "${WORK}/data" ]] || ln -s "${ROOT}/data" "${WORK}/data"
# Writable and per-run: results land in the real tree so delivery is unchanged,
# but .snakemake and logs must not be shared or the lock problem returns.
[[ -e "${WORK}/results" ]] || ln -s "${ROOT}/results" "${WORK}/results"
mkdir -p "${WORK}/logs" "${WORK}/slurm_outs" "${WORK}/.snakemake"

cd "${WORK}"
export DATA_TABLE="${TABLE}" INPUTPATH="${INPUT}"
# What this run was started with, so rerun_final.sh can find it again.
printf 'DATA_TABLE=%q\nINPUTPATH=%q\n' "${TABLE}" "${INPUT}" > run.env

# A detached launch (dispatch_detached in frontend/scripts/common.sh) names an
# attempt. The driver then outlives the ssh session that started it, so its exit
# code is left in attempts/<attempt>.exit for the launcher to poll. Written
# after the driver has finished, so the log is complete once the file exists.
if [[ -n "${TBRC_ATTEMPT:-}" ]]; then
    mkdir -p attempts
    set +e
    ./snakemake.sh "${@:4}"
    status=$?
    echo "${status}" > "attempts/${TBRC_ATTEMPT}.exit"
    exit "${status}"
fi
exec ./snakemake.sh "${@:4}"
