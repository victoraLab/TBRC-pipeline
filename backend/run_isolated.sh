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
# Usage: run_isolated.sh <run-name> <DATA_TABLE> <INPUTPATH>
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
exec ./snakemake.sh
