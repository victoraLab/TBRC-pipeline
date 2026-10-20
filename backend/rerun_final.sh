#!/bin/bash
# Rerun the final stage of a finished run -- receptor check, IgBLAST, clonality,
# report, archive -- from its assembled FASTA, without assembling again.
#
# Until this existed the only way to re-annotate, re-report or re-package a run
# was to recompute it from FASTQ, because the final stage's cleanup deletes the
# per-well tree its Snakemake rule depends on (handoff 2026-09-17, section 2).
# It is what an IgBLAST database or panel fix needs, and what RUN15_kappa_vsearch
# was rerun from scratch for.
#
# The previous archive, report, tables and QC are copied to
# <results>/previous-<timestamp>/ first: nothing of the earlier delivery is
# moved or removed, and the new outputs are written in place. The new archive
# stays on the cluster unless --push is given; with it, it is pushed like the
# original and replaces the archive of the same name on the storage server.
#
# Usage: rerun_final.sh <user>-<run> [--push]
#   <user>-<run> is the run's working directory under .parallel/, as the
#   launchers name it (`ls .parallel/` lists them).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KEY="${1:?usage: rerun_final.sh <user>-<run> [--push]}"
PUSH="${2:-}"
WORK="${ROOT}/.parallel/${KEY}"

[[ -d "${WORK}" ]] || { echo "No run directory ${WORK}." >&2; exit 1; }
# Runs started before run.env existed: the launchers name the table
# data/runs/<user>/<run>/sample.json, and <user>-<run> is the key.
if [[ -f "${WORK}/run.env" ]]; then
    # shellcheck source=/dev/null
    source "${WORK}/run.env"
else
    user="${KEY%%-*}"; run="${KEY#*-}"
    DATA_TABLE="data/runs/${user}/${run}/sample.json"
    INPUTPATH="data/runs/${user}/${run}/"
fi
[[ -f "${ROOT}/${DATA_TABLE}" ]] || { echo "No submission at ${ROOT}/${DATA_TABLE}." >&2; exit 1; }

DEST="$(python3 -c 'import json, sys; d = json.load(open(sys.argv[1])); d = d[0] if isinstance(d, list) else d; print(d["dest"])' "${ROOT}/${DATA_TABLE}")"
results_dirs=("${ROOT}/${DEST}"/*/)
if (( ${#results_dirs[@]} != 1 )) || [[ ! -d "${results_dirs[0]}" ]]; then
    echo "Expected one sample folder under ${ROOT}/${DEST}; found ${#results_dirs[@]}." >&2
    exit 1
fi
RESULT="${results_dirs[0]%/}"
SAMPLE="$(basename "${RESULT}")"
[[ -f "${RESULT}/${SAMPLE}.final.fasta" ]] || { echo "No ${SAMPLE}.final.fasta in ${RESULT}; nothing to rerun from." >&2; exit 1; }

KEEP="${RESULT}/previous-$(date +%Y%m%d-%H%M%S)"
mkdir -p "${KEEP}"
for item in "${SAMPLE}".zip "${SAMPLE}".tar.gz "${SAMPLE}".report.html "${SAMPLE}".report.pdf \
            "${SAMPLE}".imgt.ready.fasta igblast clonality qc; do
    [[ -e "${RESULT}/${item}" ]] && cp -a "${RESULT}/${item}" "${KEEP}/"
done
echo "Previous outputs copied to ${KEEP}."

extra=(--forcerun refinal --config "refinal_push=$([[ "${PUSH}" == "--push" ]] && echo 1 || echo 0)")
cd "${ROOT}"
exec ./run_isolated.sh "${KEY}" "${DATA_TABLE}" "${INPUTPATH}" refinal "${extra[@]}"
