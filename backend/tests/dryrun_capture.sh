#!/bin/bash
# Print, for each of the four pipelines, every shell command Snakemake would run.
#
#   bash backend/tests/dryrun_capture.sh BACKEND_DIR OUT_DIR
#
# For each pipeline in the list below this runs `snakemake -n -p --cores 1` over a fixture
# sample.json (backend/tests/fixtures/<pipeline>.sample.json) and fake FASTQs, inside the TBRC
# backend image, and writes the output to OUT_DIR/<pipeline>.txt. Nothing is executed: a dry run
# only prints the commands, which is what compare_dryrun.py reads.
#
# Why it exists: quoting every value in the rules must change nothing the commands do. Run this
# against a tree from before the change and a tree after it, then compare_dryrun.py the two folders.
# BACKEND_DIR is only read (mounted read-only); each pipeline runs in a fresh directory inside the
# container that is discarded afterwards, the same layout run_isolated.sh builds, so the printed
# paths do not depend on where BACKEND_DIR lives.
#
# The fixtures are invented (user tbrc_test, fake sample names) and valid under
# backend/config/trusted_servers.json, so they load in a tree that checks sample.json and in one
# that does not. Several carry empty values on purpose: an empty value is where quoting goes wrong.
#
# The image needs pandas and Snakemake 7.32. Override it with TBRC_BACKEND_IMAGE.
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 BACKEND_DIR OUT_DIR" >&2
    exit 2
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND="$(cd "$1" && pwd)"
mkdir -p "$2"
OUT="$(cd "$2" && pwd)"
IMAGE="${TBRC_BACKEND_IMAGE:-tbrc-backend:2026-11-06}"
PIPELINES=(classic_luka tiago_vsearch tiago_nanopore tiago_sulkseq)

# Runs inside the container. $1 is the pipeline. Two FASTQs per Illumina run (R1/R2); nanopore
# takes the files as they come, so it gets two chunks and exercises the "more than one input" paths.
read -r -d '' INSIDE <<'EOF' || true
set -euo pipefail
export HOME=/tmp
pipeline="$1"
run="/work/${pipeline}"
mkdir -p "${run}/in" "${run}/logs"
cd "${run}"
for item in workflow config data; do
    ln -s "/backend/${item}" "${item}"
done
sample="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['sample_name'])" "/fixtures/${pipeline}.sample.json")"
if [[ "${pipeline}" == "tiago_nanopore" ]]; then
    touch "in/${sample}_chunk_0.fastq.gz" "in/${sample}_chunk_1.fastq.gz"
else
    touch "in/${sample}_S1_R1_001.fastq.gz" "in/${sample}_S1_R2_001.fastq.gz"
fi
export DATA_TABLE="/fixtures/${pipeline}.sample.json" INPUTPATH="${run}/in"
snakemake -n -p --cores 1
EOF

for pipeline in "${PIPELINES[@]}"; do
    docker run --rm \
        -v "${BACKEND}:/backend:ro" \
        -v "${HERE}/fixtures:/fixtures:ro" \
        "${IMAGE}" bash -lc "${INSIDE}" _ "${pipeline}" > "${OUT}/${pipeline}.txt" 2>&1 \
        || { echo "dry run of ${pipeline} failed; see ${OUT}/${pipeline}.txt" >&2; exit 1; }
    if [[ ! -s "${OUT}/${pipeline}.txt" ]]; then
        echo "dry run of ${pipeline} printed nothing" >&2
        exit 1
    fi
done
echo "Wrote ${#PIPELINES[@]} dry runs to ${OUT}."
