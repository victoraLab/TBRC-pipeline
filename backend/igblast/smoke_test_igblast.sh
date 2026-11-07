#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BACKEND_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUN_IGBLAST_SCRIPT="${BACKEND_DIR}/workflow/scripts/run_igblast.sh"

QUERY_FASTA="${1:-}"
SPECIES="${2:-mouse}"
PANEL="${3:-ig_all}"
SAMPLE_NAME="${4:-igblast_smoke_test}"
OUTPUT_DIR="${5:-${SCRIPT_DIR}/work/smoke_test}"

if [[ -z "${QUERY_FASTA}" ]]; then
    echo "Usage: $0 <query_fasta> [species] [panel] [sample_name] [output_dir]" >&2
    exit 1
fi

if [[ ! -f "${QUERY_FASTA}" ]]; then
    echo "Query FASTA not found: ${QUERY_FASTA}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "Running IGBlast smoke test"
echo "  Query FASTA: ${QUERY_FASTA}"
echo "  Species: ${SPECIES}"
echo "  Panel: ${PANEL}"
echo "  Output dir: ${OUTPUT_DIR}"
echo "  igblastn: ${IGBLAST_BIN:-igblastn}"
echo "  IGBLAST_DATA: ${IGBLAST_DATA:-${SCRIPT_DIR}/internal_data}"
echo

"${RUN_IGBLAST_SCRIPT}" "${QUERY_FASTA}" "${OUTPUT_DIR}" "${SAMPLE_NAME}" "${SPECIES}" "${PANEL}"

OUTPUT_TSV="${OUTPUT_DIR}/${SAMPLE_NAME}.igblast.tsv"
if [[ ! -s "${OUTPUT_TSV}" ]]; then
    echo "IGBlast smoke test failed: output TSV was not created." >&2
    exit 1
fi

echo "IGBlast smoke test completed."
echo "Output: ${OUTPUT_TSV}"
echo
head -n 5 "${OUTPUT_TSV}" || true
