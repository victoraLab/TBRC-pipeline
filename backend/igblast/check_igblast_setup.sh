#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RAW_DB_ROOT="${SCRIPT_DIR}/db"
WORK_ROOT="${SCRIPT_DIR}/work"
REFS_ROOT="${SCRIPT_DIR}/refs"
BIN_ROOT="${SCRIPT_DIR}/bin"
DEFAULT_INTERNAL_DATA_ROOT="${SCRIPT_DIR}/internal_data"

SPECIES="${1:-mouse}"
PANEL="${2:-ig_all}"
IGBLAST_BIN="${IGBLAST_BIN:-igblastn}"

normalize_optional() {
    local value="${1:-}"
    local lowered
    lowered="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    if [[ -z "${value}" || "${lowered}" == "na" || "${lowered}" == "nan" || "${lowered}" == "null" || "${lowered}" == "none" ]]; then
        return 0
    fi
    printf '%s' "${value}"
}

SPECIES="$(normalize_optional "${SPECIES}")"
PANEL="$(normalize_optional "${PANEL}")"
IGBLAST_BIN="$(normalize_optional "${IGBLAST_BIN}")"
SPECIES="${SPECIES:-mouse}"
PANEL="${PANEL:-ig_all}"
IGBLAST_BIN="${IGBLAST_BIN:-igblastn}"

families_for_panel() {
    case "$1" in
        igh) echo "IGH" ;;
        ig_light) echo "IGK IGL" ;;
        ig_all) echo "IGH IGK IGL" ;;
        tra) echo "TRA" ;;
        trb) echo "TRB" ;;
        tcr_all) echo "TRA TRB TRD TRG" ;;
        all_receptors) echo "IGH IGK IGL TRA TRB TRD TRG" ;;
        *)
            echo "Unknown panel: $1" >&2
            exit 1
            ;;
    esac
}

PANEL_DIR="${WORK_ROOT}/${SPECIES}/${PANEL}"
DB_V="${PANEL_DIR}/${SPECIES}_${PANEL}_V"
DB_D="${PANEL_DIR}/${SPECIES}_${PANEL}_D"
DB_J="${PANEL_DIR}/${SPECIES}_${PANEL}_J"

if [[ -x "${BIN_ROOT}/edit_imgt_file.pl" ]]; then
    EDIT_IMGT_FILE="${BIN_ROOT}/edit_imgt_file.pl"
elif command -v edit_imgt_file.pl >/dev/null 2>&1; then
    EDIT_IMGT_FILE="$(command -v edit_imgt_file.pl)"
else
    EDIT_IMGT_FILE="${BIN_ROOT}/edit_imgt_file.pl"
fi

if [[ -n "${IGBLAST_DATA:-}" && -d "${IGBLAST_DATA}" ]]; then
    INTERNAL_DATA_ROOT="${IGBLAST_DATA}"
else
    INTERNAL_DATA_ROOT="${DEFAULT_INTERNAL_DATA_ROOT}"
fi

if [[ -f "${REFS_ROOT}/${SPECIES}.gl.aux" ]]; then
    AUX_FILE="${REFS_ROOT}/${SPECIES}.gl.aux"
elif [[ -d "${INTERNAL_DATA_ROOT}" && -f "$(cd "$(dirname "${INTERNAL_DATA_ROOT}")" && pwd)/optional_file/${SPECIES}_gl.aux" ]]; then
    AUX_FILE="$(cd "$(dirname "${INTERNAL_DATA_ROOT}")" && pwd)/optional_file/${SPECIES}_gl.aux"
else
    AUX_FILE="${REFS_ROOT}/${SPECIES}.gl.aux"
fi

echo "IGBlast binary: ${IGBLAST_BIN}"
command -v "${IGBLAST_BIN}"
echo
echo "Species: ${SPECIES}"
echo "Panel: ${PANEL}"
echo "Raw DB root: ${RAW_DB_ROOT}/${SPECIES}"
echo "V DB: ${DB_V}"
echo "D DB: ${DB_D}"
echo "J DB: ${DB_J}"
echo "Aux file: ${AUX_FILE}"
echo "Internal data dir: ${INTERNAL_DATA_ROOT}"
echo "edit_imgt_file.pl: ${EDIT_IMGT_FILE}"
echo

missing=0
if [[ ! -x "${EDIT_IMGT_FILE}" ]]; then
    echo "edit_imgt_file.pl not found or not executable at ${EDIT_IMGT_FILE}." >&2
    missing=1
fi

if [[ ! -d "${RAW_DB_ROOT}/${SPECIES}" ]]; then
    echo "Raw FASTA directory not found for species ${SPECIES}: ${RAW_DB_ROOT}/${SPECIES}" >&2
    missing=1
else
    for family in $(families_for_panel "${PANEL}"); do
        if compgen -G "${RAW_DB_ROOT}/${SPECIES}/${family}"'*.fasta' >/dev/null; then
            echo "Raw FASTA present for ${family}."
        else
            echo "Raw FASTA missing for ${family} in ${RAW_DB_ROOT}/${SPECIES}." >&2
            missing=1
        fi
    done
fi

for db_prefix in "${DB_V}" "${DB_J}"; do
    for extension in nhr nin nsq; do
        if [[ ! -f "${db_prefix}.${extension}" ]]; then
            echo "Missing DB index: ${db_prefix}.${extension}" >&2
            missing=1
        fi
    done
done

if [[ -f "${DB_D}.nhr" ]]; then
    echo "D DB present."
else
    echo "D DB not present for this panel."
fi

if [[ -f "${AUX_FILE}" ]]; then
    echo "Aux file present."
else
    echo "Aux file not found for ${SPECIES}." >&2
fi

if [[ -d "${INTERNAL_DATA_ROOT}" ]]; then
    echo "internal_data directory present."
else
    echo "internal_data directory not found." >&2
    missing=1
fi

echo
"${IGBLAST_BIN}" -help >/dev/null
echo "igblastn executable responds to -help."

if [[ "${missing}" -ne 0 ]]; then
    exit 1
fi

echo "IGBlast setup looks ready for ${SPECIES}/${PANEL}."
