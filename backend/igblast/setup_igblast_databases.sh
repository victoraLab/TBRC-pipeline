#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RAW_DB_ROOT="${SCRIPT_DIR}/db"
BUILD_ROOT="${SCRIPT_DIR}/work"
BIN_ROOT="${SCRIPT_DIR}/bin"

if [[ -x "${BIN_ROOT}/edit_imgt_file.pl" ]]; then
    EDIT_IMGT_FILE="${BIN_ROOT}/edit_imgt_file.pl"
elif command -v edit_imgt_file.pl >/dev/null 2>&1; then
    EDIT_IMGT_FILE="$(command -v edit_imgt_file.pl)"
else
    echo "edit_imgt_file.pl was not found in ${BIN_ROOT} or PATH." >&2
    echo "Place the IgBlast helper script at backend/igblast/bin/edit_imgt_file.pl or install it in PATH." >&2
    exit 1
fi

if ! command -v makeblastdb >/dev/null 2>&1; then
    echo "makeblastdb was not found in PATH." >&2
    exit 1
fi

usage() {
    echo "Usage: $0 [species] [panel]"
    echo "  species: human | mouse | all (default: all)"
    echo "  panel: igh | ig_light | ig_all | tra | trb | tcr_all | all_receptors | all (default: all)"
    exit 1
}

SPECIES_ARG="${1:-all}"
PANEL_ARG="${2:-all}"

case "${SPECIES_ARG}" in
    all) SPECIES_LIST=("human" "mouse") ;;
    human|mouse) SPECIES_LIST=("${SPECIES_ARG}") ;;
    *) usage ;;
esac

case "${PANEL_ARG}" in
    all) PANEL_LIST=("igh" "ig_light" "ig_all" "tra" "trb" "tcr_all" "all_receptors") ;;
    igh|ig_light|ig_all|tra|trb|tcr_all|all_receptors) PANEL_LIST=("${PANEL_ARG}") ;;
    *) usage ;;
esac

families_for_panel() {
    case "$1" in
        igh) echo "IGH" ;;
        ig_light) echo "IGK IGL" ;;
        ig_all) echo "IGH IGK IGL" ;;
        tra) echo "TRA" ;;
        trb) echo "TRB" ;;
        tcr_all) echo "TRA TRB TRD TRG" ;;
        all_receptors) echo "IGH IGK IGL TRA TRB TRD TRG" ;;
        *) return 1 ;;
    esac
}

build_segment() {
    local species="$1"
    local panel="$2"
    local segment="$3"
    local panel_dir="${BUILD_ROOT}/${species}/${panel}"
    local raw_combined_fasta="${panel_dir}/${species}_${panel}_${segment}.raw.fasta"
    local processed_fasta_tmp="${panel_dir}/${species}_${panel}_${segment}.edited.fasta"
    local processed_fasta="${panel_dir}/${species}_${panel}_${segment}.fasta"
    local db_prefix="${panel_dir}/${species}_${panel}_${segment}"
    local found_any=0

    mkdir -p "${panel_dir}"
    : > "${raw_combined_fasta}"

    for family in $(families_for_panel "${panel}"); do
        local source_fasta="${RAW_DB_ROOT}/${species}/${family}${segment}.fasta"
        if [[ -f "${source_fasta}" ]]; then
            cat "${source_fasta}" >> "${raw_combined_fasta}"
            found_any=1
        fi
    done

    if [[ "${found_any}" -eq 0 ]]; then
        rm -f "${raw_combined_fasta}" "${processed_fasta_tmp}" "${processed_fasta}"
        return 0
    fi

    perl "${EDIT_IMGT_FILE}" "${raw_combined_fasta}" > "${processed_fasta_tmp}"

    awk '
        BEGIN {
            RS=">";
            ORS="";
        }
        NR == 1 {
            next;
        }
        {
            split($0, lines, "\n");
            header = lines[1];
            seq_id = header;
            sub(/[[:space:]].*$/, "", seq_id);

            if (!(seq_id in seen)) {
                seen[seq_id] = 1;
                printf(">%s", $0);
            } else {
                printf("Skipping duplicate seq_id during BLAST DB build: %s\n", seq_id) > "/dev/stderr";
            }
        }
    ' "${processed_fasta_tmp}" > "${processed_fasta}"

    makeblastdb -parse_seqids -dbtype nucl -in "${processed_fasta}" -out "${db_prefix}" >/dev/null
}

for species in "${SPECIES_LIST[@]}"; do
    if [[ ! -d "${RAW_DB_ROOT}/${species}" ]]; then
        echo "Skipping species without raw FASTA folder: ${species}" >&2
        continue
    fi

    for panel in "${PANEL_LIST[@]}"; do
        echo "Building ${species}/${panel}"
        build_segment "${species}" "${panel}" "V"
        build_segment "${species}" "${panel}" "D"
        build_segment "${species}" "${panel}" "J"
    done
done

echo "IGBlast databases built under ${BUILD_ROOT}"
