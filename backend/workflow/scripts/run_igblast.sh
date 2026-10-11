#!/bin/bash

set -euo pipefail

QUERY_FASTA="$1"
OUTPUT_DIR="$2"
SAMPLE_NAME="$3"
SPECIES="${4:-human}"
PANEL="${5:-ig_all}"
DB_V="${6:-}"
DB_D="${7:-}"
DB_J="${8:-}"
AUX_PATH="${9:-}"
ORGANISM="${10:-human}"
IGBLAST_BIN="${IGBLAST_BIN:-igblastn}"
# Exit status for "this host is not set up for IgBLAST", which the caller treats
# as a skip rather than a failure. It used to be 2, and BLAST+ also exits 2 on
# any database error, so a run whose IgBLAST aborted mid-way was recorded as
# "skipped: databases not configured" and delivered with an empty igblast/
# folder. CMP15_kappa, September 2026, is the worked example. 64 is EX_USAGE
# from sysexits.h; BLAST+ uses 0-4 and 255.
SKIP_EXIT=64
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BACKEND_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
IGBLAST_ROOT="${BACKEND_DIR}/igblast"
DB_BUILD_ROOT="${IGBLAST_ROOT}/work"
REFS_ROOT="${IGBLAST_ROOT}/refs"
IGBLAST_DATA_DIR="${IGBLAST_ROOT}/internal_data"

resolve_aux_path() {
    local requested_path="${1:-}"
    local species_name="${2:-}"
    local organism_name="${3:-}"
    local search_dir=""
    local candidate=""
    local -a candidate_names=(
        "${species_name}.gl.aux"
        "${species_name}_gl.aux"
        "${organism_name}.gl.aux"
        "${organism_name}_gl.aux"
    )

    # Allow either an exact aux file or a directory containing species-specific
    # aux files such as mouse_gl.aux and human_gl.aux.
    if [[ -n "${requested_path}" ]]; then
        if [[ -f "${requested_path}" ]]; then
            printf '%s' "${requested_path}"
            return 0
        fi
        if [[ -d "${requested_path}" ]]; then
            search_dir="${requested_path}"
        fi
    fi

    if [[ -z "${search_dir}" && -f "${REFS_ROOT}/${species_name}.gl.aux" ]]; then
        printf '%s' "${REFS_ROOT}/${species_name}.gl.aux"
        return 0
    fi

    # The NCBI IgBlast tarball ships internal_data/ and optional_file/
    # side-by-side, so search next to IGBLAST_DATA when available.
    if [[ -z "${search_dir}" && -n "${IGBLAST_DATA:-}" ]]; then
        candidate="$(cd "$(dirname "${IGBLAST_DATA}")" && pwd)/optional_file"
        if [[ -d "${candidate}" ]]; then
            search_dir="${candidate}"
        fi
    fi

    if [[ -z "${search_dir}" && -d "${IGBLAST_ROOT}/optional_file" ]]; then
        search_dir="${IGBLAST_ROOT}/optional_file"
    fi

    if [[ -n "${search_dir}" ]]; then
        for candidate in "${candidate_names[@]}"; do
            if [[ -f "${search_dir}/${candidate}" ]]; then
                printf '%s' "${search_dir}/${candidate}"
                return 0
            fi
        done
    fi

    return 0
}

seqtype_for_panel() {
    local panel_name="${1:-}"
    case "${panel_name}" in
        igh|ig_light|ig_all|ig_b6)
            printf '%s' "Ig"
            ;;
        tra|trb|tcr_all)
            printf '%s' "TCR"
            ;;
        *)
            return 0
            ;;
    esac
}

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
DB_V="$(normalize_optional "${DB_V}")"
DB_D="$(normalize_optional "${DB_D}")"
DB_J="$(normalize_optional "${DB_J}")"
AUX_PATH="$(normalize_optional "${AUX_PATH}")"
ORGANISM="$(normalize_optional "${ORGANISM}")"
IGBLAST_BIN="$(normalize_optional "${IGBLAST_BIN}")"

SPECIES="${SPECIES:-human}"
PANEL="${PANEL:-ig_all}"
IGBLAST_BIN="${IGBLAST_BIN:-igblastn}"

if [[ "${ORGANISM}" == "NA" || -z "${ORGANISM}" ]]; then
    ORGANISM="${SPECIES}"
fi

PANEL_DIR="${DB_BUILD_ROOT}/${SPECIES}/${PANEL}"
IG_SEQTYPE="$(seqtype_for_panel "${PANEL}")"

# Panels are normally built as <species>_<panel>_<segment>, but some were built
# by hand with the "ig_" prefix dropped (mouse/ig_b6 holds mouse_b6_V and so on).
# Accept either spelling so a hand-built panel is usable without renaming files.
resolve_db() {
    local segment="$1"
    local primary="${PANEL_DIR}/${SPECIES}_${PANEL}_${segment}"
    local alternate="${PANEL_DIR}/${SPECIES}_${PANEL#ig_}_${segment}"

    if [[ -f "${primary}.nhr" ]]; then
        printf '%s' "${primary}"
    elif [[ -f "${alternate}.nhr" ]]; then
        printf '%s' "${alternate}"
    else
        printf '%s' "${primary}"
    fi
}

if [[ -z "${DB_V}" ]]; then
    DB_V="$(resolve_db V)"
fi

if [[ -z "${DB_D}" ]]; then
    DB_D="$(resolve_db D)"
fi

if [[ -z "${DB_J}" ]]; then
    DB_J="$(resolve_db J)"
fi

# igblastn insists on a D database for every Ig query, even a light chain, which
# has no D segment. When the panel ships none (ig_light), leaving -germline_db_D
# out does not skip the D search: igblastn falls back to its built-in default
# name (<organism>_gl_D), fails to find it, and aborts. Lend it this species'
# heavy-chain D set instead; light-chain reads simply do not align to it.
# Measured on CMP15_kappa: the V calls and the locus column are identical with
# and without the lent set, and the abort is gone.
if [[ -z "${DB_D}" || ! -f "${DB_D}.nhr" ]] && [[ "${IG_SEQTYPE}" == "Ig" ]]; then
    for fallback_panel in igh ig_all; do
        candidate="${DB_BUILD_ROOT}/${SPECIES}/${fallback_panel}/${SPECIES}_${fallback_panel}_D"
        if [[ -f "${candidate}.nhr" ]]; then
            echo "IGBlast D DB: panel ${PANEL} has no D segments; using ${candidate}" >&2
            DB_D="${candidate}"
            break
        fi
    done
fi

if [[ -z "${DB_V}" || -z "${DB_J}" ]]; then
    echo "IGBlast skipped: missing V/J database paths." >&2
    exit "${SKIP_EXIT}"
fi

if ! command -v "${IGBLAST_BIN}" >/dev/null 2>&1; then
    echo "IGBlast skipped: igblastn binary not found (${IGBLAST_BIN})." >&2
    exit "${SKIP_EXIT}"
fi

if [[ -n "${IGBLAST_DATA:-}" && -d "${IGBLAST_DATA}" ]]; then
    export IGBLAST_DATA
elif [[ -d "${IGBLAST_DATA_DIR}" ]]; then
    export IGBLAST_DATA="${IGBLAST_DATA_DIR}"
else
    echo "IGBlast skipped: internal_data directory was not found. Place the IgBlast internal_data folder under ${IGBLAST_ROOT}/internal_data or set IGBLAST_DATA." >&2
    exit "${SKIP_EXIT}"
fi

AUX_PATH="$(resolve_aux_path "${AUX_PATH}" "${SPECIES}" "${ORGANISM}")"

mkdir -p "${OUTPUT_DIR}"

IGBLAST_OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}.igblast.tsv"
# AIRR outfmt 19 gives downstream-friendly tabular output for optional
# clonality processing.
IGBLAST_ARGS=(
    -germline_db_V "${DB_V}"
    -germline_db_J "${DB_J}"
    -organism "${ORGANISM}"
    -query "${QUERY_FASTA}"
    -outfmt 19
    -out "${IGBLAST_OUTPUT}"
)

if [[ -n "${DB_D}" && "${DB_D}" != "NA" && -f "${DB_D}.nhr" ]]; then
    IGBLAST_ARGS+=(-germline_db_D "${DB_D}")
fi

if [[ -n "${IG_SEQTYPE}" ]]; then
    IGBLAST_ARGS+=(-ig_seqtype "${IG_SEQTYPE}")
fi

if [[ -n "${AUX_PATH}" && "${AUX_PATH}" != "NA" ]]; then
    IGBLAST_ARGS+=(-auxiliary_data "${AUX_PATH}")
fi

if [[ ! -f "${DB_V}.nhr" || ! -f "${DB_V}.nin" || ! -f "${DB_V}.nsq" || ! -f "${DB_J}.nhr" || ! -f "${DB_J}.nin" || ! -f "${DB_J}.nsq" ]]; then
    echo "IGBlast skipped: built databases were not found for species=${SPECIES}, panel=${PANEL}. Run backend/igblast/setup_igblast_databases.sh first." >&2
    exit "${SKIP_EXIT}"
fi

echo "IGBlast binary: ${IGBLAST_BIN}" >&2
echo "IGBlast species/panel: ${SPECIES}/${PANEL}" >&2
echo "IGBlast seqtype: ${IG_SEQTYPE:-auto}" >&2
echo "IGBlast V DB: ${DB_V}" >&2
echo "IGBlast J DB: ${DB_J}" >&2
echo "IGBlast aux: ${AUX_PATH:-none}" >&2
echo "IGBlast internal_data: ${IGBLAST_DATA}" >&2

"${IGBLAST_BIN}" "${IGBLAST_ARGS[@]}"
echo "IGBlast output: ${IGBLAST_OUTPUT}"
