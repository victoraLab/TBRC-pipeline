#!/bin/bash

set -euo pipefail

# Usage function
usage() {
    echo "Usage: $0 [-i INPUT] [-o OUTPUT] [-d DESTINY] [-s SAMPLE] [-t TRIMMING_SEQ] [-y SHOULD_TRIM] [-L PREPROCESS_METHOD] [-u SERVER_USER] [-n SERVER_NAME] [-e NOTIFICATION_EMAIL] [-k SHOULD_KEEP] [-p SERVER_RESULTS_PATH] [-j SERVER_STORAGE] [-r SERVER_SSH_PORT] [-m SERVER_TRANSFER_MODE] [-g RUN_IGBLAST] [-q RUN_CLONALITY] [-f ARCHIVE_FORMAT] [-z IGBLAST_SPECIES] [-l IGBLAST_PANEL] [-b IGBLAST_BIN] [-c IGBLAST_ORGANISM] [-v IGBLAST_DB_V] [-w IGBLAST_DB_D] [-x IGBLAST_DB_J] [-I IGBLAST_DATA] [-a IGBLAST_AUX] [-A ARCHIVE_SCOPE] [-R RECEPTOR_CHECK] [-G RUN_GROUP] [-F PRIMER_FOOTPRINT] [-T LIGHT_TAIL_3P]"
    exit 1
}

# Parsing arguments using getopts
while getopts ":i:o:s:d:t:y:L:u:e:k:p:j:r:m:g:q:f:z:l:b:c:v:w:x:I:a:A:P:B:R:G:F:T:" opt; do
    case ${opt} in
        i) INPUT="${OPTARG}" ;;
        o) OUTPUT="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        d) DESTINY="${OPTARG}" ;;
        t) TRIMMING_SEQ="${OPTARG}" ;;
        y) SHOULD_TRIM="${OPTARG}" ;;
        L) PREPROCESS_METHOD="${OPTARG}" ;;
        u) SERVER_USER="${OPTARG}" ;;
        e) NOTIFICATION_EMAIL="${OPTARG}" ;;
        k) SHOULD_KEEP="${OPTARG}" ;;
        p) SERVER_RESULTS_PATH="${OPTARG}" ;;
        j) SERVER_STORAGE="${OPTARG}" ;;
        r) SERVER_SSH_PORT="${OPTARG}" ;;
        m) SERVER_TRANSFER_MODE="${OPTARG}" ;;
        g) RUN_IGBLAST="${OPTARG}" ;;
        q) RUN_CLONALITY="${OPTARG}" ;;
        f) ARCHIVE_FORMAT="${OPTARG}" ;;
        z) IGBLAST_SPECIES="${OPTARG}" ;;
        l) IGBLAST_PANEL="${OPTARG}" ;;
        b) IGBLAST_BIN="${OPTARG}" ;;
        c) IGBLAST_ORGANISM="${OPTARG}" ;;
        v) IGBLAST_DB_V="${OPTARG}" ;;
        w) IGBLAST_DB_D="${OPTARG}" ;;
        x) IGBLAST_DB_J="${OPTARG}" ;;
        I) IGBLAST_DATA="${OPTARG}" ;;
        a) IGBLAST_AUX="${OPTARG}" ;;
        A) ARCHIVE_SCOPE="${OPTARG}" ;;
        P) PIPELINE_NAME="${OPTARG}" ;;
        B) PRIMER_SET_NAME="${OPTARG}" ;;
        R) RECEPTOR_CHECK="${OPTARG}" ;;
        G) RUN_GROUP="${OPTARG}" ;;
        F) PRIMER_FOOTPRINT="${OPTARG}" ;;
        T) LIGHT_TAIL_3P="${OPTARG}" ;;
        *) usage ;;
    esac
done

PREPROCESS_METHOD="${PREPROCESS_METHOD:-standard}"
# Every optional flag gets a default before use. The script runs under `set -u`,
# so a value the caller omitted -- or one lost because an empty argument collapsed
# on the command line -- would otherwise abort with "unbound variable" partway
# through, after the workflow had already done its work.
: "${INPUT:=}"; : "${OUTPUT:=}"; : "${SAMPLE:=}"; : "${DESTINY:=}"
: "${TRIMMING_SEQ:=}"; : "${SHOULD_TRIM:=false}"; : "${PREPROCESS_METHOD:=standard}"
: "${SERVER_USER:=}"; : "${NOTIFICATION_EMAIL:=}"; : "${SHOULD_KEEP:=false}"
: "${SERVER_RESULTS_PATH:=}"; : "${SERVER_STORAGE:=}"; : "${SERVER_SSH_PORT:=}"
: "${SERVER_TRANSFER_MODE:=push}"; : "${RUN_IGBLAST:=false}"; : "${RUN_CLONALITY:=false}"
: "${ARCHIVE_FORMAT:=zip}"; : "${ARCHIVE_SCOPE:=}"
: "${PIPELINE_NAME:=}"; : "${PRIMER_SET_NAME:=}"; : "${RECEPTOR_CHECK:=auto}"; : "${RUN_GROUP:=none}"
: "${PRIMER_FOOTPRINT:=}"; : "${LIGHT_TAIL_3P:=CACTGACACTCCCACGAGCTC}"
[[ "${PRIMER_FOOTPRINT}" == "none" ]] && PRIMER_FOOTPRINT=""
: "${IGBLAST_SPECIES:=}"; : "${IGBLAST_PANEL:=}"; : "${IGBLAST_BIN:=igblastn}"
: "${IGBLAST_ORGANISM:=}"; : "${IGBLAST_DB_V:=}"; : "${IGBLAST_DB_D:=}"
: "${IGBLAST_DB_J:=}"; : "${IGBLAST_DATA:=}"; : "${IGBLAST_AUX:=}"

INHERITED_IGBLAST_DATA="${IGBLAST_DATA:-}"
mkdir -p "${DESTINY}"
DESTINY_ABS="$(cd "${DESTINY}" && pwd)"
# The final stage can be re-entered after a partial failure, so keep both the
# original handoff path and the already-moved result path in play.
FINAL_FASTA_PATH="${DESTINY}/split1/split2/collapsed/head/final/final.fasta"
RESULT_FASTA_PATH="${DESTINY}/${SAMPLE}.final.fasta"
IMGT_READY_FASTA_PATH="${DESTINY}/${SAMPLE}.imgt.ready.fasta"
QC_DIR="${DESTINY}/qc"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
QC_SCRIPT="${SCRIPT_DIR}/fasta_qc.R"
PLATE_QC_SCRIPT="${SCRIPT_DIR}/plate_stats.R"
PLATE_DEPTH_QC_SCRIPT="${SCRIPT_DIR}/plate_depth_qc.R"
PLATE_CLONOTYPE_QC_SCRIPT="${SCRIPT_DIR}/plate_clonotype_qc.R"
RUN_REPORT_SCRIPT="${SCRIPT_DIR}/run_report.R"
IGBLAST_SCRIPT="${SCRIPT_DIR}/run_igblast.sh"
CLONALITY_SCRIPT="${SCRIPT_DIR}/run_clonality.R"
RECEPTOR_SCRIPT="${SCRIPT_DIR}/detect_receptor.py"
PER_LOCUS_SCRIPT="${SCRIPT_DIR}/annotate_per_locus.py"
FOOTPRINT_SCRIPT="${SCRIPT_DIR}/trim_primer_footprint.py"
PLATE_QC_DIR="${QC_DIR}/plates"
IGBLAST_DIR="${DESTINY}/igblast"
CLONALITY_DIR="${DESTINY}/clonality"
RECEPTOR_DIR="${QC_DIR}/receptor"
ANNOTATION_DIR="${QC_DIR}/annotation"
ARCHIVE_NAME="$(basename "${OUTPUT}")"
ARCHIVE_PATH="${DESTINY_ABS}/${ARCHIVE_NAME}"
ARCHIVE_FORMAT="${ARCHIVE_FORMAT:-zip}"
PACKAGE_ROOT="${DESTINY}/.package_${SAMPLE}"
PACKAGE_DIR="${PACKAGE_ROOT}/${SAMPLE}"

normalize_optional() {
    local value="${1:-}"
    local lowered
    lowered="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    if [[ -z "${value}" || "${lowered}" == "na" || "${lowered}" == "nan" || "${lowered}" == "null" || "${lowered}" == "none" ]]; then
        return 0
    fi
    printf '%s' "${value}"
}

ARCHIVE_FORMAT="$(normalize_optional "${ARCHIVE_FORMAT}")"
IGBLAST_BIN="$(normalize_optional "${IGBLAST_BIN}")"
IGBLAST_SPECIES="$(normalize_optional "${IGBLAST_SPECIES}")"
IGBLAST_PANEL="$(normalize_optional "${IGBLAST_PANEL}")"
IGBLAST_DB_V="$(normalize_optional "${IGBLAST_DB_V}")"
IGBLAST_DB_D="$(normalize_optional "${IGBLAST_DB_D}")"
IGBLAST_DB_J="$(normalize_optional "${IGBLAST_DB_J}")"
IGBLAST_DATA="$(normalize_optional "${IGBLAST_DATA}")"
IGBLAST_AUX="$(normalize_optional "${IGBLAST_AUX}")"
IGBLAST_ORGANISM="$(normalize_optional "${IGBLAST_ORGANISM}")"

ARCHIVE_FORMAT="${ARCHIVE_FORMAT:-zip}"
EFFECTIVE_IGBLAST_DATA="${IGBLAST_DATA:-${INHERITED_IGBLAST_DATA}}"

# What goes INTO the archive. This is deliberately a different question from
# SHOULD_KEEP, which only decides whether the bulky intermediates survive on
# disk in the run directory (see the cleanup block near the end of this file).
#   light -> everyday outputs only: the imgt.ready FASTA, the untrimmed
#            assembled FASTA, every QC/annotation table, and the multi-page QC
#            PDFs.
#   full  -> everything light has, plus the per-plate PNG expansions of the
#            figures already in those PDFs.
# The untrimmed FASTA was full-only from 24 Aug to Oct 2026, on the mistaken
# premise that archives never carried it; RUN12-RUN14's did, users relied on it,
# and at ~3% of a light archive it is the only delivered file that shows what
# the primer trim removed.
# An older sample.json has no archive_scope field at all. Rather than silently
# changing what such a run produces in an arbitrary direction, fall back to the
# only related intent it does express: a submitter who asked to keep
# intermediates wanted the heavy bundle, everybody else gets the light one.
ARCHIVE_SCOPE="$(normalize_optional "${ARCHIVE_SCOPE:-}")"
ARCHIVE_SCOPE="$(printf '%s' "${ARCHIVE_SCOPE}" | tr '[:upper:]' '[:lower:]')"
if [[ "${ARCHIVE_SCOPE}" != "light" && "${ARCHIVE_SCOPE}" != "full" ]]; then
    SHOULD_KEEP_LOWER="${SHOULD_KEEP:-}"
    SHOULD_KEEP_LOWER="${SHOULD_KEEP_LOWER,,}"
    if [[ "${SHOULD_KEEP_LOWER}" == "true" ]]; then
        ARCHIVE_SCOPE="full"
    else
        ARCHIVE_SCOPE="light"
    fi
    echo "No archive scope supplied; derived archive_scope=${ARCHIVE_SCOPE} from keep_raw=${SHOULD_KEEP:-unset}."
fi
echo "Archive scope: ${ARCHIVE_SCOPE}"

mkdir -p "${QC_DIR}"
echo "Working directory: $(pwd)"
echo "Result directory: ${DESTINY}"

if [[ -f "${FINAL_FASTA_PATH}" ]]; then
    mv "${FINAL_FASTA_PATH}" "${RESULT_FASTA_PATH}"
elif [[ -f "${RESULT_FASTA_PATH}" ]]; then
    # A previous final-stage run may already have consumed final.fasta before
    # failing later during QC, packaging, or transfer.
    echo "Reusing existing final FASTA from a previous partial final-stage run: ${RESULT_FASTA_PATH}"
else
    echo "Expected assembled FASTA was not found at ${FINAL_FASTA_PATH}." >&2
    echo "The barcode splitting/assembly stage likely finished without producing final.fasta. Check logs/${SAMPLE}.genfasta.err.log." >&2
    exit 1
fi

# luka_light preprocessing already stripped the 5' handle inside
# split_fasta.sh, so running the same cutadapt again here can only destroy
# data: --discard-untrimmed would throw away every read whose handle is
# already gone. The Jan-2026 CMu/RUN10 run is the worked example -- 3,411
# assembled sequences went in, 417 came out, and exactly one of those was
# longer than 200 nt. That still holds for classic_luka and tiago_sulkseq.
# tiago_vsearch now clusters on trimmed copies but delivers the reads whole,
# so its luka_light FASTA is trimmed here like a standard one, plus the 3' tail
# luka_light always removed, and the primer footprint is cut per locus after
# annotation instead of a fixed 23 nt from every read.
PRETRIMMED="false"
if [[ "${PREPROCESS_METHOD,,}" == "luka_light" && "${PIPELINE_NAME:-}" != "tiago_vsearch" ]]; then
    PRETRIMMED="true"
fi
# tiago_nanopore's demultiplexer cuts through the 5' handle itself, so a
# cutadapt pass for it could only discard every sequence, and with no trim
# sequence set it stopped the run ("Adapter sequence is empty"). The primer
# after the handle is still on, so the per-locus footprint cut still applies.
HANDLE_REMOVED="false"
[[ "${PIPELINE_NAME:-}" == "tiago_nanopore" ]] && HANDLE_REMOVED="true"
EFFECTIVE_TRIM="${SHOULD_TRIM}"
FOOTPRINT_TRIM="${SHOULD_TRIM,,}"
if [[ "${SHOULD_TRIM,,}" == "true" && "${PRETRIMMED}" == "true" ]]; then
    EFFECTIVE_TRIM="false"
    FOOTPRINT_TRIM="false"
    echo "Primer trimming skipped: luka_light preprocessing already removed the 5' handle (${TRIMMING_SEQ}) during barcode splitting."
elif [[ "${SHOULD_TRIM,,}" == "true" && "${HANDLE_REMOVED}" == "true" ]]; then
    EFFECTIVE_TRIM="false"
    echo "5' handle trim skipped: the nanopore demultiplexer already removed it. The primer footprint is still cut per locus after annotation."
fi

if [[ "${EFFECTIVE_TRIM,,}" == "true" ]] && command -v Rscript >/dev/null 2>&1; then
    Rscript "${QC_SCRIPT}" "${RESULT_FASTA_PATH}" "${QC_DIR}" "pre_trim"
fi

# Trim fasta file if specified
if [[ "${EFFECTIVE_TRIM,,}" == "true" ]]; then
    # -O 15 raises cutadapt's minimum overlap from its default of 3. With the
    # default, a mismatched trim sequence still "matches" a few percent of
    # reads on 3 coincidental bases and the run delivers a handful of ~8 nt
    # fragments that look like a real result. A 15-base floor keeps every
    # genuine match (measured: 2007/2007 kappa, 2586/2589 heavy) and takes the
    # mismatched case to zero, where the check below can catch it.
    cutadapt -g "${TRIMMING_SEQ}" -e 0.2 -O 15 --discard-untrimmed -o "${IMGT_READY_FASTA_PATH}" "${RESULT_FASTA_PATH}"
    if [[ "${PREPROCESS_METHOD,,}" == "luka_light" ]]; then
        # The light-chain 3' tail and the plate barcode after it, as luka_light
        # has always removed them; settings as in its rescue step.
        TAIL_TMP="$(mktemp "${TMPDIR:-/tmp}/tbrc_tail.XXXXXX")"
        cutadapt -a "${LIGHT_TAIL_3P}" -e 0.2 -O 10 --minimum-length 1 -o "${TAIL_TMP}" "${IMGT_READY_FASTA_PATH}"
        cat "${TAIL_TMP}" > "${IMGT_READY_FASTA_PATH}"
    fi

    # A trim sequence that does not belong to this library does not error --
    # it silently deletes the library. Report the retention explicitly and
    # refuse to package an empty result.
    TRIM_IN="$(grep -c '^>' "${RESULT_FASTA_PATH}" || true)"
    TRIM_OUT="$(grep -c '^>' "${IMGT_READY_FASTA_PATH}" || true)"
    echo "Primer trim retention: ${TRIM_OUT}/${TRIM_IN} sequences kept."
    if [[ "${TRIM_OUT}" -eq 0 ]]; then
        echo "The trim sequence ${TRIMMING_SEQ} matched none of the ${TRIM_IN} assembled sequences, so the IMGT-ready FASTA would be empty." >&2
        echo "Check that the Trim sequence matches the primers used for this library, or set Trim primers to FALSE." >&2
        exit 1
    fi
    if [[ "${TRIM_IN}" -gt 0 ]] && (( TRIM_OUT * 2 < TRIM_IN )); then
        echo "WARNING: the trim sequence ${TRIMMING_SEQ} discarded $(( TRIM_IN - TRIM_OUT )) of ${TRIM_IN} sequences (>50%)." >&2
        echo "WARNING: this usually means the Trim sequence does not match the primers used for this barcode set." >&2
    fi
else
    cp "${RESULT_FASTA_PATH}" "${IMGT_READY_FASTA_PATH}"
fi

if command -v Rscript >/dev/null 2>&1; then
    if [[ "${EFFECTIVE_TRIM,,}" == "true" ]]; then
        Rscript "${QC_SCRIPT}" "${IMGT_READY_FASTA_PATH}" "${QC_DIR}" "post_trim"
    else
        Rscript "${QC_SCRIPT}" "${IMGT_READY_FASTA_PATH}" "${QC_DIR}" "final_output"
    fi
fi

# The two plate plots are diagnostics, and under `set -e` an unguarded failure here
# stops the script before the archive is built. ggplate's colour scale needs at
# least two wells with data, so a sparse plate (a small subsample, a failed
# library) made a run that had assembled and annotated deliver nothing. Skip the
# figure, say so, and keep going, as the clonotype QC below already does.
if [[ -f "${DESTINY}/BC.stats.txt" && -f "${DESTINY}/collapsed.stats.txt" ]] && command -v Rscript >/dev/null 2>&1; then
    if ! Rscript "${PLATE_QC_SCRIPT}" "${DESTINY}/" "${PLATE_QC_DIR}" "${SAMPLE}"; then
        echo "Plate QC failed and was skipped; the run continues without its plate figures." >&2
    fi
fi

if [[ -f "${IMGT_READY_FASTA_PATH}" ]] && command -v Rscript >/dev/null 2>&1; then
    if Rscript -e "quit(status = if (requireNamespace('ggplate', quietly = TRUE) && requireNamespace('clonality', quietly = TRUE)) 0 else 1)" >/dev/null 2>&1; then
        if ! Rscript "${PLATE_DEPTH_QC_SCRIPT}" "${IMGT_READY_FASTA_PATH}" "${PLATE_QC_DIR}"; then
            echo "Plate read-depth QC failed and was skipped; the run continues without it." >&2
        fi
    else
        echo "Plate read-depth QC skipped because ggplate and/or clonality are not installed." >&2
    fi
fi

# How many barcode plates this library actually used. Asking the submitter was
# the obvious way to get this and the wrong one -- it is a number nobody
# remembers at submission time, and it is sitting right here in the delivered
# headers once the run is done.
#
# The plate token is not always in the same place: barcode set B writes
# P05B04 (plate first), set B2 writes A01P01 (plate last), and the two-pass
# demultiplexer concatenates the barcode names in whichever order the set
# defines. So match P## anywhere in the well token rather than anchoring it.
PLATES_USED=0
if [[ -f "${IMGT_READY_FASTA_PATH}" ]]; then
    PLATES_USED="$(grep '^>' "${IMGT_READY_FASTA_PATH}" 2>/dev/null \
        | sed "s/^>//; s/^${SAMPLE}//" \
        | grep -oE 'P[0-9]{1,2}' \
        | sort -u | wc -l | tr -d ' ')"
    # A single-plate library carries no plate barcode at all, so its wells are
    # just row+column and the count above is legitimately zero.
    if [[ "${PLATES_USED}" == "0" ]] && grep -q '^>' "${IMGT_READY_FASTA_PATH}" 2>/dev/null; then
        PLATES_USED=1
    fi
fi
# The frontend reads this off the backend log the same way it reads input size.
echo "TBRC_PLATES_USED=${PLATES_USED}"

# ---------------------------------------------------------------- receptor type
# The panel is a choice made on a submission form; the sequences are the
# evidence. When they disagree nothing fails: IgBLAST aligns T-cell reads
# against an immunoglobulin germline database, finds nothing, and writes a table
# that is syntactically valid and biologically empty. Clonality then finds no
# clones in it. Both exit zero and the run is delivered. This is where that is
# caught, before either step runs.
RECEPTOR_VERDICT="not_run"
RECEPTOR_APPLIED="${IGBLAST_PANEL:-}"
SPECIES_VERDICT="not_run"
SPECIES_APPLIED="${IGBLAST_SPECIES:-}"
LOCUS_PLAN=""
if [[ "${RUN_IGBLAST,,}" == "true" && -f "${IMGT_READY_FASTA_PATH}" && -f "${RECEPTOR_SCRIPT}" ]] \
   && command -v python3 >/dev/null 2>&1; then
    RECEPTOR_MODE="$(printf '%s' "${RECEPTOR_CHECK}" | tr '[:upper:]' '[:lower:]')"
    case "${RECEPTOR_MODE}" in auto|warn|off) ;; *) RECEPTOR_MODE="auto" ;; esac

    # A pinned database is a more specific statement than a panel name. Changing
    # the panel underneath one would quietly ignore what the submitter pinned,
    # so in that case the check reports and does not act.
    if [[ -n "${IGBLAST_DB_V}" || -n "${IGBLAST_DB_J}" ]] && [[ "${RECEPTOR_MODE}" == "auto" ]]; then
        RECEPTOR_MODE="warn"
        echo "Receptor check: this submission pinned germline databases, so the panel will be reported but not changed."
    fi

    mkdir -p "${RECEPTOR_DIR}"
    if RECEPTOR_OUT="$(IGBLAST_DATA="${EFFECTIVE_IGBLAST_DATA}" python3 "${RECEPTOR_SCRIPT}" \
            --fasta "${IMGT_READY_FASTA_PATH}" \
            --igblast-script "${IGBLAST_SCRIPT}" \
            --species "${IGBLAST_SPECIES:-mouse}" \
            --organism "${IGBLAST_ORGANISM:-}" \
            --declared-panel "${IGBLAST_PANEL:-ig_all}" \
            --aux "${IGBLAST_AUX:-}" \
            --igblast-bin "${IGBLAST_BIN:-igblastn}" \
            --igblast-data "${EFFECTIVE_IGBLAST_DATA:-}" \
            --mode "${RECEPTOR_MODE}" \
            --report "${RECEPTOR_DIR}/locus_check.tsv" \
            --composition "${RECEPTOR_DIR}/locus_composition.tsv" \
            2>"${RECEPTOR_DIR}/detect.log")"; then
        RECEPTOR_APPLIED="$(printf '%s\n' "${RECEPTOR_OUT}" | sed -n 's/^APPLIED_PANEL=//p' | tail -1)"
        RECEPTOR_VERDICT="$(printf '%s\n' "${RECEPTOR_OUT}" | sed -n 's/^RECEPTOR_VERDICT=//p' | tail -1)"
        SPECIES_APPLIED="$(printf '%s\n' "${RECEPTOR_OUT}" | sed -n 's/^APPLIED_SPECIES=//p' | tail -1)"
        SPECIES_VERDICT="$(printf '%s\n' "${RECEPTOR_OUT}" | sed -n 's/^SPECIES_VERDICT=//p' | tail -1)"
        LOCUS_PLAN="$(printf '%s\n' "${RECEPTOR_OUT}" | sed -n 's/^LOCUS_PLAN=//p' | tail -1)"
        case "${RECEPTOR_VERDICT}" in
            agree)
                echo "Receptor check: sequences match the declared panel (${IGBLAST_PANEL:-ig_all})."
                ;;
            corrected)
                echo "Receptor check: the library is not what was declared. Panel changed from '${IGBLAST_PANEL:-ig_all}' to '${RECEPTOR_APPLIED}'." >&2
                IGBLAST_PANEL="${RECEPTOR_APPLIED}"
                ;;
            widened)
                # A single-locus panel on a plate that carries another locus.
                # RUN15 was submitted as igh and was 62% kappa; every kappa
                # sequence was forced onto an IGHV gene and delivered.
                echo "Receptor check: the plate carries loci the declared panel does not. Panel widened from '${IGBLAST_PANEL:-ig_all}' to '${RECEPTOR_APPLIED}'." >&2
                IGBLAST_PANEL="${RECEPTOR_APPLIED}"
                ;;
            narrow_panel)
                echo "Receptor check: the plate carries loci the declared panel does not, and the panel was left as submitted. Those sequences will be forced onto the wrong genes." >&2
                ;;
            mixed)
                if [[ -n "${RECEPTOR_APPLIED}" && "${RECEPTOR_APPLIED}" != "${IGBLAST_PANEL}" ]]; then
                    echo "Receptor check: both receptor classes are present. Panel changed from '${IGBLAST_PANEL}' to '${RECEPTOR_APPLIED}'." >&2
                    IGBLAST_PANEL="${RECEPTOR_APPLIED}"
                else
                    echo "Receptor check: both receptor classes are present; the declared panel already covers both." >&2
                fi
                ;;
            contradicted)
                echo "Receptor check: the sequences contradict the declared panel, and the panel was left as submitted. Annotation will be largely empty." >&2
                ;;
            skipped)
                : ;;
            *)
                echo "Receptor check: inconclusive; the panel was left as submitted. See qc/receptor/locus_check.tsv." >&2
                ;;
        esac
    else
        echo "Receptor check could not run; continuing with the declared panel. See qc/receptor/detect.log." >&2
        RECEPTOR_VERDICT="failed"
        SPECIES_VERDICT="failed"
    fi

    # The species, judged the same way. A plate whose heavy chains were human
    # went out annotated as mouse at 75% V identity with "completed" on the
    # report (RUN15_heavy). When another species wins decisively the check
    # names it and, in auto mode, the run is annotated against it. When the
    # loci need different sets -- one locus per species, or a knocked-in V on
    # the host's D and J -- the check returns a plan per locus and the
    # annotation step below follows it (annotate_per_locus.py). A mixed result
    # without a plan keeps the declared species and is reported.
    case "${SPECIES_VERDICT}" in
        agree)
            echo "Species check: sequences score best against ${IGBLAST_SPECIES:-mouse} germlines, as declared."
            ;;
        corrected)
            if [[ -n "${SPECIES_APPLIED}" && "${SPECIES_APPLIED}" != "${IGBLAST_SPECIES:-}" ]]; then
                echo "Species check: the library is not ${IGBLAST_SPECIES:-mouse}. Species changed to '${SPECIES_APPLIED}' for annotation." >&2
                IGBLAST_SPECIES="${SPECIES_APPLIED}"
                IGBLAST_ORGANISM="${SPECIES_APPLIED}"
                # An auxiliary file pinned for the declared species would give
                # the new one wrong CDR3 frames. Hand run_igblast.sh the
                # directory instead and let it pick the file by species.
                if [[ -n "${IGBLAST_AUX:-}" && -f "${IGBLAST_AUX}" ]]; then
                    IGBLAST_AUX="$(dirname "${IGBLAST_AUX}")"
                fi
            fi
            ;;
        contradicted)
            echo "Species check: the sequences score best against another species' germlines, and the species was left as submitted. V identity and gene calls will be wrong." >&2
            ;;
        per_locus)
            echo "Species check: the loci need different germline sets (${LOCUS_PLAN}); each locus is annotated against its own." >&2
            ;;
        mixed)
            echo "Species check: no single species accounts for the plate; the declared species was kept. See qc/receptor/locus_composition.tsv for which locus is which." >&2
            ;;
        *)
            : ;;
    esac
fi

# The report has never been able to say whether annotation ran. A user opening
# an archive with no igblast/ folder cannot tell whether it was skipped by
# choice, skipped because a database was missing, or silently broken -- and the
# report showed plate maps and read depth as though nothing were absent. Record
# the outcome here so the report can state it plainly.
IGBLAST_STATE="not_requested"
IGBLAST_REASON="IgBLAST was not requested for this run."
CLONALITY_STATE="not_requested"
CLONALITY_REASON="Clonality was not requested for this run."

# How well the delivered sequences match the germline set they were annotated
# with. A wrong species or a wrong strain does not fail: it delivers a table
# at 75% identity with every record "annotated" (RUN15_heavy, human heavy
# chains against mouse germlines). Real repertoires sit above 90% even with
# heavy somatic mutation, so a median this low is a wrong reference, not
# biology, and the report has to say so.
IDENTITY_FLOOR="${TBRC_IDENTITY_FLOOR:-85}"
IGBLAST_MEDIAN_IDENTITY=""
IGBLAST_IDENTITY_VERDICT=""
measure_v_identity() {
    local table="$1"
    [[ -s "${table}" ]] || return 0
    command -v python3 >/dev/null 2>&1 || return 0
    IGBLAST_MEDIAN_IDENTITY="$(python3 - "${table}" <<'PYEOF'
import csv, sys
vals = []
with open(sys.argv[1], newline="", errors="replace") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        raw = (row.get("v_identity") or "").strip()
        if not raw:
            continue
        try:
            vals.append(float(raw))
        except ValueError:
            pass
if vals:
    vals.sort()
    print("{:.1f}".format(vals[len(vals) // 2]))
PYEOF
)"
    if [[ -n "${IGBLAST_MEDIAN_IDENTITY}" ]]; then
        if python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) < float(sys.argv[2]) else 1)" "${IGBLAST_MEDIAN_IDENTITY}" "${IDENTITY_FLOOR}"; then
            IGBLAST_IDENTITY_VERDICT="low"
            echo "WARNING: median V identity is ${IGBLAST_MEDIAN_IDENTITY}%, below the ${IDENTITY_FLOOR}% floor. The germline set (${IGBLAST_SPECIES:-}/${IGBLAST_PANEL:-}) does not match these sequences; the gene calls in igblast/ are not trustworthy." >&2
        else
            IGBLAST_IDENTITY_VERDICT="ok"
        fi
        echo "Median V identity: ${IGBLAST_MEDIAN_IDENTITY}% (floor ${IDENTITY_FLOOR}%)."
    fi
}

write_annotation_status() {
    mkdir -p "${ANNOTATION_DIR}"
    # tiago_nanopore and tiago_sulkseq annotate inside the workflow, because
    # their own downstream steps need V, J and junction before this script runs.
    # For those the "Run IgBLAST" box is about the OPTIONAL second pass on the
    # delivered FASTA, so reporting "not annotated" would contradict the clone
    # counts on the same page.
    local pipeline_airr=""
    for candidate in \
        "${QC_DIR}/sulkseq/${SAMPLE}.airr.tsv" \
        "${QC_DIR}/nanopore/${SAMPLE}.long_consensus.airr.tsv"; do
        [[ -s "${candidate}" ]] && pipeline_airr="${candidate}" && break
    done
    if [[ -n "${pipeline_airr}" && "${IGBLAST_STATE}" == "not_requested" ]]; then
        IGBLAST_STATE="completed_by_pipeline"
        IGBLAST_REASON="Annotated by the ${PIPELINE_NAME:-pipeline} workflow, which needs V, J and junction for its own steps."
    fi
    {
        printf 'metric\tvalue\n'
        printf 'plates_used\t%s\n' "${PLATES_USED:-0}"
        printf 'igblast_state\t%s\n' "${IGBLAST_STATE}"
        printf 'igblast_reason\t%s\n' "${IGBLAST_REASON}"
        printf 'igblast_panel\t%s\n' "${IGBLAST_PANEL:-}"
        printf 'igblast_species\t%s\n' "${IGBLAST_SPECIES:-}"
        printf 'igblast_median_v_identity\t%s\n' "${IGBLAST_MEDIAN_IDENTITY:-}"
        printf 'igblast_identity_floor\t%s\n' "${IDENTITY_FLOOR}"
        printf 'igblast_identity_verdict\t%s\n' "${IGBLAST_IDENTITY_VERDICT:-}"
        printf 'species_verdict\t%s\n' "${SPECIES_VERDICT:-}"
        printf 'luka_light_heavy_records\t%s\n' "${LUKA_LIGHT_HEAVY:-}"
        printf 'clonality_state\t%s\n' "${CLONALITY_STATE}"
        printf 'clonality_reason\t%s\n' "${CLONALITY_REASON}"
        printf 'annotated_sequences\t%s\n' \
            "$(if [[ -f "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv" ]]; then
                   echo $(( $(wc -l < "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv") - 1 ))
               elif [[ -n "${pipeline_airr}" ]]; then
                   echo $(( $(wc -l < "${pipeline_airr}") - 1 ))
               else echo 0; fi)"
    } > "${ANNOTATION_DIR}/status.tsv"
}

if [[ "${RUN_IGBLAST,,}" == "true" ]]; then
    echo "IGBlast requested."
    echo "IGBlast submission fields:"
    echo "  species=${IGBLAST_SPECIES:-}"
    echo "  panel=${IGBLAST_PANEL:-}"
    echo "  bin=${IGBLAST_BIN:-}"
    echo "  organism=${IGBLAST_ORGANISM:-}"
    echo "  db_v=${IGBLAST_DB_V:-}"
    echo "  db_d=${IGBLAST_DB_D:-}"
    echo "  db_j=${IGBLAST_DB_J:-}"
    echo "  data=${EFFECTIVE_IGBLAST_DATA:-}"
    echo "  aux=${IGBLAST_AUX:-}"
    # Capture stderr separately so skipped-vs-failed IGBlast cases can be
    # handled without losing the detailed diagnostics. One file per call: it
    # used to be a fixed /tmp/tbrc_igblast.stderr, which two runs' final jobs
    # share whenever SLURM puts them on the same node, so one run could report
    # the other run's IgBLAST error as its own.
    IGBLAST_STDERR="$(mktemp "${TMPDIR:-/tmp}/tbrc_igblast.XXXXXX")"
    annotate() {
        if [[ -n "${LOCUS_PLAN}" && "${SPECIES_VERDICT}" == "per_locus" ]]; then
            IGBLAST_DATA="${EFFECTIVE_IGBLAST_DATA}" IGBLAST_BIN="${IGBLAST_BIN:-igblastn}" python3 "${PER_LOCUS_SCRIPT}" \
                --fasta "${IMGT_READY_FASTA_PATH}" --out-dir "${IGBLAST_DIR}" --sample "${SAMPLE}" \
                --plan "${LOCUS_PLAN}" --declared-species "${IGBLAST_SPECIES:-human}" \
                --panel "${IGBLAST_PANEL:-ig_all}" --organism "${IGBLAST_ORGANISM:-}" \
                --aux "${IGBLAST_AUX:-}" --igblast-script "${IGBLAST_SCRIPT}"
        else
            IGBLAST_DATA="${EFFECTIVE_IGBLAST_DATA}" IGBLAST_BIN="${IGBLAST_BIN:-igblastn}" "${IGBLAST_SCRIPT}" "${IMGT_READY_FASTA_PATH}" "${IGBLAST_DIR}" "${SAMPLE}" "${IGBLAST_SPECIES:-human}" "${IGBLAST_PANEL:-ig_all}" "${IGBLAST_DB_V:-}" "${IGBLAST_DB_D:-}" "${IGBLAST_DB_J:-}" "${IGBLAST_AUX:-}" "${IGBLAST_ORGANISM:-human}"
        fi
    }
    if annotate 2>"${IGBLAST_STDERR}"; then
        echo "IGBlast completed."
        IGBLAST_STATE="completed"
        if [[ -n "${LOCUS_PLAN}" && "${SPECIES_VERDICT}" == "per_locus" ]]; then
            IGBLAST_REASON="Annotated per locus with the ${IGBLAST_PANEL:-ig_all} panel (${LOCUS_PLAN}, as V species/D-J species); other loci as ${IGBLAST_SPECIES:-human}."
        else
            IGBLAST_REASON="Annotated with the ${IGBLAST_PANEL:-ig_all} germline panel for ${IGBLAST_SPECIES:-human}."
        fi
        measure_v_identity "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv"
        # The primer footprint, cut per locus now that the locus is known
        # (trim_primer_footprint.py says why), then annotated again so the
        # delivered table describes the delivered sequences.
        if [[ -n "${PRIMER_FOOTPRINT}" && "${FOOTPRINT_TRIM}" == "true" ]]; then
            FOOTPRINT_TMP="$(mktemp "${TMPDIR:-/tmp}/tbrc_footprint.XXXXXX")"
            FOOTPRINT_OUT="$(python3 "${FOOTPRINT_SCRIPT}" --fasta "${IMGT_READY_FASTA_PATH}" \
                --airr "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv" --footprint "${PRIMER_FOOTPRINT}" --out "${FOOTPRINT_TMP}")"
            echo "${FOOTPRINT_OUT}"
            if [[ "$(sed -n 's/^FOOTPRINT_CUT=//p' <<< "${FOOTPRINT_OUT}")" != "0" ]]; then
                cat "${FOOTPRINT_TMP}" > "${IMGT_READY_FASTA_PATH}"
                if ! annotate 2>>"${IGBLAST_STDERR}"; then
                    IGBLAST_STATUS=$?
                    cat "${IGBLAST_STDERR}" >&2 || true
                    IGBLAST_STATE="failed"
                    IGBLAST_REASON="IgBLAST failed on the primer-trimmed sequences (exit ${IGBLAST_STATUS})."
                    write_annotation_status
                    exit 1
                fi
                IGBLAST_REASON="${IGBLAST_REASON} Primer footprint cut per locus (${PRIMER_FOOTPRINT}) before the final annotation."
                measure_v_identity "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv"
            fi
        fi
        # classic_luka and tiago_sulkseq cut a fixed 23 nt after the 5' handle
        # before assembly. That is the F_106 Vk primer on a kappa read, but the
        # heavy primer is the handle itself, so on a heavy read the 23 nt are
        # template: FR1 nt 24-46. A mixed plate lets heavy reads through the
        # handle filter (RUN15_kappa: 424 records). Count them so the report can
        # say so. tiago_vsearch no longer has the problem (above).
        if [[ "${PRETRIMMED}" == "true" ]]; then
            LUKA_LIGHT_HEAVY="$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "v_call") c = i; next }
                c && $c ~ /^IGHV/ { n++ } END { print n + 0 }' "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv")"
        fi
    else
        IGBLAST_STATUS=$?
        cat "${IGBLAST_STDERR}" >&2 || true
        # 64 is run_igblast.sh's own "not set up on this host" code. It used to
        # be 2, which BLAST+ also returns for a database error, so a real abort
        # was filed under "skipped" and the run was delivered without its
        # annotation (RUN15_kappa, September 2026).
        if [[ "${IGBLAST_STATUS}" -eq 64 ]]; then
            echo "IGBlast skipped because server-side dependencies or database paths are not configured." >&2
            IGBLAST_STATE="skipped"
            IGBLAST_REASON="$(tail -1 "${IGBLAST_STDERR}" 2>/dev/null | tr -d '\t' || true)"
            IGBLAST_REASON="${IGBLAST_REASON:-The germline databases or the IgBLAST binary were not found on the server.}"
        else
            IGBLAST_STATE="failed"
            IGBLAST_REASON="IgBLAST exited with status ${IGBLAST_STATUS}."
            write_annotation_status
            exit "${IGBLAST_STATUS}"
        fi
    fi
fi

# Per-locus well verdicts, from the delivered records and their annotation
# (well_verdicts_by_locus.py says why the per-well verdict misleads on a mixed
# plate). Same thresholds as rank_wells.py; a failure only costs the table.
if [[ "${IGBLAST_STATE}" == "completed" && -f "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv" ]]; then
    mkdir -p "${QC_DIR}/assembly"
    python3 "${SCRIPT_DIR}/well_verdicts_by_locus.py" --fasta "${IMGT_READY_FASTA_PATH}" \
        --airr "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv" --sample "${SAMPLE}" \
        --out "${QC_DIR}/assembly/well_verdicts_by_locus.tsv" \
        || echo "Per-locus well verdicts failed and were skipped." >&2
fi

# CLONALITY_LIB points at an R library that holds a newer clonality than the
# conda env's own copy, so the package can be updated without rebuilding the
# production env (the 0.10 -> 0.15 move was validated this way on five finished
# runs). Unset, the env's copy is used exactly as before; unsetting it again is
# the whole rollback. It is prepended to R_LIBS for the clonality calls only.
run_with_clonality_lib() {
    if [[ -n "${CLONALITY_LIB:-}" ]]; then
        R_LIBS="${CLONALITY_LIB}${R_LIBS:+:${R_LIBS}}" "$@"
    else
        "$@"
    fi
}

if [[ "${RUN_CLONALITY,,}" == "true" ]]; then
    IGBLAST_OUTPUT="${IGBLAST_DIR}/${SAMPLE}.igblast.tsv"

    CLONALITY_STATE="skipped"
    CLONALITY_REASON="Clonality could not run."
    if [[ ! -f "${IGBLAST_OUTPUT}" ]]; then
        echo "Clonality skipped because IGBlast output was not found at ${IGBLAST_OUTPUT}." >&2
        CLONALITY_REASON="Clonality needs an IgBLAST table, and none was produced."
    elif command -v Rscript >/dev/null 2>&1; then
        if run_with_clonality_lib Rscript -e "quit(status = if (requireNamespace('clonality', quietly = TRUE)) 0 else 1)" >/dev/null 2>&1; then
            CLONALITY_VERSION="$(run_with_clonality_lib Rscript -e 'cat(as.character(packageVersion("clonality")))' 2>/dev/null || true)"
            if run_with_clonality_lib Rscript "${CLONALITY_SCRIPT}" "${IGBLAST_OUTPUT}" "${CLONALITY_DIR}" "${SAMPLE}" "${IGBLAST_PANEL:-ig_all}"; then
                echo "Clonality completed (clonality ${CLONALITY_VERSION:-unknown})."
                CLONALITY_STATE="completed"
                CLONALITY_REASON="Clones called from the ${IGBLAST_PANEL:-ig_all} annotation with clonality ${CLONALITY_VERSION:-unknown}."
            else
                CLONALITY_STATUS=$?
                if [[ "${CLONALITY_STATUS}" -eq 2 ]]; then
                    echo "Clonality skipped for the selected IGBlast receptor scope." >&2
                else
                    exit "${CLONALITY_STATUS}"
                fi
            fi
        else
            echo "Clonality skipped because the clonality R package is not installed." >&2
            CLONALITY_REASON="The clonality R package is not installed on the server."
        fi
    else
        echo "Clonality skipped because Rscript is not available." >&2
        CLONALITY_REASON="Rscript is not available on the server."
    fi
fi

write_annotation_status

CLONALITY_OUTPUT="${CLONALITY_DIR}/${SAMPLE}.clonality.tsv"
PLATE_CLONOTYPE_INPUT=""
if [[ -f "${CLONALITY_OUTPUT}" ]]; then
    PLATE_CLONOTYPE_INPUT="${CLONALITY_OUTPUT}"
elif [[ -f "${IGBLAST_DIR}/${SAMPLE}.igblast.tsv" ]]; then
    PLATE_CLONOTYPE_INPUT="${IGBLAST_DIR}/${SAMPLE}.igblast.tsv"
fi

if [[ -n "${PLATE_CLONOTYPE_INPUT}" ]] && command -v Rscript >/dev/null 2>&1; then
    if Rscript -e "quit(status = if (requireNamespace('ggplate', quietly = TRUE) && requireNamespace('ggplot2', quietly = TRUE)) 0 else 1)" >/dev/null 2>&1; then
        if ! Rscript "${PLATE_CLONOTYPE_QC_SCRIPT}" "${PLATE_CLONOTYPE_INPUT}" "${PLATE_QC_DIR}"; then
            echo "Plate clonotype contamination QC failed and was skipped." >&2
        fi
    else
        echo "Plate clonotype contamination QC skipped because ggplate and/or ggplot2 are not installed." >&2
    fi
fi

# One report per run, built from the QC tables written above. It has to run
# after every QC step and before the archive is staged, so the HTML and PDF are
# packaged with everything else rather than left behind on the cluster.
#
# Never fatal. A run whose report fails still delivers its sequences; losing the
# summary is an inconvenience, losing the archive is not acceptable.
if [[ -f "${RUN_REPORT_SCRIPT}" ]] && command -v Rscript >/dev/null 2>&1; then
    if Rscript -e "quit(status = if (requireNamespace('ggplot2', quietly = TRUE) && requireNamespace('dplyr', quietly = TRUE)) 0 else 1)" >/dev/null 2>&1; then
        echo "Building run report..."
        if Rscript "${RUN_REPORT_SCRIPT}" \
              -s "${SAMPLE}" \
              -d "${DESTINY}" \
              -f "${IMGT_READY_FASTA_PATH}" \
              -o "${DESTINY}/${SAMPLE}" \
              -p "${PIPELINE_NAME:-unknown}" \
              -b "${PRIMER_SET_NAME}" \
              -u "${SERVER_USER}" \
              -N "$([[ "${PIPELINE_NAME}" == "tiago_nanopore" ]] && echo 1 || echo 0)" \
              -S "${SCRIPT_DIR}"; then
            echo "Run report written."
        else
            echo "Run report failed; continuing without it." >&2
        fi
    else
        echo "Run report skipped: ggplot2/dplyr not available." >&2
    fi
fi

# Drop the per-plate PNG expansions from a staged package. Every one of them is
# one page of a sibling multi-page PDF (plate_stats.R and the panel QC scripts
# write both), so this removes bytes, not information -- and only ever when the
# PDF that carries the same figure is actually present in the package.
prune_redundant_plate_pngs() {
    local plates_dir="$1"
    local png base pdf
    [[ -d "${plates_dir}" ]] || return 0
    while IFS= read -r -d '' png; do
        base="${png##*/}"
        if [[ "${base}" == *.plate_*.png ]]; then
            pdf="${plates_dir}/${base%%.plate_*}.pdf"
        else
            pdf="${plates_dir}/${base%.png}.pdf"
        fi
        if [[ -f "${pdf}" ]]; then
            rm -f "${png}"
        fi
    done < <(find "${plates_dir}" -type f -name '*.png' -print0)
}

write_archive_manifest() {
    local manifest="${PACKAGE_DIR}/ARCHIVE_CONTENTS.txt"
    {
        echo "TBRC run ${SAMPLE}"
        echo "Archive scope: ${ARCHIVE_SCOPE}"
        echo
        echo "Always included:"
        printf '  %-28s %s\n' "${SAMPLE}.report.html" "the run report -- open this first; self-contained, needs no network"
        printf '  %-28s %s\n' "${SAMPLE}.report.pdf" "the same report as a printable multi-page PDF"
        if [[ -f "${DESTINY}/nanopore/${SAMPLE}.vdj_per_well.fasta" ]]; then
            printf '  %-28s %s\n' "${SAMPLE}.consensus_per_well.fasta" "one long consensus per well: full amplicon including the constant-region tail"
            printf '  %-28s %s\n' "${SAMPLE}.vdj_per_well.fasta" "one VDJ consensus per well, called more strictly for mutation analysis"
        fi
        printf '  %-28s %s\n' "${SAMPLE}.imgt.ready.fasta" "analysis-ready sequences (trimmed when primer trimming was on)"
        if [[ "${PRETRIMMED}" == "true" ]]; then
            printf '  %-28s %s\n' "${SAMPLE}.final.fasta" "assembled sequences; luka_light trims during assembly here, so these are trimmed too"
        elif [[ "${HANDLE_REMOVED}" == "true" ]]; then
            printf '  %-28s %s\n' "${SAMPLE}.final.fasta" "one consensus per well; demultiplexing removed barcodes and handle, a primer after the handle is still on"
        else
            printf '  %-28s %s\n' "${SAMPLE}.final.fasta" "assembled sequences before primer trimming: barcodes, handle and primers still on"
        fi
        printf '  %-28s %s\n' "qc/" "QC tables (.tsv) and multi-page figures (.pdf)"
        printf '  %-28s %s\n' "igblast/, clonality/" "annotation and clonality tables, when those steps ran"
        echo
        if [[ "${ARCHIVE_SCOPE}" == "full" ]]; then
            echo "Included because this run asked for the full archive:"
            printf '  %-28s %s\n' "qc/plates/*.plate_*.png" "one PNG per plate; the same pages are in the plate PDFs"
            printf '  %-28s %s\n' "qc/plates/*.panel.png" "PNG panel; the same figure is in the matching panel PDF"
        else
            echo "Left out of this light archive (re-submit with the full archive option to get them):"
            printf '  %-28s %s\n' "qc/plates/*.plate_*.png" "per-plate PNGs; every page is already in the plate PDFs"
            printf '  %-28s %s\n' "qc/plates/*.panel.png" "panel PNGs; the same figure is in the matching panel PDF"
        fi
        echo
        echo "Never archived, in either scope: the demultiplexed per-well FASTA tree"
        echo "(split1/), the raw assembled ${SAMPLE}.fasta and the trimmed FASTQ. Those stay in"
        echo "the run directory on the compute host, and only when the run kept intermediates"
        echo "(keep_raw=true)."
    } > "${manifest}"
}

rm -f "${OUTPUT}" "${ARCHIVE_PATH}"
rm -rf "${PACKAGE_ROOT}"
mkdir -p "${PACKAGE_DIR}"
(
    # Package outputs under a run-named folder so extracted archives do not mix
    # files from different runs in a user downloads folder.
    cp "${IMGT_READY_FASTA_PATH}" "${PACKAGE_DIR}/"
    # The run report, in both scopes. It is the one file a user should open
    # first, so it belongs in even the lightest archive; the two together
    # come to well under a megabyte.
    for report in "${DESTINY}/${SAMPLE}.report.html" "${DESTINY}/${SAMPLE}.report.pdf"; do
        [[ -f "${report}" ]] && cp "${report}" "${PACKAGE_DIR}/"
    done
    # A nanopore run delivers two sequences per well: the long consensus carries
    # the constant-region tail for locus-swap analysis, the VDJ consensus is
    # called strictly for mutation analysis. The manifest names both, so both
    # must actually be staged.
    for np in "${DESTINY}/nanopore/${SAMPLE}.consensus_per_well.fasta" \
              "${DESTINY}/nanopore/${SAMPLE}.vdj_per_well.fasta"; do
        [[ -f "${np}" ]] && cp "${np}" "${PACKAGE_DIR}/"
    done
    # Untrimmed, in both scopes (see the scope note near the top).
    [[ -f "${RESULT_FASTA_PATH}" ]] && cp "${RESULT_FASTA_PATH}" "${PACKAGE_DIR}/"
    cp -R "${QC_DIR}" "${PACKAGE_DIR}/"
    if [[ -d "${IGBLAST_DIR}" ]]; then
        cp -R "${IGBLAST_DIR}" "${PACKAGE_DIR}/"
    fi
    if [[ -d "${CLONALITY_DIR}" ]]; then
        cp -R "${CLONALITY_DIR}" "${PACKAGE_DIR}/"
    fi
    if [[ "${ARCHIVE_SCOPE}" != "full" ]]; then
        prune_redundant_plate_pngs "${PACKAGE_DIR}/qc/plates"
    fi
    write_archive_manifest
    cd "${PACKAGE_ROOT}"
    if [[ "${ARCHIVE_FORMAT}" == "tar.gz" ]]; then
        tar -czf "${ARCHIVE_PATH}" "${SAMPLE}"
    else
        zip -r "${ARCHIVE_PATH}" "${SAMPLE}"
    fi
)
rm -rf "${PACKAGE_ROOT}"

echo "Packaged output target: ${OUTPUT}"
echo "Packaged archive path: ${ARCHIVE_PATH}"
if [[ -f "${ARCHIVE_PATH}" ]]; then
    ls -lh "${ARCHIVE_PATH}"
else
    echo "Packaged archive not found at ${ARCHIVE_PATH}" >&2
    exit 1
fi

# Remove bulky intermediates before transfer so a later rsync failure does not
# leave the run directory in an unexpectedly half-cleaned state.
if [[ "${SHOULD_KEEP,,}" == "false" ]]; then
    rm -rf "${DESTINY}/split1/"
    rm -f "${DESTINY}/${SAMPLE}.fasta"
    rm -f "${DESTINY}/${SAMPLE}_S1_R1_001.fastq"
    rm -f "${DESTINY}/${SAMPLE}_S1_R2_001.fastq"
    rm -f "${DESTINY}/BC.stats.txt"
    rm -f "${DESTINY}/collapsed.stats.txt"
fi

# Transfer to the configured results destination. Push mode uses rsync over
# SSH; pull mode simply leaves the archive staged on the HPC.
if [[ "${SERVER_TRANSFER_MODE:-push}" == "pull" ]]; then
    echo "Transfer mode: pull"
    echo "Remote rsync skipped. Synology should pull: ${ARCHIVE_PATH}"
elif [[ -n "${SERVER_RESULTS_PATH:-}" && "${SERVER_RESULTS_PATH}" != "NA" ]]; then
    REMOTE_RESULTS_DIR="${SERVER_RESULTS_PATH%/}/${SERVER_USER}/"
    # The Snakefile only lets a plain path segment through, as 'none' when unset.
    if [[ -n "${RUN_GROUP}" && "${RUN_GROUP}" != "none" ]]; then
        REMOTE_RESULTS_DIR="${REMOTE_RESULTS_DIR}${RUN_GROUP}/"
    fi
    REMOTE_HOST="${REMOTE_RESULTS_DIR%%:*}"
    REMOTE_PATH="${REMOTE_RESULTS_DIR#*:}"
    SSH_ARGS=()
    RSYNC_SSH=()

    if [[ "${REMOTE_HOST}" == "${REMOTE_PATH}" ]]; then
        echo "Results destination is not in host:path format: ${SERVER_RESULTS_PATH}" >&2
        exit 1
    fi

    if [[ -n "${SERVER_SSH_PORT:-}" && "${SERVER_SSH_PORT}" != "NA" ]]; then
        SSH_ARGS=(-p "${SERVER_SSH_PORT}")
        RSYNC_SSH=(-e "ssh -p ${SERVER_SSH_PORT}")
    fi

    echo "Transfer destination: ${REMOTE_RESULTS_DIR}"
    echo "Transfer ssh port: ${SERVER_SSH_PORT:-default}"
    ssh "${SSH_ARGS[@]}" "${REMOTE_HOST}" "mkdir -p '${REMOTE_PATH%/}'"
    rsync -avP "${RSYNC_SSH[@]}" "${ARCHIVE_PATH}" "${REMOTE_RESULTS_DIR}"
    REMOTE_FILE_PATH="${REMOTE_PATH%/}/$(basename "${ARCHIVE_PATH}")"
    echo "Remote file check: ${REMOTE_HOST}:${REMOTE_FILE_PATH}"
    ssh "${SSH_ARGS[@]}" "${REMOTE_HOST}" "ls -lh '${REMOTE_FILE_PATH}'"
fi

echo "${SERVER_STORAGE} ok."

send_notification_email() {
    local recipient="${1:-}"
    local subject="$2"
    local body="$3"

    if [[ -z "${recipient}" ]]; then
        echo "Notification email skipped: no recipient provided." >&2
        return 0
    fi
    if [[ ! "${recipient}" =~ ^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}$ ]]; then
        # The app checks this too. Re-checked here because a script that hands a
        # value to the mail command should not rely on its caller being correct.
        echo "Notification email skipped: recipient is not a valid address." >&2
        return 0
    fi

    echo "Notification recipient: ${recipient}"

    if command -v mail >/dev/null 2>&1; then
        if mail -s "${subject}" "${recipient}" <<< "${body}"; then
            echo "Notification email accepted by mail command."
            return 0
        fi
        echo "Notification email failed through mail command." >&2
        return 1
    fi

    if command -v sendmail >/dev/null 2>&1; then
        if printf 'To: %s\nSubject: %s\n\n%s\n' "${recipient}" "${subject}" "${body}" | sendmail -t; then
            echo "Notification email accepted by sendmail."
            return 0
        fi
        echo "Notification email failed through sendmail." >&2
        return 1
    fi

    echo "Notification email skipped: neither mail nor sendmail is available on this host." >&2
    return 0
}

EMAIL_SUBJECT="Your Clonality Run ${SAMPLE} is completed"
EMAIL_BODY="Your clonality run is completed successfully and it was uploaded to the folder: pipelines/tbrc/results/ of ${SERVER_STORAGE}."

if ! send_notification_email "${NOTIFICATION_EMAIL}" "${EMAIL_SUBJECT}" "${EMAIL_BODY}"; then
    echo "Notification delivery did not complete successfully." >&2
fi
