#!/bin/bash
# Check a run before it is submitted, not after it fails.
#
# Every scenario below is one somebody has actually hit, or one the pipeline
# handles so quietly that the user would be left guessing. The pattern in the
# historical failures is not that the pipeline lacks guards -- it has some -- but
# that they fire after dispatch, so the person waits, gets a job failure, and has
# to read a cluster log to learn the folder name was misspelled.
#
# Output is one finding per line:
#     SEVERITY|CODE|MESSAGE|HINT
# SEVERITY is FAIL (submission should be blocked), WARN (proceed, but the user
# should know) or OK. The caller decides what to do with them; this script only
# reports.
set -uo pipefail
source "$(cd "$(dirname "$0")" && pwd)/common.sh"

STORAGE_SERVER_NAME=""; INPUT_FOLDER=""; PIPELINE=""; PRIMER_SET=""; TRIM_SEQ=""; RUN_NAME=""; USER_ID=""
while getopts "s:i:p:b:t:n:u:" flag; do
  case "${flag}" in
    s) STORAGE_SERVER_NAME="${OPTARG}" ;;
    i) INPUT_FOLDER="${OPTARG}" ;;
    p) PIPELINE="${OPTARG}" ;;
    b) PRIMER_SET="${OPTARG}" ;;
    t) TRIM_SEQ="${OPTARG}" ;;
    n) RUN_NAME="${OPTARG}" ;;
    u) USER_ID="${OPTARG}" ;;
    *) ;;
  esac
done

emit() { printf '%s|%s|%s|%s\n' "$1" "$2" "$3" "$4"; }

# ---------------------------------------------------------------- local checks

if [[ -z "${INPUT_FOLDER}" ]]; then
  emit FAIL no_folder "No source folder given." \
       "Type the folder name exactly as it appears on the storage server."
  exit 0
fi

# validate_path_segment exits rather than returns (goServer.bash wants that), so
# call it in a subshell: called directly it ended this script before the FAIL
# line below was printed, and the app saw an empty preflight instead of the
# reason. The folder name goes into remote commands further down, so a bad
# one also stops the checks here.
if ! ( validate_path_segment "${INPUT_FOLDER}" ) 2>/dev/null; then
  emit FAIL bad_folder_name "Folder name contains characters that are not allowed: ${INPUT_FOLDER}" \
       "Use only letters, numbers, dot, dash and underscore. Spaces and slashes will not work."
  exit 0
fi

if [[ -n "${RUN_NAME}" ]] && ! ( validate_path_segment "${RUN_NAME}" ) 2>/dev/null; then
  emit FAIL bad_run_name "Run name contains characters that are not allowed: ${RUN_NAME}" \
       "Use only letters, numbers, dot, dash and underscore."
fi

RAWDATA_SERVER="$(json_get_storage_field "${STORAGE_SERVER_NAME}" rawdata_path 2>/dev/null)"
SSH_PORT="$(json_get_storage_field "${STORAGE_SERVER_NAME}" ssh_port 2>/dev/null)"
if [[ -z "${RAWDATA_SERVER}" || "${RAWDATA_SERVER}" == "null" ]]; then
  emit FAIL unknown_server "Storage server not recognised: ${STORAGE_SERVER_NAME:-<none>}" \
       "Pick a server from the list. If it should be there, the server configuration needs updating."
  exit 0
fi

HOST="${RAWDATA_SERVER%%:*}"
BASE="${RAWDATA_SERVER#*:}"
REMOTE_PATH="${BASE%/}/${INPUT_FOLDER}"
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=10)
[[ -n "${SSH_PORT}" && "${SSH_PORT}" != "null" ]] && SSH_OPTS+=(-p "${SSH_PORT}")

# ---------------------------------------------------------------- reachability

if ! ssh -n "${SSH_OPTS[@]}" "${HOST}" true >/dev/null 2>&1; then
  emit FAIL server_unreachable "Cannot reach the storage server ${HOST}." \
       "It may be off, or you may be off the network. Nothing is wrong with your run settings."
  exit 0
fi

# ------------------------------------------------------- the folder and its files
# One remote call gathers everything: a round trip per check would be slow enough
# that people stop using the button.
PROBE=$(ssh -n "${SSH_OPTS[@]}" "${HOST}" "
  p=$(shq "${REMOTE_PATH}")
  if [ ! -d \"\$p\" ]; then echo 'MISSING'; exit 0; fi
  echo 'EXISTS'
  echo \"ALL:\$(find \"\$p\" -maxdepth 1 -type f ! -name '.*' 2>/dev/null | wc -l)\"
  echo \"FQ:\$(find \"\$p\" -maxdepth 1 -type f \\( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \\) 2>/dev/null | wc -l)\"
  echo \"EMPTY:\$(find \"\$p\" -maxdepth 1 -type f -size 0 2>/dev/null | wc -l)\"
  echo \"BYTES:\$(du -sb \"\$p\" 2>/dev/null | cut -f1)\"
  find \"\$p\" -maxdepth 1 -type f \\( -name '*.fastq*' -o -name '*.fq*' \\) -printf 'F:%f:%s\n' 2>/dev/null | head -20
  find \"\$p\" -maxdepth 1 -type f ! -name '*.fastq*' ! -name '*.fq*' ! -name '.*' -printf 'O:%f\n' 2>/dev/null | head -6
" 2>/dev/null)

if [[ -z "${PROBE}" ]] || grep -q '^MISSING$' <<<"${PROBE}"; then
  emit FAIL folder_missing "No folder called '${INPUT_FOLDER}' on ${STORAGE_SERVER_NAME}." \
       "Check the spelling and the capitalisation. Folder names are case sensitive."
  exit 0
fi

ALL=$(grep -oP '(?<=^ALL:)\d+' <<<"${PROBE}" | head -1); ALL=${ALL:-0}
FQ=$(grep -oP '(?<=^FQ:)\d+'  <<<"${PROBE}" | head -1); FQ=${FQ:-0}
EMPTY=$(grep -oP '(?<=^EMPTY:)\d+' <<<"${PROBE}" | head -1); EMPTY=${EMPTY:-0}
BYTES=$(grep -oP '(?<=^BYTES:)\d+' <<<"${PROBE}" | head -1); BYTES=${BYTES:-0}

if [[ "${ALL}" -eq 0 ]]; then
  emit FAIL folder_empty "The folder '${INPUT_FOLDER}' exists but has no files in it." \
       "The upload may still be running, or the data may be in a subfolder."
  exit 0
fi

if [[ "${FQ}" -eq 0 ]]; then
  others=$(grep '^O:' <<<"${PROBE}" | sed 's/^O://' | paste -sd', ' - | cut -c1-90)
  emit FAIL no_fastq "The folder has ${ALL} file(s) but none of them are FASTQ." \
       "Found instead: ${others:-unknown}. The pipeline needs .fastq or .fastq.gz files."
  exit 0
fi

if [[ "${EMPTY}" -gt 0 ]]; then
  emit FAIL empty_files "${EMPTY} file(s) in the folder are 0 bytes." \
       "A zero-byte file usually means an interrupted copy. Re-upload before submitting."
fi

# ------------------------------------------------------- pairing and platform
# The Snakefile requires exactly one R1 and one R2 for the Illumina pipelines,
# and a single FASTQ for nanopore. Getting this wrong fails deep inside
# snakemake with a ValueError about file counts, which tells the user nothing.
R1=$(grep -cE '^F:.*([_.-]R1[_.-])' <<<"${PROBE}" || true)
R2=$(grep -cE '^F:.*([_.-]R2[_.-])' <<<"${PROBE}" || true)

if [[ "${PIPELINE}" == "tiago_nanopore" ]]; then
  if [[ "${R1}" -gt 0 && "${R2}" -gt 0 ]]; then
    emit WARN nanopore_paired "These look like paired Illumina reads (R1 and R2), but the nanopore pipeline is selected." \
         "Nanopore runs are a single FASTQ. Did you mean tiago_vsearch?"
  fi
else
  if [[ "${FQ}" -eq 1 ]]; then
    emit FAIL single_file "Only one FASTQ found, but this pipeline needs a read pair." \
         "Illumina runs need both R1 and R2. For a single long-read FASTQ, choose tiago_nanopore."
  elif [[ "${R1}" -ne 1 || "${R2}" -ne 1 ]]; then
    emit FAIL pairing "Cannot identify exactly one R1 and one R2 (found R1=${R1}, R2=${R2} among ${FQ} FASTQ files)." \
         "File names must contain _R1_ and _R2_. Extra FASTQs in the folder also cause this."
  fi
fi

# ------------------------------------------------------------ file integrity
# A truncated gzip is the classic interrupted-transfer failure. It survives every
# name-based check and only surfaces when a tool tries to read it.
#
# There is no cheap version of this test. Truncation can only be found by
# decompressing to the end of the stream: reading the tail alone reports "not in
# gzip format" for any mid-stream fragment, healthy or not, so a tail-based
# check silently passes everything. That was tried here and it did exactly that.
#
# So the test is the real one, under a wall-clock budget, and the report says
# what was actually verified. The first version had no budget at all: on a 19 GB
# folder it held a submission for six minutes, and because Shiny is single
# threaded it froze the app for every connected user, not just the submitter.
INTEGRITY_BUDGET="${TBRC_PREFLIGHT_GZIP_SECONDS:-20}"
INTEGRITY_PER_FILE="${TBRC_PREFLIGHT_GZIP_FILE_SECONDS:-8}"
INTEGRITY_OUT=$(ssh -n "${SSH_OPTS[@]}" "${HOST}" "
  cd $(shq "${REMOTE_PATH}") 2>/dev/null || exit 0
  have_timeout=0; command -v timeout >/dev/null 2>&1 && have_timeout=1
  deadline=\$(( \$(date +%s) + ${INTEGRITY_BUDGET} ))
  ok=0; skipped=0
  for f in *.gz; do
    [ -e \"\$f\" ] || continue
    # The overall budget stops us starting another file; the per-file timeout
    # stops one huge file consuming the whole budget by itself. A 10 GB FASTQ
    # takes minutes to decompress, so without the second bound the first bound
    # never gets a chance to apply.
    if [ \$(date +%s) -ge \$deadline ]; then skipped=\$(( skipped + 1 )); continue; fi
    if [ \$have_timeout -eq 1 ]; then
      timeout ${INTEGRITY_PER_FILE} gzip -t -- \"\$f\" 2>/dev/null
      rc=\$?
    else
      gzip -t -- \"\$f\" 2>/dev/null; rc=\$?
    fi
    if [ \$rc -eq 0 ]; then ok=\$(( ok + 1 ))
    elif [ \$rc -eq 124 ]; then skipped=\$(( skipped + 1 ))
    else echo \"BAD \$f\"
    fi
  done
  echo \"OK \$ok\"; echo \"SKIPPED \$skipped\"
" 2>/dev/null)

BAD=$(printf '%s\n' "${INTEGRITY_OUT}" | sed -n 's/^BAD //p' | head -5)
VERIFIED=$(printf '%s\n' "${INTEGRITY_OUT}" | sed -n 's/^OK //p' | tail -1)
SKIPPED=$(printf '%s\n' "${INTEGRITY_OUT}" | sed -n 's/^SKIPPED //p' | tail -1)

if [[ -n "${BAD}" ]]; then
  emit FAIL truncated "Some files are corrupt or truncated: $(paste -sd', ' <<<"${BAD}" | cut -c1-80)" \
       "This is what an interrupted upload looks like. Copy them again and re-check."
elif [[ "${SKIPPED:-0}" -gt 0 ]]; then
  emit WARN integrity_unchecked \
       "Verified ${VERIFIED:-0} file(s) as intact; ${SKIPPED} were too large to finish checking in time." \
       "The run was not blocked. If it fails while reading a FASTQ, an interrupted upload is the first thing to suspect."
fi

# ---------------------------------------------------------------- size sanity
if [[ "${BYTES}" -gt 0 && "${BYTES}" -lt 10485760 ]]; then
  emit WARN tiny_input "The folder holds only $(( BYTES / 1024 )) KB of data." \
       "That is far smaller than a normal run and usually means a partial upload."
fi

# ------------------------------------------------- would this overwrite results?
if [[ -n "${RUN_NAME}" && -n "${USER_ID}" ]]; then
  RES_ROOT="$(json_get_cluster_field pipeline_root 2>/dev/null)"
  DTN="$(json_get_cluster_field dtn_host 2>/dev/null)"
  if [[ -n "${RES_ROOT}" && -n "${DTN}" ]]; then
    if ssh -n -o BatchMode=yes -o ConnectTimeout=10 "${DTN}" \
        "test -d $(shq "${RES_ROOT}/results/${USER_ID}/${INPUT_FOLDER}/${RUN_NAME}")" 2>/dev/null; then
      emit WARN overwrite "A run called '${RUN_NAME}' already exists for ${USER_ID}/${INPUT_FOLDER}." \
           "Submitting again will replace those results. Use a different run name to keep both."
    fi
  fi
fi

# ------------------------------------------------ does the trim sequence match?
# The most damaging failure in this pipeline's history was silent: six T-cell
# runs were trimmed with the B-cell handle. Nothing errored -- cutadapt simply
# discarded almost everything, and the delivered archive looked like a poor
# library rather than a wrong parameter. Checking the handle against the actual
# reads costs one command and would have caught all six.
if [[ -n "${TRIM_SEQ}" && "${#TRIM_SEQ}" -ge 12 ]]; then
  PROBE_SEQ="${TRIM_SEQ:0:12}"
  ONE=$(grep -m1 -oP '(?<=^F:)[^:]+' <<<"${PROBE}")
  if [[ -n "${ONE}" ]]; then
    HITS=$(ssh -n "${SSH_OPTS[@]}" "${HOST}" "
      f=$(shq "${REMOTE_PATH}/${ONE}")
      case \"\$f\" in *.gz) zcat -- \"\$f\" 2>/dev/null;; *) cat -- \"\$f\" 2>/dev/null;; esac \
        | head -80000 | awk 'NR%4==2' | grep -c -- $(shq "${PROBE_SEQ}") || true
    " 2>/dev/null | tr -d '[:space:]')
    HITS=${HITS:-0}
    if [[ "${HITS}" -eq 0 ]]; then
      emit FAIL trim_mismatch "The trim sequence does not appear in the reads at all (checked 20,000 reads)." \
           "This is almost always the wrong primer handle for this library -- for example a B-cell handle on a T-cell run. Trimming would discard nearly everything."
    elif [[ "${HITS}" -lt 200 ]]; then
      emit WARN trim_rare "The trim sequence appears in only ${HITS} of 20,000 reads." \
           "Check the handle matches this library. A correct handle usually matches thousands."
    fi
  fi
fi

emit OK checked "Checked ${FQ} FASTQ file(s), $(( BYTES / 1048576 )) MB, on ${STORAGE_SERVER_NAME}." ""
