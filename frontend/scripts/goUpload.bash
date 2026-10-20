#!/bin/bash

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

require_tool rsync
require_tool ssh

USER_ID=""
DATA_PATH_1=""
DATA_PATH_2=""
FILE_NAME_1=""
FILE_NAME_2=""

while getopts "u:d:e:n:m:" flag; do
    case "${flag}" in
        u) USER_ID="${OPTARG}" ;;
        d) DATA_PATH_1="${OPTARG}" ;;
        e) DATA_PATH_2="${OPTARG}" ;;
        n) FILE_NAME_1="${OPTARG}" ;;
        m) FILE_NAME_2="${OPTARG}" ;;
        *) echo "Invalid option" >&2; exit 1 ;;
    esac
done

if [[ -z "${USER_ID}" || -z "${DATA_PATH_1}" || -z "${DATA_PATH_2}" || -z "${FILE_NAME_1}" || -z "${FILE_NAME_2}" ]]; then
    echo "Missing required arguments." >&2
    exit 1
fi

# FILE_NAME_1 and FILE_NAME_2 are the names of the uploaded files, chosen
# entirely by whoever uploaded them. They reach a directory path, a cp target and
# -- via DATA_TABLE and INPUTPATH -- a command run over ssh on the login node. A
# name containing a single quote used to close the quoting in that command and
# run the remainder as user on the cluster.
#
# The app checks these too. This check is here as well because a script that
# executes remote commands should not depend on its only caller being correct.
validate_upload_filename() {
    local name="$1"
    if [[ "${name}" != "$(basename -- "${name}")" ]]; then
        echo "Rejected file name (path separator): ${name}" >&2; exit 1
    fi
    if [[ ! "${name}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Rejected file name (illegal character): ${name}" >&2; exit 1
    fi
    if [[ "${name}" == .* ]]; then
        echo "Rejected file name (leading dot): ${name}" >&2; exit 1
    fi
    if [[ ! "${name}" =~ [.](fastq|fq)([.]gz)?$ ]]; then
        echo "Rejected file name (not a FASTQ): ${name}" >&2; exit 1
    fi
}
validate_upload_filename "${FILE_NAME_1}"
validate_upload_filename "${FILE_NAME_2}"
validate_path_segment "${USER_ID}" || { echo "Rejected user id: ${USER_ID}" >&2; exit 1; }

DTN_HOST="$(json_get_cluster_field dtn_host)"
SSH_HOST="$(json_get_cluster_field ssh_host)"
PIPELINE_ROOT="$(json_get_cluster_field pipeline_root)"
SNAKEMAKE_ENV="$(json_get_cluster_field snakemake_env)"
PATH_PREFIX="$(json_get_cluster_field path_prefix)"

RUN_DIR="${RUNS_DIR}/${USER_ID}/${FILE_NAME_1%.*}"
mkdir -p "${RUN_DIR}"
# The app now stages uploads here itself, so the source may already be the
# target; cp refuses to copy a file onto itself.
[[ "${DATA_PATH_1}" -ef "${RUN_DIR}/${FILE_NAME_1}" ]] || cp "${DATA_PATH_1}" "${RUN_DIR}/${FILE_NAME_1}"
[[ "${DATA_PATH_2}" -ef "${RUN_DIR}/${FILE_NAME_2}" ]] || cp "${DATA_PATH_2}" "${RUN_DIR}/${FILE_NAME_2}"
INPUT_BYTES="$(
    wc -c < "${RUN_DIR}/${FILE_NAME_1}" 2>/dev/null
    wc -c < "${RUN_DIR}/${FILE_NAME_2}" 2>/dev/null
)"
INPUT_BYTES="$(printf '%s\n' "${INPUT_BYTES}" | awk '{sum+=$1} END {print sum+0}')"
echo "TBRC_INPUT_BYTES=${INPUT_BYTES:-0}"

rsync -a -P "${RUNS_DIR}/${USER_ID}" "${DTN_HOST}:${PIPELINE_ROOT}/data/runs/"

DATA_TABLE="data/runs/${USER_ID}/${FILE_NAME_1%.*}/sample.json"
INPUTPATH="data/runs/${USER_ID}/${FILE_NAME_1%.*}"

# Every value is quoted with shq() (printf %q) rather than wrapped in literal
# single quotes, so nothing a value contains can end its own quoting.
# Own working directory per run, as in goServer.bash, so concurrent
# submissions do not refuse each other; started detached, so the run survives
# this connection.
dispatch_detached "${SSH_HOST}" "${PIPELINE_ROOT}" "${PATH_PREFIX}" "${SNAKEMAKE_ENV}" \
  "${USER_ID}-${FILE_NAME_1%.*}" "${DATA_TABLE}" "${INPUTPATH}"
