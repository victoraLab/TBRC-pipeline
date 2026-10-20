#!/bin/bash

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

require_tool jq
require_tool rsync
require_tool ssh

STORAGE_SERVER_NAME=""
USER_ID=""
INPUT_FOLDER=""
RUN_NAME=""
LOCAL_RUN=0

while getopts ":s:u:f:n:l" opt; do
  case "${opt}" in
    s) STORAGE_SERVER_NAME="${OPTARG}" ;;
    u) USER_ID="${OPTARG}" ;;
    # Trimmed here as well as in the frontend. A script that builds a remote
    # path from this value should not depend on its only caller having cleaned
    # it: whitespace produced a "Folder not found" for a folder that existed.
    f) INPUT_FOLDER="$(printf '%s' "${OPTARG}" | tr -d '\r\n' | sed 's#^[[:space:]/]*##; s#[[:space:]/]*$##')" ;;
    n) RUN_NAME="${OPTARG}" ;;
    l) LOCAL_RUN=1 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    :) echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [[ -z "${STORAGE_SERVER_NAME}" || -z "${USER_ID}" || -z "${INPUT_FOLDER}" || -z "${RUN_NAME}" ]]; then
  echo "Missing required arguments." >&2
  exit 1
fi

# These three values originate in a web form and are interpolated into shell
# commands executed on the DTN and on the storage server. Validate before use.
validate_path_segment "${USER_ID}" "User name"
validate_path_segment "${INPUT_FOLDER}" "Folder in server"
validate_path_segment "${RUN_NAME}" "Run name"

RAWDATA_SERVER="$(json_get_storage_field "${STORAGE_SERVER_NAME}" rawdata_path)"
STORAGE_SERVER_SSH_PORT="$(json_get_storage_field "${STORAGE_SERVER_NAME}" ssh_port)"
DTN_HOST="$(json_get_cluster_field dtn_host)"
SSH_HOST="$(json_get_cluster_field ssh_host)"
PIPELINE_ROOT="$(json_get_cluster_field pipeline_root)"
SNAKEMAKE_ENV="$(json_get_cluster_field snakemake_env)"
PATH_PREFIX="$(json_get_cluster_field path_prefix)"

if [[ -z "${RAWDATA_SERVER}" || "${RAWDATA_SERVER}" == "null" ]]; then
  echo "Server not recognized or missing configuration: ${STORAGE_SERVER_NAME}" >&2
  exit 1
fi

REMOTE_STORAGE_HOST="${RAWDATA_SERVER%%:*}"
REMOTE_STORAGE_BASE_PATH="${RAWDATA_SERVER#*:}"
REMOTE_INPUT_PATH="${REMOTE_STORAGE_BASE_PATH%/}/${INPUT_FOLDER}"

if [[ "${REMOTE_STORAGE_HOST}" == "${REMOTE_STORAGE_BASE_PATH}" ]]; then
  echo "Storage server path is not in host:path format: ${RAWDATA_SERVER}" >&2
  exit 1
fi

# sample.json is written locally by the app, then copied to the HPC run folder
# before Snakemake starts.
DATA_TABLE="data/runs/${USER_ID}/${RUN_NAME}/sample.json"
INPUTPATH="data/runs/${USER_ID}/${RUN_NAME}/"
# Each run gets its own working directory on the backend
# (backend/run_isolated.sh), so a second submission is no longer refused by the
# first one's Snakemake lock. User plus run name, because two users may pick
# the same run name; resubmitting the same run while it is live is still
# refused, by that directory's own driver guard.
RUN_KEY="${USER_ID}-${RUN_NAME}"

if [[ "${LOCAL_RUN}" == "1" ]]; then
  # The app gates this toggle on capability, but guard here too: a deploy that
  # ships only frontend/ has no sibling backend/, and the bare "cd" failure this
  # replaces was the single largest cause of run failures.
  if [[ ! -d "${BACKEND_DIR}" ]]; then
    echo "Cannot run locally: no backend directory at ${BACKEND_DIR}." >&2
    echo "This host received only the frontend. Submit to the HPC cluster instead." >&2
    exit 1
  fi

  for tool in snakemake pandaseq cutadapt fastx_trimmer; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
      echo "Cannot run locally: ${tool} is not installed on this host." >&2
      echo "Submit to the HPC cluster instead." >&2
      exit 1
    fi
  done

  cd "${BACKEND_DIR}"
  export DATA_TABLE="${DATA_TABLE}"
  export INPUTPATH="${INPUTPATH}"

  # PATH_PREFIX and the named Conda environment describe the RU cluster. In a
  # container the tools are already on PATH in Conda's base environment, and
  # forcing the cluster's PATH would hide them. Apply each only where it exists.
  if [[ -n "${PATH_PREFIX}" && -d "${PATH_PREFIX%%:*}" ]]; then
    export PATH="${PATH_PREFIX}"
  fi

  if [[ -n "${SNAKEMAKE_ENV}" && "${SNAKEMAKE_ENV}" != "base" ]] \
     && command -v conda >/dev/null 2>&1 \
     && conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${SNAKEMAKE_ENV}"; then
    conda activate "${SNAKEMAKE_ENV}"
  fi

  ./run_isolated.sh "${RUN_KEY}" "${DATA_TABLE}" "${INPUTPATH}"
  exit 0
fi

RSYNC_SSH_OPTION="$(rsync_transport_args "${STORAGE_SERVER_SSH_PORT}")"
STORAGE_SSH_CHECK=""
if [[ -n "${STORAGE_SERVER_SSH_PORT}" && "${STORAGE_SERVER_SSH_PORT}" != "null" ]]; then
  STORAGE_SSH_CHECK="-p ${STORAGE_SERVER_SSH_PORT}"
fi

# Where this host can reach the storage server itself, validate the folder over a
# single ssh hop. That is one less shell to quote for, and it answers faster.
# Epsilon cannot resolve the NAS hostnames at all, so it falls back to proxying
# the check through the data-transfer node.
if ssh -n -o BatchMode=yes -o ConnectTimeout=8 ${STORAGE_SSH_CHECK} "${REMOTE_STORAGE_HOST}" true >/dev/null 2>&1; then
  STORAGE_DIRECT=1
else
  STORAGE_DIRECT=0
fi

if [[ "${STORAGE_DIRECT}" == "1" ]]; then
  echo "Storage server reachable directly from this host."
  if ! ssh -n ${STORAGE_SSH_CHECK} "${REMOTE_STORAGE_HOST}" "test -d $(shq "${REMOTE_INPUT_PATH}")"; then
    echo "Folder not found on ${STORAGE_SERVER_NAME}: ${INPUT_FOLDER}" >&2
    exit 1
  fi
else
  # Two shell hops (DTN, then the storage server) means quoting twice: the inner
  # command is quoted for the storage server, then the whole thing for the DTN.
  REMOTE_TEST_COMMAND="test -d $(shq "${REMOTE_INPUT_PATH}")"
  REMOTE_CHECK_COMMAND="ssh ${STORAGE_SSH_CHECK} $(shq "${REMOTE_STORAGE_HOST}") $(shq "${REMOTE_TEST_COMMAND}")"
  if ! ssh "${DTN_HOST}" "${REMOTE_CHECK_COMMAND}"; then
    echo "Folder not found on ${STORAGE_SERVER_NAME}: ${INPUT_FOLDER}" >&2
    exit 1
  fi
fi

# Refuse to start a transfer that cannot finish. Three runs in August 2026 died
# partway through the pull with a bare "Disk quota exceeded (122)" from rsync,
# which tells the submitting user nothing they can act on.
REMOTE_SIZE_COMMAND="du -sk $(shq "${REMOTE_INPUT_PATH}") 2>/dev/null | cut -f1"
if [[ "${STORAGE_DIRECT}" == "1" ]]; then
  SOURCE_SIZE_KB="$(ssh -n ${STORAGE_SSH_CHECK} "${REMOTE_STORAGE_HOST}" "${REMOTE_SIZE_COMMAND}" 2>/dev/null || echo "")"
else
  SOURCE_SIZE_KB="$(ssh "${DTN_HOST}" "ssh ${STORAGE_SSH_CHECK} $(shq "${REMOTE_STORAGE_HOST}") $(shq "${REMOTE_SIZE_COMMAND}")" 2>/dev/null || echo "")"
fi
AVAIL_KB="$(ssh "${DTN_HOST}" "df -Pk $(shq "${PIPELINE_ROOT}") | awk 'NR==2 {print \$4}'" 2>/dev/null || echo "")"

if [[ "${SOURCE_SIZE_KB}" =~ ^[0-9]+$ && "${AVAIL_KB}" =~ ^[0-9]+$ ]]; then
  # Ask for headroom: the workflow writes trimmed FASTQ, an assembled FASTA and
  # per-well splits on top of the raw input.
  NEEDED_KB=$(( SOURCE_SIZE_KB * 3 ))
  if (( AVAIL_KB < NEEDED_KB )); then
    echo "Not enough space on the cluster to run ${RUN_NAME}." >&2
    echo "The run needs about $(( NEEDED_KB / 1024 )) MB of working space but only $(( AVAIL_KB / 1024 )) MB is free at ${PIPELINE_ROOT}." >&2
    echo "Free space under the pipeline results folder, or ask the HPC team to raise the quota, then submit again." >&2
    exit 1
  fi
  echo "Space check: $(( SOURCE_SIZE_KB / 1024 )) MB input, $(( AVAIL_KB / 1024 )) MB free."
else
  echo "Space check skipped: could not read source size or free space." >&2
fi

# Pull the raw run folder onto the DTN first, excluding Synology metadata
# folders that otherwise confuse FASTQ discovery downstream.
REMOTE_PULL_COMMAND="mkdir -p $(shq "${PIPELINE_ROOT}/data/runs/${USER_ID}") && rsync -a -P --exclude='@eaDir/' --exclude='@eaDir'"
if [[ -n "${RSYNC_SSH_OPTION}" ]]; then
  REMOTE_PULL_COMMAND="${REMOTE_PULL_COMMAND} ${RSYNC_SSH_OPTION}"
fi
REMOTE_PULL_COMMAND="${REMOTE_PULL_COMMAND} $(shq "${RAWDATA_SERVER}${INPUT_FOLDER}") $(shq "${PIPELINE_ROOT}/data/runs/${USER_ID}/")"

ssh "${DTN_HOST}" "${REMOTE_PULL_COMMAND}"
rsync -a -P "${RUNS_DIR}/${USER_ID}/${RUN_NAME}/" "${DTN_HOST}:${PIPELINE_ROOT}/data/runs/${USER_ID}/${RUN_NAME}/"
# Emit a lightweight marker that the Shiny app can parse to improve the
# runtime estimator without opening or recomputing FASTQ sizes itself.
# Keep the remote size probe POSIX-friendly. The earlier grouped `find (...)`
# expression was valid bash syntax but could still trip plain `sh` on some
# hosts once it passed through SSH quoting layers.
INPUT_BYTES="$(ssh "${DTN_HOST}" "find '${PIPELINE_ROOT}/data/runs/${USER_ID}/${RUN_NAME}' -maxdepth 2 -type f ! -path '*/@eaDir/*' -printf '%f\t%s\n' | awk -F '\t' '\$1 ~ /\\.fastq(\\.gz)?$/ {sum += \$2} END {print sum + 0}'")"
echo "TBRC_INPUT_BYTES=${INPUT_BYTES:-0}"
dispatch_detached "${SSH_HOST}" "${PIPELINE_ROOT}" "${PATH_PREFIX}" "${SNAKEMAKE_ENV}" \
  "${RUN_KEY}" "${DATA_TABLE}" "${INPUTPATH}"
