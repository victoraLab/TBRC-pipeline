#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEPLOY_CONFIG="${SCRIPT_DIR}/deploy_config.sh"

if [[ ! -f "${DEPLOY_CONFIG}" ]]; then
  echo "Missing local deploy config: ${DEPLOY_CONFIG}" >&2
  echo "Copy scripts/deploy_config.example.sh to scripts/deploy_config.sh and fill in your private targets." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "${DEPLOY_CONFIG}"

: "${SHINY_TARGET:?SHINY_TARGET is required in scripts/deploy_config.sh}"
: "${HPC_TARGET:?HPC_TARGET is required in scripts/deploy_config.sh}"
: "${SHINY_SSH_PORT:=6123}"

rsync -aP -e "ssh -p ${SHINY_SSH_PORT}" \
  "${ROOT}/frontend/" \
  "${SHINY_TARGET}"

rsync -aP \
  "${ROOT}/backend/" \
  "${HPC_TARGET}"
