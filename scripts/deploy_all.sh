#!/bin/bash
set -euo pipefail

# Deploys the backend to the cluster. The frontend is not deployed from here: the app runs
# from the checkout (or the image) it was started from, so update that and restart the app.
# With the Docker setup in deploy/docker that is a git pull and a restart of the container.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEPLOY_CONFIG="${SCRIPT_DIR}/deploy_config.sh"

# The cluster refuses a sample.json whose servers or program paths are not in
# backend/config/trusted_servers.json, which is generated from
# frontend/server_config.json. Stop before copying anything if the two have
# drifted apart (fix: python3 scripts/make_trusted_servers.py, then commit).
#
# The list is optional. An installation that has not generated one (every fresh
# clone of the public repository) has no such file, and the cluster then checks
# formats only; it is not an error, so say so once and go on. python3
# scripts/make_trusted_servers.py turns the list on, and from then on this check
# runs on every deploy.
TRUSTED_LIST="${ROOT}/backend/config/trusted_servers.json"
if [[ -f "${TRUSTED_LIST}" ]]; then
  python3 "${ROOT}/scripts/make_trusted_servers.py" --check
else
  echo "Trusted-server list is off: no backend/config/trusted_servers.json, so the cluster checks formats only (python3 scripts/make_trusted_servers.py turns it on)."
fi

if [[ ! -f "${DEPLOY_CONFIG}" ]]; then
  echo "Missing local deploy config: ${DEPLOY_CONFIG}" >&2
  echo "Copy scripts/deploy_config.example.sh to scripts/deploy_config.sh and fill in your private targets." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "${DEPLOY_CONFIG}"

: "${HPC_TARGET:?HPC_TARGET is required in scripts/deploy_config.sh}"

# --dry-run lists what would change and copies nothing.
DRY=()
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY=(--dry-run --itemize-changes)
fi

# Code only. Each host owns its state: this used to copy everything under
# backend/, so a deploy from a working copy carried that copy's runs, logs and
# results to the cluster. Files are only ever added or updated on the far side,
# never removed, so a file dropped here has to be dropped there by hand.
# The cluster's IgBLAST databases, internal_data and binaries are built there;
# config/site.env is host-local by design (see config/site.env.example).
BACKEND_EXCLUDES=(
  --exclude '.snakemake/' --exclude '.parallel/' --exclude 'data/runs/'
  --exclude 'results/' --exclude 'logs/' --exclude 'slurm_outs/'
  --exclude 'config/site.env'
  --exclude 'igblast/work/' --exclude 'igblast/internal_data/'
  --exclude 'igblast/bin/' --exclude 'igblast/refs/'
  --exclude '__pycache__/' --exclude 'Rplots.pdf'
)

# ${DRY[@]+...}: macOS's bash 3.2 calls an empty array unbound under set -u, so
# a plain "${DRY[@]}" aborted every real deploy and only --dry-run ever ran.
rsync -aP ${DRY[@]+"${DRY[@]}"} "${BACKEND_EXCLUDES[@]}" \
  "${ROOT}/backend/" \
  "${HPC_TARGET}"
