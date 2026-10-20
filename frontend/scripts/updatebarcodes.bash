#!/bin/bash

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

require_tool rsync

DTN_HOST="$(json_get_cluster_field dtn_host)"
PIPELINE_ROOT="$(json_get_cluster_field pipeline_root)"
# The app writes uploaded sets to TBRC_BC_DIR when it is set (a volume in the
# the Docker host container) and to frontend/bc otherwise. This used to read
# frontend/bc unconditionally, so in the container a new set was saved and
# never reached the cluster, and its runs failed for want of barcodes.
SOURCE_DIRECTORY="${TBRC_BC_DIR:-${FRONTEND_DIR}/bc}/"
TARGET_DIRECTORY="${DTN_HOST}:${PIPELINE_ROOT}/data/barcodes/"

# Add new sets only; never rewrite one the cluster already has. The app's tree
# and the cluster's are not identical (CBC is stored reverse-complemented with
# BC1 and BC2 swapped in frontend/bc; B2's BC1 lacked its final newline), and a
# plain rsync after any upload would have replaced the cluster's copies, the
# ones runs actually read. The app refuses an upload under an existing folder
# name, so a new set is always a new folder.
rsync -av --ignore-existing --exclude="bc_names.txt" "${SOURCE_DIRECTORY}" "${TARGET_DIRECTORY}"
