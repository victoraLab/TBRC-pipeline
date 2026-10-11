#!/bin/bash

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

require_tool rsync

DTN_HOST="$(json_get_cluster_field dtn_host)"
PIPELINE_ROOT="$(json_get_cluster_field pipeline_root)"
SOURCE_DIRECTORY="${FRONTEND_DIR}/bc/"
TARGET_DIRECTORY="${DTN_HOST}:${PIPELINE_ROOT}/data/barcodes/"

rsync -av --exclude="bc_names.txt" "${SOURCE_DIRECTORY}" "${TARGET_DIRECTORY}"
