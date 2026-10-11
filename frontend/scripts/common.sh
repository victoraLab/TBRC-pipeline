#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRONTEND_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${FRONTEND_DIR}/.." && pwd)"
BACKEND_DIR="${REPO_ROOT}/backend"
CONFIG_FILE="${FRONTEND_DIR}/server_config.json"

# Where run metadata is staged before it is shipped to the cluster. The Shiny app
# resolves the same location from the same variable, so the two cannot disagree.
RUNS_DIR="${TBRC_RUNS_DIR:-runs}"

require_tool() {
    local tool_name="$1"

    if ! command -v "${tool_name}" >/dev/null 2>&1; then
        echo "Required command not found: ${tool_name}" >&2
        exit 1
    fi
}

json_get_storage_field() {
    local server_name="$1"
    local field_name="$2"

    jq -r --arg server_name "${server_name}" --arg field_name "${field_name}" '
        .storage_server[$server_name] as $entry
        | if $entry == null then
              empty
          elif ($entry | type) == "string" then
              if $field_name == "rawdata_path" then
                  $entry
              elif $field_name == "results_path" then
                  ($entry | sub("rawdata/?$"; "results/"))
              else
                  empty
              end
          else
              ($entry[$field_name] // empty)
          end
    ' "${CONFIG_FILE}"
}

json_get_cluster_field() {
    local field_name="$1"
    local value

    value="$(jq -r --arg field_name "${field_name}" '.cluster[$field_name] // empty' "${CONFIG_FILE}")"
    # The app runs as the account whose ssh config maps these hosts to the
    # cluster user. A terminal user on another machine may not have that, so
    # ./tbrc exports TBRC_SSH_USER and a bare host gets it here. Unset, nothing
    # changes.
    case "${field_name}" in
        ssh_host|dtn_host)
            if [[ -n "${TBRC_SSH_USER:-}" && -n "${value}" && "${value}" != *@* ]]; then
                value="${TBRC_SSH_USER}@${value}"
            fi
            ;;
    esac
    printf '%s' "${value}"
}

rsync_transport_args() {
    local ssh_port="$1"

    if [[ -n "${ssh_port}" && "${ssh_port}" != "null" ]]; then
        printf -- "-e 'ssh -p %s'" "${ssh_port}"
    fi
}

ensure_remote_directory() {
    local remote_spec="$1"
    local remote_parent="${remote_spec%/}"
    local remote_host="${remote_parent%%:*}"
    local remote_path="${remote_parent#*:}"

    if [[ "${remote_host}" == "${remote_path}" ]]; then
        echo "Remote destination is missing host:path format: ${remote_spec}" >&2
        exit 1
    fi

    ssh "${remote_host}" "mkdir -p '${remote_path}'"
}

# Quote a value for safe interpolation into a shell command string that will be
# evaluated on a remote host. Every user-supplied value that ends up inside an
# ssh command string must pass through this, once per shell hop.
shq() {
    printf '%q' "${1:-}"
}

# Run names, user ids and folder names become path segments on three different
# machines and are interpolated into remote shell commands. Restrict them to a
# conservative character set and reject anything that could traverse or escape.
validate_path_segment() {
    local value="${1:-}"
    local label="${2:-value}"

    if [[ -z "${value}" ]]; then
        echo "${label} cannot be empty." >&2
        exit 1
    fi

    if [[ ! "${value}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "${label} may only contain letters, numbers, dot, underscore and hyphen: ${value}" >&2
        exit 1
    fi

    if [[ "${value}" == "." || "${value}" == ".." || "${value}" == -* ]]; then
        echo "${label} is not a valid folder name: ${value}" >&2
        exit 1
    fi
}
