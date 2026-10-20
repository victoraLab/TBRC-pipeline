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

# Start a run's Snakemake driver on the login node detached from this ssh
# session, then follow it until it exits, relaying its output as it grows and
# exiting with its status. The driver used to run in the foreground of the ssh
# command, so a dropped connection or a restart of the app host killed a run
# that was otherwise healthy. run_isolated.sh writes the exit code to
# .parallel/<key>/attempts/<attempt>.exit after the driver finishes; that file,
# not the ssh session, says when the run is over. A failed poll (a network
# blip, a login node rebooting) is not a verdict on the run and is retried. A
# driver that is gone without leaving that file (killed with the login node,
# say) would be waited for forever, so three polls in a row that find neither a
# live driver nor an exit code end the wait as a failure.
#
# Usage: dispatch_detached SSH_HOST PIPELINE_ROOT PATH_PREFIX SNAKEMAKE_ENV RUN_KEY DATA_TABLE INPUTPATH
dispatch_detached() {
    local ssh_host="$1" root="$2" path_prefix="$3" env_name="$4" key="$5" table="$6" inputs="$7"
    local attempt work log exit_file seen=0 poll="${TBRC_POLL_SECONDS:-60}"
    attempt="$(date +%Y%m%d-%H%M%S)-$$"
    work=".parallel/${key}"
    log="${work}/attempts/${attempt}.log"
    exit_file="${work}/attempts/${attempt}.exit"

    ssh "${ssh_host}" "cd $(shq "${root}") \
        && export PATH=$(shq "${path_prefix}") \
        && conda activate $(shq "${env_name}") \
        && mkdir -p $(shq "${work}/attempts") \
        && { TBRC_ATTEMPT=$(shq "${attempt}") setsid nohup ./run_isolated.sh $(shq "${key}") \
               $(shq "${table}") $(shq "${inputs}") > $(shq "${log}") 2>&1 < /dev/null & } \
        && echo \"Run ${key} started on ${ssh_host}, attempt ${attempt}.\""

    local chunk body code="" gone=0 pidfile="${work}/.snakemake/tbrc-driver.pid"
    while [[ -z "${code}" ]]; do
        sleep "${poll}"
        # The exit file is read before the log: it is written after the last
        # line of output, so a log read after seeing it is complete.
        if chunk="$(ssh -o BatchMode=yes -o ConnectTimeout=20 "${ssh_host}" "cd $(shq "${root}") \
                && e=\$(cat $(shq "${exit_file}") 2>/dev/null); tail -n +$((seen + 1)) $(shq "${log}") 2>/dev/null; \
                [ -n \"\$e\" ] && echo __TBRC_EXIT__=\$e; \
                p=$(shq "${pidfile}"); [ -f \"\$p\" ] && kill -0 \$(cat \"\$p\") 2>/dev/null && echo __TBRC_ALIVE__; true")"; then
            body="$(printf '%s\n' "${chunk}" | grep -v '^__TBRC_EXIT__=\|^__TBRC_ALIVE__$' || true)"
            if [[ -n "${body}" ]]; then
                printf '%s\n' "${body}"
                seen=$((seen + $(printf '%s\n' "${body}" | wc -l)))
            fi
            code="$(printf '%s\n' "${chunk}" | sed -n 's/^__TBRC_EXIT__=//p' | tail -1)"
            if [[ -z "${code}" ]] && ! grep -q '^__TBRC_ALIVE__$' <<< "${chunk}"; then
                gone=$((gone + 1))
                if (( gone >= 3 )); then
                    echo "The driver for ${key} is no longer running and left no exit status (attempt ${attempt}); treating the run as failed. Its log: ${root}/${log}" >&2
                    return 70
                fi
            else
                gone=0
            fi
        else
            echo "Could not reach ${ssh_host} to check on ${key}; trying again." >&2
        fi
    done
    return "${code}"
}
