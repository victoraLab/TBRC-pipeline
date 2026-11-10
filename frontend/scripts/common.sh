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

    ssh "${remote_host}" "mkdir -p $(shq "${remote_path}")"
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

# ---- the cluster gate -----------------------------------------------------------
# With cluster.use_gate true in server_config.json, the cluster's login and transfer nodes
# accept from this host only the named actions of deploy/cluster/README.md, and an rsync
# push in three fixed shapes. The launchers then ask for an action by name and never
# compose a command for a remote shell. Without it (the Docker setup in deploy/docker, the public export)
# they compose the commands as before.

# On only when cluster.use_gate is exactly true (a missing key, false, null, a string, a
# number or a config that does not parse all mean off), and never when TBRC_USE_GATE is
# exactly 0. The tbrc command line sets that for its server-mode call: it reaches the
# cluster with the user's own key, which the gate does not cover, so it must ask for the
# commands themselves. Any other value of TBRC_USE_GATE, or none, leaves it to the config.
gate_enabled() {
    if [[ "${TBRC_USE_GATE:-}" == "0" ]]; then
        return 1
    fi
    jq -e '.cluster.use_gate == true' "${CONFIG_FILE}" >/dev/null 2>&1
}

# A path on the cluster: relative to the gate's own root when the gate is on (the gate maps
# it under the pipeline root it was installed with), under the configured pipeline_root
# otherwise. Usage: cluster_path data/runs/
cluster_path() {
    if gate_enabled; then
        printf '%s' "$1"
    else
        printf '%s/%s' "$(json_get_cluster_field pipeline_root)" "$1"
    fi
}

# Ask the cluster for one named action: ssh HOST ACTION ARG... Each ARG is sent as a bare
# word. The gate splits what it receives on single spaces and refuses anything that is not
# printable ASCII, so quoting added here would only make it refuse. A value that is not a
# plain word (empty, a space, a quote, a newline) is refused here instead, as the gate would
# refuse it and with its status, 3, so that nothing odd is sent and the caller sees the same
# failure as for any remote call that fails. GATE_CONNECT_TIMEOUT sets ssh's ConnectTimeout.
# Usage: gate_call SSH_HOST ACTION [ARG...]
gate_call() {
    local host="$1" arg
    shift
    for arg in "$@"; do
        if [[ ! "${arg}" =~ ^[A-Za-z0-9._-]+$ ]]; then
            echo "tbrc-gate: not sent, a value is not a plain word" >&2
            return 3
        fi
    done
    ssh -n -o BatchMode=yes -o ConnectTimeout="${GATE_CONNECT_TIMEOUT:-20}" "${host}" "$@"
}

# ---- following a run on the login node ------------------------------------------------
# A run's Snakemake driver is started on the login node detached from the ssh session,
# then followed until it exits, relaying its output as it grows and exiting with its
# status. The driver used to run in the foreground of the ssh command, so a dropped
# connection or a restart of the app host killed a run that was otherwise healthy.
# run_isolated.sh writes the exit code to .parallel/<key>/attempts/<attempt>.exit after
# the driver finishes; that file, not the ssh session, says when the run is over. A failed
# poll (a network blip, a login node rebooting) is not a verdict on the run and is
# retried. A driver that is gone without leaving that file (killed with the login node,
# say) would be waited for forever, so three polls in a row that find neither a live
# driver nor an exit code end the wait as a failure.
#
# Both ways of starting a run (dispatch_detached, dispatch_via_gate) follow it with the
# loop below; they differ in how they start it and in how they ask about it.

# How long a run's wrapper puts up with the gate refusing its status check before it
# gives up on the wait: about ten minutes of polls in a row that end with status 3 (the
# gate refused the request) or 4 (the gate is not installed properly). Asking again does
# not help either, unlike a dropped connection, so a permanent refusal must not hold the
# wrapper, and with it the run's name, for ever. The run itself is not touched.
GATE_POLL_REFUSED_SECONDS=600

# How many polls in a row the gate may refuse: GATE_POLL_REFUSED_SECONDS divided by the
# pause between polls (TBRC_POLL_SECONDS, 60 unless set), rounded up, at least 1.
# TBRC_GATE_POLL_REFUSED_LIMIT sets it outright (the tests do).
gate_poll_refused_limit() {
    local interval="${TBRC_POLL_SECONDS:-60}"
    if [[ "${TBRC_GATE_POLL_REFUSED_LIMIT:-}" =~ ^[1-9][0-9]*$ ]]; then
        echo "${TBRC_GATE_POLL_REFUSED_LIMIT}"
        return 0
    fi
    [[ "${interval}" =~ ^[1-9][0-9]*$ ]] || interval=1
    echo $(( (GATE_POLL_REFUSED_SECONDS + interval - 1) / interval ))
}

# Follow an attempt until it ends. POLL_COMMAND [ARG...] SEEN is run for each look; it
# prints the log lines after the first SEEN, then __TBRC_EXIT__=<n> once the run is over
# and __TBRC_ALIVE__ while its driver is running, and fails if it could not ask. Any
# failed look is retried, except that with REFUSALS_END_IT=1 (the gate's poll) enough
# looks in a row that fail with status 3 or 4 end the wait with 70; an answer, or a
# failure of another kind (255: no connection), starts that count again.
# Usage: follow_attempt SSH_HOST RUN_KEY ATTEMPT LOG_HINT REFUSALS_END_IT POLL_COMMAND [ARG...]
follow_attempt() {
    local ssh_host="$1" key="$2" attempt="$3" log_hint="$4" refusals_end_it="$5"
    shift 5
    local chunk body code="" gone=0 seen=0 poll="${TBRC_POLL_SECONDS:-60}"
    local poll_rc refused=0 refused_limit=0
    if [[ "${refusals_end_it}" == "1" ]]; then
        refused_limit="$(gate_poll_refused_limit)"
    fi
    while [[ -z "${code}" ]]; do
        sleep "${poll}"
        # The exit file is read before the log: it is written after the last
        # line of output, so a log read after seeing it is complete.
        if chunk="$("$@" "${seen}")"; then
            refused=0
            body="$(printf '%s\n' "${chunk}" | grep -v '^__TBRC_EXIT__=\|^__TBRC_ALIVE__$' || true)"
            if [[ -n "${body}" ]]; then
                printf '%s\n' "${body}"
                seen=$((seen + $(printf '%s\n' "${body}" | wc -l)))
            fi
            code="$(printf '%s\n' "${chunk}" | sed -n 's/^__TBRC_EXIT__=//p' | tail -1)"
            if [[ -z "${code}" ]] && ! grep -q '^__TBRC_ALIVE__$' <<< "${chunk}"; then
                gone=$((gone + 1))
                if (( gone >= 3 )); then
                    echo "The driver for ${key} is no longer running and left no exit status (attempt ${attempt}); treating the run as failed. Its log: ${log_hint}" >&2
                    return 70
                fi
            else
                gone=0
            fi
        else
            poll_rc=$?
            if (( refused_limit > 0 )) && [[ "${poll_rc}" == 3 || "${poll_rc}" == 4 ]]; then
                refused=$((refused + 1))
                if (( refused >= refused_limit )); then
                    echo "The cluster gate refused the status check for ${key} ${refused} times in a row (attempt ${attempt}); giving up on following it. The run may still be running on the cluster. Its log: ${log_hint}" >&2
                    return 70
                fi
            else
                refused=0
            fi
            echo "Could not reach ${ssh_host} to check on ${key}; trying again." >&2
        fi
    done
    return "${code}"
}

# One look at a run, by a shell command composed here. SEEN is last, as follow_attempt passes it.
# Usage: poll_shell SSH_HOST PIPELINE_ROOT EXIT_FILE LOG PIDFILE SEEN
poll_shell() {
    local ssh_host="$1" root="$2" exit_file="$3" log="$4" pidfile="$5" seen="$6"
    ssh -o BatchMode=yes -o ConnectTimeout=20 "${ssh_host}" "cd $(shq "${root}") \
                && e=\$(cat $(shq "${exit_file}") 2>/dev/null); tail -n +$((seen + 1)) $(shq "${log}") 2>/dev/null; \
                [ -n \"\$e\" ] && echo __TBRC_EXIT__=\$e; \
                p=$(shq "${pidfile}"); [ -f \"\$p\" ] && kill -0 \$(cat \"\$p\") 2>/dev/null && echo __TBRC_ALIVE__; true"
}

# One look at a run, by the gate's poll action. SEEN is last, as follow_attempt passes it.
# Usage: poll_via_gate SSH_HOST USER RUN ATTEMPT SEEN
poll_via_gate() {
    gate_call "$1" poll "$2" "$3" "$4" "$5"
}

# Start a run's Snakemake driver on the login node detached from this ssh session, then
# follow it (follow_attempt).
#
# Usage: dispatch_detached SSH_HOST PIPELINE_ROOT PATH_PREFIX SNAKEMAKE_ENV RUN_KEY DATA_TABLE INPUTPATH
dispatch_detached() {
    local ssh_host="$1" root="$2" path_prefix="$3" env_name="$4" key="$5" table="$6" inputs="$7"
    local attempt work log exit_file
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
        && echo $(shq "Run ${key} started on ${ssh_host}, attempt ${attempt}.")"

    follow_attempt "${ssh_host}" "${key}" "${attempt}" "${root}/${log}" 0 \
        poll_shell "${ssh_host}" "${root}" "${exit_file}" "${log}" "${work}/.snakemake/tbrc-driver.pid"
}

# The same through the gate: its launch action starts run_isolated.sh for data/runs/<user>/<run>/
# (table and inputs are fixed there, as the launchers pass them) and its poll action answers.
# Usage: dispatch_via_gate SSH_HOST USER RUN
dispatch_via_gate() {
    local ssh_host="$1" user="$2" run="$3" key attempt
    key="${user}-${run}"
    attempt="$(date +%Y%m%d-%H%M%S)-$$"
    gate_call "${ssh_host}" launch "${user}" "${run}" "${attempt}" || return $?
    follow_attempt "${ssh_host}" "${key}" "${attempt}" ".parallel/${key}/attempts/${attempt}.log (under the pipeline root)" 1 \
        poll_via_gate "${ssh_host}" "${user}" "${run}" "${attempt}"
}

# ---- runs on this host ("Run on this server") ---------------------------------
# A cluster run gets its folder (sample.json, a MIBS design, uploaded reads)
# shipped to <pipeline_root>/data/runs/<user>/<run>/ and its raw data pulled
# there by the DTN. A run on this host needs the same two things in this host's
# backend/data/runs/<user>/<run>/. The local branch of goServer.bash used to
# start Snakemake on that folder without putting anything in it, so every local
# run stopped at DAG build ("Expected exactly one R1 and one R2 FASTQ.gz").

# Refuse before any data moves when this host cannot run the workflow. The app
# hides the switch in that case; this guards a direct call.
require_local_backend() {
    local tool
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
}

# Usage: stage_run_locally USER RUN -- the app's run folder into the workflow's.
stage_run_locally() {
    local target="${BACKEND_DIR}/data/runs/$1/$2"
    mkdir -p "${target}"
    rsync -a "${RUNS_DIR}/$1/$2/" "${target}/"
}

# Bytes of FASTQ in a folder, for the app's runtime estimate (TBRC_INPUT_BYTES).
local_fastq_bytes() {
    local f total=0
    while IFS= read -r -d '' f; do
        total=$(( total + $(wc -c < "${f}") ))
    done < <(find "$1" -maxdepth 2 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) \
                 ! -path '*/@eaDir/*' -print0 2>/dev/null)
    echo "${total}"
}

# Run the workflow here, in the foreground; the app already runs this script in
# the background. Usage: run_locally RUN_KEY DATA_TABLE INPUTPATH PATH_PREFIX SNAKEMAKE_ENV
run_locally() {
    local key="$1" table="$2" inputs="$3" path_prefix="$4" env_name="$5"
    cd "${BACKEND_DIR}"
    # PATH_PREFIX and the named Conda environment describe the cluster. In a
    # container the tools are already on PATH in Conda's base environment, and
    # forcing the cluster's PATH would hide them. Apply each only where it exists.
    if [[ -n "${path_prefix}" && -d "${path_prefix%%:*}" ]]; then
        export PATH="${path_prefix}"
    fi
    if [[ -n "${env_name}" && "${env_name}" != "base" ]] \
       && command -v conda >/dev/null 2>&1 \
       && conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${env_name}"; then
        conda activate "${env_name}"
    fi
    # Jobs run on this host, not through sbatch (backend/snakemake.sh).
    export DATA_TABLE="${table}" INPUTPATH="${inputs}" TBRC_EXECUTOR=local
    ./run_isolated.sh "${key}" "${table}" "${inputs}"
}
