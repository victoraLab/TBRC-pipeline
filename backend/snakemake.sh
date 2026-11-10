#!/bin/bash
set -uo pipefail

# Launchers call this through run_isolated.sh, which gives every run its own
# working directory. Until October 2026 every submission ran here in the
# pipeline root instead and shared one .snakemake state directory. Two failure
# modes followed from that, and between them they caused 12 of the 30 recorded
# run failures:
#
#   * an interrupted run leaves output files marked incomplete, and the NEXT
#     user's run refuses to start until someone passes --rerun-incomplete;
#   * an interrupted run leaves a lock behind, and every later run refuses to
#     start with "Directory cannot be locked" until someone unlocks it.
#
# Neither is the submitting user's fault and neither is fixable from the web UI,
# so clear both here before starting. Both still apply inside a per-run
# directory: a resubmitted run reuses its directory and its leftovers.

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "${WORKDIR}"

# Is another run already going? This guard used to ask
# `pgrep -f "snakemake.*${WORKDIR}"`, which never matched: snakemake is invoked
# without --directory, so the working directory does not appear on its command
# line. The guard therefore always concluded "no run active", declared the live
# run's lock stale, cleared it, and started a second Snakemake on the same
# directory -- two drivers dispatching the same rules to the same output paths.
#
# A PID file says exactly who holds the directory, and cannot be fooled by a
# pattern that happens to match this script's own command line.
PIDFILE="${WORKDIR}/.snakemake/tbrc-driver.pid"
mkdir -p "${WORKDIR}/.snakemake"

if [[ -f "${PIDFILE}" ]]; then
    OTHER_PID="$(cat "${PIDFILE}" 2>/dev/null || true)"
    if [[ -n "${OTHER_PID}" ]] && kill -0 "${OTHER_PID}" 2>/dev/null; then
        echo "Another Snakemake run (pid ${OTHER_PID}) is active in ${WORKDIR}. Refusing to start a second one." >&2
        echo "Wait for the running submission to finish, then try again." >&2
        exit 1
    fi
    echo "Removing a PID file left by an interrupted run (pid ${OTHER_PID:-unknown} is gone)."
    rm -f "${PIDFILE}"
fi

echo "$$" > "${PIDFILE}"
# Release the directory however this script exits, including on kill.
trap 'rm -f "${PIDFILE}"' EXIT INT TERM

if [[ -d .snakemake/locks ]] && [[ -n "$(ls -A .snakemake/locks 2>/dev/null)" ]]; then
    echo "Clearing a Snakemake lock left by an interrupted run."
    snakemake --unlock --cores 1 >/dev/null 2>&1 || rm -f .snakemake/locks/*
fi

# cluster.json has always declared time, mem and cpus_per_task, but the --cluster
# line below never passed them to sbatch, so every rule inherited the partition
# default instead: DefaultTime on the cluster partition is TEN MINUTES. Any step that ran
# longer was killed with "DUE TO TIME LIMIT" no matter what cluster.json said.
# RUN13's demultiplexing takes about ten and a half minutes, which made it a coin
# flip -- it passed once and timed out on the next identical run. Larger runs
# would have failed outright and looked like pipeline errors.
#
# Host-local settings that must not live in git (config/site.env is not tracked).
# Today that is CLONALITY_LIB, an R library holding a newer clonality than the
# conda env's own; sbatch passes the submitting environment on to every rule.
if [[ -f "${WORKDIR}/config/site.env" ]]; then
    # shellcheck source=/dev/null
    source "${WORKDIR}/config/site.env"
    export CLONALITY_LIB
fi

# --rerun-incomplete: redo anything a previous interrupted run left half-written.
# Extra arguments go to snakemake as they are; rerun_final.sh uses them to ask
# for the refinal target instead of the whole run.
#
# A run on the app's own host ("Run on this server": run_locally in
# frontend/scripts/common.sh) has no SLURM. It sets TBRC_EXECUTOR=local and the
# jobs run here, on TBRC_LOCAL_CORES cores (default: every core this host or
# container has). Until this existed, every such run died at its first job with
# "sbatch: not found". Only the explicit setting selects it, never a missing
# sbatch, so a cluster host with a broken PATH fails loudly instead of running
# every rule on the login node.
if [[ "${TBRC_EXECUTOR:-cluster}" == "local" ]]; then
    snakemake \
        --cores "${TBRC_LOCAL_CORES:-$(nproc 2>/dev/null || echo 4)}" \
        --rerun-incomplete \
        "$@"
else
    # The SLURM account to charge is __default__.account in config/cluster.json. Empty or absent
    # means none: a cluster whose partitions need no account gets no --account, and sbatch uses the
    # user's default. This lab's file names its account, so its sbatch line is what it always was.
    # Until November 2026 the account was written into the line below.
    #
    # Snakemake runs the --cluster string through a shell, so the value must not reach it as it
    # stands. A real account name is letters, digits, dot, underscore and hyphen, so anything else
    # stops the run here, naming the field and not the value; the name is quoted as well, which does
    # nothing to a valid one. A file that cannot be read stops the run too: carrying on would drop
    # the account silently and the first sbatch would be refused far from the cause.
    if ! CLUSTER_ACCOUNT="$(python3 -c '
import json, sys
try:
    with open(sys.argv[1], encoding="utf-8") as handle:
        value = json.load(handle)["__default__"].get("account")
except (OSError, ValueError, KeyError, TypeError, AttributeError) as err:
    sys.exit("cannot read __default__ from %s: %s" % (sys.argv[1], err))
if value is not None and not isinstance(value, str):
    sys.exit("__default__.account in %s must be a string" % sys.argv[1])
sys.stdout.write(value or "")
' config/cluster.json)"; then
        echo "snakemake.sh: config/cluster.json must be a JSON file with a __default__ section (see above)." >&2
        exit 1
    fi
    SBATCH_ACCOUNT=""
    if [[ -n "${CLUSTER_ACCOUNT}" ]]; then
        if [[ ! "${CLUSTER_ACCOUNT}" =~ ^[A-Za-z0-9_][A-Za-z0-9._-]*$ ]]; then
            echo "snakemake.sh: __default__.account in config/cluster.json must be a SLURM account name: letters, digits, dot, underscore and hyphen only, not starting with a hyphen." >&2
            exit 1
        fi
        SBATCH_ACCOUNT="--account $(printf '%q' "${CLUSTER_ACCOUNT}") "
    fi
    snakemake \
        -j 9999 \
        --rerun-incomplete \
        --cluster-config config/cluster.json \
        --cluster "sbatch ${SBATCH_ACCOUNT}--partition {cluster.partition} --job-name {cluster.name} -o {cluster.output} -e {cluster.error} --ntasks-per-node 1 --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus_per_task}" \
        "$@"
fi
