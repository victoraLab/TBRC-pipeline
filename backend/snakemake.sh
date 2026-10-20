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
# CMP13's demultiplexing takes about ten and a half minutes, which made it a coin
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
snakemake \
    -j 9999 \
    --rerun-incomplete \
    --cluster-config config/cluster.json \
    --cluster "sbatch --partition {cluster.partition} --job-name {cluster.name} -o {cluster.output} -e {cluster.error} --ntasks-per-node 1 --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus_per_task}" \
    "$@"
