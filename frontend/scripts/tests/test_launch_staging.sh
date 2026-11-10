#!/bin/bash
# Checks that the app's launchers put a run's sample.json and FASTQ pair where the
# workflow looks for them, in all four ways a run starts:
#
#   upload  + local     goUpload.bash -l    (this host runs the pipeline)
#   upload  + cluster   goUpload.bash       (shipped to the cluster)
#   storage + local     goServer.bash -l
#   storage + cluster   goServer.bash
#
# Why it exists: the local branch of goServer.bash started Snakemake on
# data/runs/<user>/<run>/ without ever putting anything there, and goUpload.bash
# had no local branch and pointed the cluster at a folder named after the R1
# file while the app wrote sample.json under the run name. Both failed only at
# run time, on a real submission.
#
# A fifth part covers the same launchers, and preflight.bash, with cluster.use_gate
# on: the cluster accepts only named actions (deploy/cluster/README.md), so every
# ssh to it must begin with one, nothing may name the storage server, and the
# rsync pushes must keep the exact shape the gate's receive rule accepts.
#
# Nothing leaves this machine. The launchers run from a throwaway copy of the
# repository with a stub backend/run_isolated.sh that reports what it was given,
# stub pipeline tools, and stub ssh/rsync: a "storage server" is a local folder,
# and every cluster command is logged instead of sent. The stub ssh also plays the
# gate: it takes the named actions, refuses what the real gate refuses, and answers
# from that local folder.
#
# Run where bash >= 4 and jq exist, e.g. in the TBRC image from the repo root:
#   docker run --rm -v "$PWD:/repo" -w /repo tbrc-unified:2026-11-06 bash frontend/scripts/tests/test_launch_staging.sh
set -uo pipefail

SRC="$(cd "$(dirname "$0")/.." && pwd)"            # the real frontend/scripts
REAL_RSYNC="$(command -v rsync)"
T="$(mktemp -d)"
trap 'rm -rf "${T}"' EXIT
FAILS=0
pass() { echo "  PASS  $*"; }
fail() { echo "  FAIL  $*"; FAILS=$((FAILS + 1)); }

mkdir -p "${T}/repo/frontend" "${T}/repo/backend" "${T}/bin" "${T}/nas/raw" "${T}/runs"
cp -R "${SRC}" "${T}/repo/frontend/scripts"
cat > "${T}/repo/frontend/server_config.json" <<'JSON'
{ "storage_server": { "nas1": { "rawdata_path": "tester@nas.invalid:/raw/", "results_path": "tester@nas.invalid:/results/",
                                "ssh_port": 2222, "transfer_mode": "pull" } },
  "cluster": { "dtn_host": "tester@dtn.invalid", "ssh_host": "login.invalid", "pipeline_root": "/cluster/TBRC",
               "snakemake_env": "snakemake", "path_prefix": "/nonexistent/bin:/usr/bin" } }
JSON

# The workflow entry point: reports, then checks what the Snakefile would need.
cat > "${T}/repo/backend/run_isolated.sh" <<'SH'
#!/bin/bash
echo "RUN_ISOLATED key=$1 table=$2 input=$3 cwd=$(basename "$PWD")"
echo "EXECUTOR=${TBRC_EXECUTOR:-unset}"
[[ -f "$2" ]] && echo "TABLE_OK" || echo "TABLE_MISSING $2"
r1=$(ls "$3" 2>/dev/null | grep -c '_R1_.*\.fastq\.gz$'); r2=$(ls "$3" 2>/dev/null | grep -c '_R2_.*\.fastq\.gz$')
[[ "$r1" == 1 && "$r2" == 1 ]] && echo "READS_OK" || echo "READS_BAD R1=$r1 R2=$r2 in $3"
SH
chmod +x "${T}/repo/backend/run_isolated.sh"

for t in snakemake pandaseq cutadapt fastx_trimmer; do printf '#!/bin/sh\nexit 0\n' > "${T}/bin/${t}"; done
# What the stub cluster says when a launcher polls a detached run, by run key. The
# default run is over at once (exit 0). "chatty" shows lines, then ends; "failing"
# ends with exit 3; "vanished" leaves neither an exit code nor a live driver.
# Some runs also answer with a failed poll (poll_status, counted per run key):
# "refused" and "refusedfour" are refused (status 3, 4) every time, "slowgate" is
# refused five times and then fine, "flaky" is refused twice, answers, is refused
# twice more and answers again, "netdown" is unreachable (status 255) six times.
cat > "${T}/poll_answer.sh" <<'SH'
poll_status() {  # RUN_KEY: sets POLL_N, returns the status this poll is answered with
  local f="${T_DIR}/polls_$1"
  POLL_N=0; [[ -f "$f" ]] && POLL_N=$(cat "$f")
  POLL_N=$((POLL_N + 1)); echo "${POLL_N}" > "$f"
  case "$1" in
    *-refused) return 3 ;;
    *-refusedfour) return 4 ;;
    *-slowgate) (( POLL_N <= 5 )) && return 3 ;;
    *-flaky) (( POLL_N == 1 || POLL_N == 2 || POLL_N == 4 || POLL_N == 5 )) && return 3 ;;
    *-netdown) (( POLL_N <= 6 )) && return 255 ;;
  esac
  return 0
}
poll_answer() {  # RUN_KEY LINES_ALREADY_SEEN
  case "$1" in
    *-flaky) if (( POLL_N == 3 )); then echo "flaky line"; echo "__TBRC_ALIVE__"; else echo "__TBRC_EXIT__=0"; fi ;;
    *-chatty)
      case "$2" in
        0) echo "step one"; echo "step two"; echo "__TBRC_ALIVE__" ;;
        2) echo "step three"; echo "__TBRC_EXIT__=0" ;;
        *) echo "WRONG_SEEN=$2"; echo "__TBRC_EXIT__=0" ;;
      esac ;;
    *-failing) echo "boom"; echo "__TBRC_EXIT__=3" ;;
    *-vanished) ;;
    *) echo "__TBRC_EXIT__=0" ;;
  esac
}
SH
# The cluster's gate, as far as the launchers can tell: sourced by the stub ssh with the
# words ssh was given after the host. Like deploy/cluster/tbrc-gate.py it joins them with
# single spaces, takes a named action with exactly its arguments and refuses anything
# else with "tbrc-gate: refused" and exit 3. It answers from the fake NAS (nas/raw) and
# returns, doing nothing, when the words are not a named action (the older shell commands).
cat > "${T}/gate_stub.sh" <<'SH'
declare -A GATE_ARGS=([launch]=3 [poll]=4 [free-space]=0 [run-bytes]=2 [results-exists]=3 [pull-raw]=3
                      [nas-reachable]=1 [nas-dir]=2 [nas-size]=2 [nas-probe]=2 [nas-integrity]=4 [nas-trim-hits]=3)
gate_refuse() { echo "tbrc-gate: refused" >&2; exit 3; }
gate_segment() { [[ "$1" =~ ^[A-Za-z0-9._-]+$ && "$1" != . && "$1" != .. && "$1" != -* ]] || gate_refuse; }
gate_stub() {
  local joined="$*" act="${1:-}" re='^[a-z-]+( [A-Za-z0-9._-]+)*$' words nas="${T_DIR}/nas/raw"
  act="${joined%% *}"
  [[ -n "${act}" && -n "${GATE_ARGS[${act}]+x}" ]] || return 0
  [[ "${joined}" =~ ${re} ]] || gate_refuse
  read -ra words <<< "${joined}"
  (( ${#words[@]} - 1 == GATE_ARGS[${act}] )) || gate_refuse
  case "${act}" in
    launch)
      gate_segment "${words[1]}"; gate_segment "${words[2]}"
      [[ "${words[3]}" =~ ^[0-9]{8}-[0-9]{6}-[0-9]+$ ]] || gate_refuse
      echo "Run ${words[1]}-${words[2]} started on stub, attempt ${words[3]}." ;;
    poll)
      gate_segment "${words[1]}"; gate_segment "${words[2]}"
      [[ "${words[3]}" =~ ^[0-9]{8}-[0-9]{6}-[0-9]+$ && "${words[4]}" =~ ^[0-9]{1,9}$ ]] || gate_refuse
      . "${T_DIR}/poll_answer.sh"
      poll_status "${words[1]}-${words[2]}"; rc=$?
      if (( rc == 3 )); then gate_refuse; fi
      if (( rc != 0 )); then echo "tbrc-gate: stub status ${rc}" >&2; exit ${rc}; fi
      poll_answer "${words[1]}-${words[2]}" "${words[4]}" ;;
    free-space) cat "${T_DIR}/free_kb" ;;
    run-bytes) gate_segment "${words[1]}"; gate_segment "${words[2]}"; echo 4242 ;;
    results-exists) gate_segment "${words[1]}"; gate_segment "${words[2]}"; gate_segment "${words[3]}"
      [[ -d "${T_DIR}/cluster_results/${words[1]}/${words[2]}/${words[3]}" ]]; exit $? ;;
    pull-raw) gate_segment "${words[1]}"; gate_segment "${words[2]}"; gate_segment "${words[3]}"
      [[ -d "${nas}/${words[3]}" ]]; exit $? ;;
    *)  # the checks on the storage server: nas-reachable, nas-dir, nas-size, nas-probe, nas-integrity, nas-trim-hits
      gate_segment "${words[1]}"
      [[ "${act}" == nas-reachable ]] && exit 0
      gate_segment "${words[2]}"
      local dir="${nas}/${words[2]}"
      case "${act}" in
        nas-dir) [[ -d "${dir}" ]]; exit $? ;;
        nas-size) du -sk "${dir}" 2>/dev/null | cut -f1 ;;
        nas-probe)  # the listing deploy/cluster/tbrc-gate.py's PROBE_SCRIPT prints
          local p="${dir}"
          if [ ! -d "$p" ]; then echo 'MISSING'; exit 0; fi
          echo 'EXISTS'
          echo "ALL:$(find "$p" -maxdepth 1 -type f ! -name '.*' 2>/dev/null | wc -l)"
          echo "FQ:$(find "$p" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) 2>/dev/null | wc -l)"
          echo "EMPTY:$(find "$p" -maxdepth 1 -type f -size 0 2>/dev/null | wc -l)"
          echo "BYTES:$(du -sb "$p" 2>/dev/null | cut -f1)"
          find "$p" -maxdepth 1 -type f \( -name '*.fastq*' -o -name '*.fq*' \) -printf 'F:%f:%s\n' 2>/dev/null | head -20
          find "$p" -maxdepth 1 -type f ! -name '*.fastq*' ! -name '*.fq*' ! -name '.*' -printf 'O:%f\n' 2>/dev/null | head -6 ;;
        nas-integrity)
          [[ "${words[3]}" =~ ^[1-9][0-9]{0,2}$ && "${words[3]}" -le 120 && "${words[4]}" =~ ^[1-9][0-9]{0,2}$ && "${words[4]}" -le 120 ]] || gate_refuse
          local ok=0 f
          for f in "${dir}"/*.gz; do
            [ -e "$f" ] || continue
            if gzip -t -- "$f" 2>/dev/null; then ok=$(( ok + 1 )); else echo "BAD $(basename "$f")"; fi
          done
          echo "OK ${ok}"; echo "SKIPPED 0" ;;
        nas-trim-hits)
          [[ "${words[3]}" =~ ^[ACGTRYSWKMBDHVN]{12}$ ]] || gate_refuse
          local f; f="$(find "${dir}" -maxdepth 1 -type f \( -name '*.fastq*' -o -name '*.fq*' \) 2>/dev/null | head -1)"
          [ -n "$f" ] || exit 0
          case "$f" in *.gz) zcat -- "$f" ;; *) cat -- "$f" ;; esac | head -80000 | awk 'NR%4==2' | grep -c -- "${words[3]}" || true ;;
      esac ;;
  esac
  exit 0
}
SH
echo 100000000 > "${T}/free_kb"
# ssh: storage probes succeed, a missing raw folder fails "test -d", cluster
# commands are logged, and a poll for the detached driver reports exit 0. The
# words after the host go to the stub gate first.
cat > "${T}/bin/ssh" <<SH
#!/bin/bash
echo "ssh \$*" >> "${T}/ssh.log"
T_DIR="${T}"; . "${T}/gate_stub.sh"
a=("\$@"); i=0
while [[ "\${a[i]:-}" == -* ]]; do case "\${a[i]}" in -o|-p|-l|-i) i=\$((i + 2)) ;; *) i=\$((i + 1)) ;; esac; done
gate_stub "\${a[@]:i+1}"
cmd="\${*: -1}"
case "\$cmd" in
  *"/raw/"*"test -d"*|*"test -d"*"/raw/"*) f=\$(printf '%s' "\$cmd" | grep -o '/raw/[A-Za-z0-9._-]*' | head -1); [[ -d "${T}/nas\$f" ]] ;;
  *__TBRC_EXIT__*|*".exit"*)
    . "${T}/poll_answer.sh"
    key="\$(printf '%s' "\$cmd" | grep -o '[.]parallel/[A-Za-z0-9._-]*' | head -1 | sed 's#.*/##')"
    poll_status "\$key" || exit \$?
    poll_answer "\$key" \\
                "\$(printf '%s' "\$cmd" | grep -o 'tail -n +[0-9]*' | head -1 | grep -o '[0-9]*\$' | awk '{print \$1 - 1}')" ;;
  *) exit 0 ;;
esac
SH
# rsync: host:/raw/... is served from the fake NAS; other remote transfers are
# logged; local-to-local copies are real.
cat > "${T}/bin/rsync" <<SH
#!/bin/bash
echo "rsync \$*" >> "${T}/rsync.log"
args=(); for a in "\$@"; do
  case "\$a" in tester@nas.invalid:*) args+=("${T}/nas\${a#tester@nas.invalid:}") ;; *) args+=("\$a") ;; esac
done
for a in "\${args[@]}"; do [[ "\$a" == *@*:* || "\$a" == *.invalid:* ]] && exit 0; done
clean=(); skip=0; for a in "\${args[@]}"; do
  if (( skip )); then skip=0; continue; fi
  [[ "\$a" == -e ]] && { skip=1; continue; }; [[ "\$a" == -e* ]] && continue; clean+=("\$a")
done
exec "${REAL_RSYNC}" "\${clean[@]}"
SH
chmod +x "${T}/bin/"*
export PATH="${T}/bin:${PATH}" TBRC_RUNS_DIR="${T}/runs" TBRC_POLL_SECONDS=0

fastq_pair() { mkdir -p "$1"; printf '@r\nACGT\n+\nIIII\n' | gzip > "$1/$2_S1_R1_001.fastq.gz"; printf '@r\nACGT\n+\nIIII\n' | gzip > "$1/$2_S1_R2_001.fastq.gz"; }
sample_json() { mkdir -p "$1"; printf '{"sample_name": "%s"}\n' "$2" > "$1/sample.json"; }
check_local() {  # label, output of the launcher, expected run dir under backend/
  local label="$1" out="$2" want="$3"
  grep -q "RUN_ISOLATED .*table=data/runs/${want}/sample.json input=data/runs/${want}" <<< "$out" && pass "$label: workflow pointed at data/runs/${want}" || fail "$label: workflow not pointed at data/runs/${want}: $(grep RUN_ISOLATED <<< "$out")"
  grep -q '^TABLE_OK$' <<< "$out" && pass "$label: sample.json in place" || fail "$label: $(grep TABLE_ <<< "$out")"
  grep -q '^READS_OK$' <<< "$out" && pass "$label: one R1 and one R2 in place" || fail "$label: $(grep READS_ <<< "$out")"
  # backend/snakemake.sh submits with sbatch unless told the jobs run here
  grep -q '^EXECUTOR=local$' <<< "$out" && pass "$label: jobs run on this host, not through sbatch" || fail "$label: $(grep EXECUTOR <<< "$out")"
}

# The app stages an upload in runs/<user>/<run>/ beside sample.json, the run name
# being the read name with _S1/_R1/_001 and the extension removed.
echo "== upload + local"
U="${T}/runs/alice/plateA"; sample_json "$U" plateA; fastq_pair "$U" plateA
out="$(bash "${T}/repo/frontend/scripts/goUpload.bash" -u alice -r plateA -l \
        -d "$U/plateA_S1_R1_001.fastq.gz" -e "$U/plateA_S1_R2_001.fastq.gz" \
        -n plateA_S1_R1_001.fastq.gz -m plateA_S1_R2_001.fastq.gz 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
check_local "upload+local" "$out" alice/plateA
grep -q '^TBRC_INPUT_BYTES=[1-9]' <<< "$out" && pass "input size reported" || fail "no TBRC_INPUT_BYTES"
grep -q 'dtn.invalid\|login.invalid' "${T}/ssh.log" "${T}/rsync.log" 2>/dev/null && fail "upload+local contacted the cluster" || pass "cluster never contacted"

echo "== upload + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
U="${T}/runs/alice/plateB"; sample_json "$U" plateB; fastq_pair "$U" plateB
out="$(bash "${T}/repo/frontend/scripts/goUpload.bash" -u alice -r plateB \
        -d "$U/plateB_S1_R1_001.fastq.gz" -e "$U/plateB_S1_R2_001.fastq.gz" \
        -n plateB_S1_R1_001.fastq.gz -m plateB_S1_R2_001.fastq.gz 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
d="$(grep -o 'run_isolated.sh [^>]*' "${T}/ssh.log" | head -1)"
grep -Eq "alice-plateB[[:space:]]+data/runs/alice/plateB/sample.json[[:space:]]+data/runs/alice/plateB" <<< "$d" && pass "cluster told data/runs/alice/plateB" || fail "cluster dispatch: ${d:-none}"
grep -q "${T}/runs/alice tester@dtn.invalid:/cluster/TBRC/data/runs/" "${T}/rsync.log" && pass "user folder shipped to the cluster" || fail "ship: $(cat "${T}/rsync.log")"
[[ -f "$U/sample.json" && -f "$U/plateB_S1_R1_001.fastq.gz" ]] && pass "shipped folder holds sample.json and the reads" || fail "shipped folder incomplete"

echo "== storage server + local"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
fastq_pair "${T}/nas/raw/plateC" plateC; sample_json "${T}/runs/bob/plateC" plateC
out="$(bash "${T}/repo/frontend/scripts/goServer.bash" -s nas1 -u bob -f plateC -n plateC -l 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
check_local "storage+local" "$out" bob/plateC
grep -q '^TBRC_INPUT_BYTES=[1-9]' <<< "$out" && pass "input size reported" || fail "no TBRC_INPUT_BYTES"
grep -q 'dtn.invalid\|login.invalid' "${T}/ssh.log" "${T}/rsync.log" && fail "storage+local contacted the cluster" || pass "cluster never contacted"
out="$(bash "${T}/repo/frontend/scripts/goServer.bash" -s nas1 -u bob -f noSuchPlate -n noSuchPlate -l 2>&1)"; rc=$?
[[ $rc != 0 ]] && grep -q 'Folder not found on nas1' <<< "$out" && pass "missing folder: refused with the app's wording" || fail "missing folder: exit $rc: $(tail -2 <<< "$out")"

echo "== storage server + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
sample_json "${T}/runs/bob/plateD" plateD; fastq_pair "${T}/nas/raw/plateD" plateD
out="$(bash "${T}/repo/frontend/scripts/goServer.bash" -s nas1 -u bob -f plateD -n plateD 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
grep -q "ssh tester@dtn.invalid mkdir -p /cluster/TBRC/data/runs/bob && rsync" "${T}/ssh.log" && pass "DTN pulls the raw folder" || fail "DTN pull: $(head -3 "${T}/ssh.log")"
grep -q "${T}/runs/bob/plateD/ tester@dtn.invalid:/cluster/TBRC/data/runs/bob/plateD/" "${T}/rsync.log" && pass "run folder shipped" || fail "ship: $(cat "${T}/rsync.log")"
grep -Eq "bob-plateD[[:space:]]+data/runs/bob/plateD/sample.json[[:space:]]+data/runs/bob/plateD/" "${T}/ssh.log" && pass "cluster told data/runs/bob/plateD" || fail "dispatch: $(grep -o 'run_isolated.sh [^>]*' "${T}/ssh.log")"
grep -q 'TBRC_EXECUTOR' "${T}/ssh.log" && fail "a cluster dispatch carries TBRC_EXECUTOR" || pass "cluster dispatch keeps sbatch (no TBRC_EXECUTOR)"

# updatebarcodes.bash has no section of its own above; it must put new sets where the cluster keeps them.
echo "== barcodes + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"; mkdir -p "${T}/bc/setX"; echo x > "${T}/bc/setX/BC1.txt"
TBRC_BC_DIR="${T}/bc" bash "${T}/repo/frontend/scripts/updatebarcodes.bash" > /dev/null 2>&1; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc"
[[ "$(cat "${T}/rsync.log")" == "rsync -av --ignore-existing --exclude=bc_names.txt ${T}/bc/ tester@dtn.invalid:/cluster/TBRC/data/barcodes/" ]] \
  && pass "new sets pushed under the pipeline root, never over existing ones" || fail "barcodes: $(cat "${T}/rsync.log")"

# ---- with the cluster gate on --------------------------------------------------
# Same launchers, config with cluster.use_gate true. The cluster's login and transfer
# nodes then accept only the named actions in deploy/cluster/README.md, each with bare-word
# arguments, and an rsync push in the three shapes the gate's receive rule knows.
echo
echo "== gate mode: the sections above ran without use_gate"
jq -e '.cluster | has("use_gate") | not' "${T}/repo/frontend/server_config.json" > /dev/null \
  && pass "use_gate absent from the config the four sections used" || fail "use_gate present in the config of the four sections"
cp -R "${T}/repo" "${T}/repo_gate"
jq '.cluster.use_gate = true' "${T}/repo/frontend/server_config.json" > "${T}/repo_gate/frontend/server_config.json"
G="${T}/repo_gate/frontend/scripts"
GATE_ACTIONS=" launch poll free-space run-bytes results-exists pull-raw nas-reachable nas-dir nas-size nas-probe nas-integrity nas-trim-hits "
# The ssh calls to the cluster, one per line as "host action arguments" (attempt ids masked).
gate_calls() {
  local line w i
  while IFS= read -r line; do
    [[ "${line}" == *dtn.invalid* || "${line}" == *login.invalid* ]] || continue
    read -ra w <<< "${line}"; i=1
    while [[ "${w[i]:-}" == -* ]]; do case "${w[i]}" in -o|-p) i=$((i + 2)) ;; *) i=$((i + 1)) ;; esac; done
    echo "${w[*]:i}"
  done < "${T}/ssh.log" | sed -E 's/[0-9]{8}-[0-9]{6}-[0-9]+/ATTEMPT/'
}
# Every ssh to the cluster begins with a named action, and nothing names the storage server.
gate_log_check() {
  local label="$1" line n=0 bad="" act
  while IFS= read -r line; do
    n=$((n + 1)); act="$(cut -d' ' -f2 <<< "${line}")"
    [[ "${GATE_ACTIONS}" == *" ${act} "* ]] || bad+=" [${line}]"
  done < <(gate_calls)
  (( n > 0 )) && [[ -z "${bad}" ]] && pass "$label: all ${n} ssh calls to the cluster name an action" || fail "$label: ${n} ssh calls, not actions:${bad}"
  [[ "$(grep -c '^ssh' "${T}/ssh.log")" == "${n}" ]] && pass "$label: no other ssh was made" || fail "$label: ssh calls to somewhere else: $(grep -v 'dtn.invalid\|login.invalid' "${T}/ssh.log")"
  grep -q 'nas\.invalid' "${T}/ssh.log" "${T}/rsync.log" && fail "$label: something named the storage server" || pass "$label: nothing names the storage server"
}

echo "== gate mode: storage server + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
fastq_pair "${T}/nas/raw/plateG" plateG; sample_json "${T}/runs/bob/plateG" plateG
out="$(bash "${G}/goServer.bash" -s nas1 -u bob -f plateG -n plateG 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
gate_log_check "storage+cluster"
want="tester@dtn.invalid nas-dir nas1 plateG
tester@dtn.invalid nas-size nas1 plateG
tester@dtn.invalid free-space
tester@dtn.invalid pull-raw nas1 bob plateG
tester@dtn.invalid run-bytes bob plateG
login.invalid launch bob plateG ATTEMPT
login.invalid poll bob plateG ATTEMPT 0"
[[ "$(gate_calls)" == "${want}" ]] && pass "the calls, in order, with bare-word arguments" || fail "calls: $(gate_calls)"
[[ "$(cat "${T}/rsync.log")" == "rsync -a -P ${T}/runs/bob/plateG/ tester@dtn.invalid:data/runs/bob/plateG/" ]] \
  && pass "run folder pushed to a relative destination, flags exactly -a -P" || fail "push: $(cat "${T}/rsync.log")"
grep -q '^Space check: 0 MB input' <<< "$out" && pass "space check ran" || fail "no space check: $(grep -i space <<< "$out")"
grep -q '^TBRC_INPUT_BYTES=4242$' <<< "$out" && pass "input size reported" || fail "no TBRC_INPUT_BYTES=4242"
grep -q '^Run bob-plateG started' <<< "$out" && pass "launch answered" || fail "no launch line"
grep -q 'TBRC_EXECUTOR' "${T}/ssh.log" && fail "a cluster dispatch carries TBRC_EXECUTOR" || pass "cluster dispatch keeps sbatch (no TBRC_EXECUTOR)"

echo "== gate mode: missing folder, too little room, local run"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
out="$(bash "${G}/goServer.bash" -s nas1 -u bob -f noSuchPlate -n noSuchPlate 2>&1)"; rc=$?
[[ $rc != 0 ]] && grep -q 'Folder not found on nas1' <<< "$out" && pass "missing folder: refused with the app's wording" || fail "missing folder: exit $rc: $(tail -2 <<< "$out")"
[[ "$(gate_calls)" == "tester@dtn.invalid nas-dir nas1 noSuchPlate" && ! -s "${T}/rsync.log" ]] && pass "missing folder: nothing pulled, pushed or started" || fail "after the miss: $(gate_calls) / $(cat "${T}/rsync.log")"
: > "${T}/ssh.log"; echo 10 > "${T}/free_kb"
out="$(bash "${G}/goServer.bash" -s nas1 -u bob -f plateG -n plateG 2>&1)"; rc=$?
echo 100000000 > "${T}/free_kb"
[[ $rc == 1 ]] && grep -q 'Not enough space on the cluster to run plateG' <<< "$out" && pass "too little room: refused with the app's wording" || fail "no room: exit $rc: $(tail -3 <<< "$out")"
[[ "$(gate_calls | cut -d' ' -f2 | paste -sd' ' -)" == "nas-dir nas-size free-space" ]] && pass "too little room: stopped before the pull" || fail "no room: $(gate_calls)"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
out="$(bash "${G}/goServer.bash" -s nas1 -u bob -f plateG -n plateG -l 2>&1)"; rc=$?
[[ $rc != 0 ]] && grep -q 'Cannot run locally' <<< "$out" && pass "local run from a storage server: refused, this host does not reach it" || fail "local: exit $rc: $(tail -2 <<< "$out")"
[[ ! -s "${T}/ssh.log" && ! -s "${T}/rsync.log" ]] && pass "local run: no ssh, no rsync" || fail "local run contacted: $(cat "${T}/ssh.log" "${T}/rsync.log")"

echo "== gate mode: upload + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
U="${T}/runs/alice/plateU"; sample_json "$U" plateU; fastq_pair "$U" plateU
out="$(bash "${G}/goUpload.bash" -u alice -r plateU \
        -d "$U/plateU_S1_R1_001.fastq.gz" -e "$U/plateU_S1_R2_001.fastq.gz" \
        -n plateU_S1_R1_001.fastq.gz -m plateU_S1_R2_001.fastq.gz 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
gate_log_check "upload+cluster"
[[ "$(gate_calls)" == "login.invalid launch alice plateU ATTEMPT
login.invalid poll alice plateU ATTEMPT 0" ]] && pass "launch and poll for data/runs/alice/plateU" || fail "calls: $(gate_calls)"
[[ "$(cat "${T}/rsync.log")" == "rsync -a -P ${T}/runs/alice tester@dtn.invalid:data/runs/" ]] \
  && pass "user folder pushed to data/runs/, flags exactly -a -P" || fail "push: $(cat "${T}/rsync.log")"

echo "== gate mode: barcodes + cluster"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
TBRC_BC_DIR="${T}/bc" bash "${G}/updatebarcodes.bash" > /dev/null 2>&1; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc"
[[ "$(cat "${T}/rsync.log")" == "rsync -av --ignore-existing --exclude=bc_names.txt ${T}/bc/ tester@dtn.invalid:data/barcodes/" ]] \
  && pass "new sets pushed to data/barcodes/, flags exactly -av --ignore-existing" || fail "barcodes: $(cat "${T}/rsync.log")"
[[ ! -s "${T}/ssh.log" ]] && pass "no ssh for a push" || fail "ssh: $(cat "${T}/ssh.log")"

echo "== gate mode: preflight"
P="${T}/nas/raw/plateP"; mkdir -p "$P"
for r in R1 R2; do printf '@r\nGGGAATTCGAGGTGCAGCTG\n+\nIIIIIIIIIIIIIIIIIIII\n' | gzip > "$P/plateP_S1_${r}_001.fastq.gz"; done
preflight() { bash "${G}/preflight.bash" -s nas1 -p tiago_vsearch -t GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG -n plateP -u bob "$@" 2>&1; }
codes() { cut -d'|' -f2 <<< "$1" | paste -sd' ' -; }
: > "${T}/ssh.log"; : > "${T}/rsync.log"
out="$(preflight -i plateP)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
[[ "$(codes "$out")" == "tiny_input trim_rare checked" ]] && pass "the usual findings, from the probe, the integrity test and the trim count" || fail "findings: $(codes "$out")"
grep -q '^OK|checked|Checked 2 FASTQ file(s), 0 MB, on nas1\.|$' <<< "$out" && pass "closing line unchanged" || fail "closing line: $(tail -1 <<< "$out")"
gate_log_check "preflight"
want="tester@dtn.invalid nas-reachable nas1
tester@dtn.invalid nas-probe nas1 plateP
tester@dtn.invalid nas-integrity nas1 plateP 20 8
tester@dtn.invalid results-exists bob plateP plateP
tester@dtn.invalid nas-trim-hits nas1 plateP GGGAATTCGAGG"
[[ "$(gate_calls)" == "${want}" ]] && pass "the calls, in order, with bare-word arguments" || fail "calls: $(gate_calls)"
mkdir -p "${T}/cluster_results/bob/plateP/plateP"
out="$(preflight -i plateP)"
[[ "$(codes "$out")" == "tiny_input overwrite trim_rare checked" ]] && pass "an existing result is reported" || fail "findings: $(codes "$out")"
: > "${T}/ssh.log"
out="$(preflight -i noSuchPlate)"; rc=$?
[[ $rc == 0 && "$(codes "$out")" == "folder_missing" ]] && grep -q '^FAIL|folder_missing|No folder called' <<< "$out" && pass "missing folder: FAIL folder_missing" || fail "missing: exit $rc: $out"
[[ "$(gate_calls | cut -d' ' -f2 | paste -sd' ' -)" == "nas-reachable nas-probe" ]] && pass "missing folder: stopped after the probe" || fail "missing: $(gate_calls)"

# A value the gate refuses is a failed remote call, as an ssh that fails is today: the
# check ends without a verdict (the app then says the check could not run), and the
# closing OK line is never printed.
echo "== gate mode: a refused value is a failed check, not a pass"
: > "${T}/ssh.log"
out="$(preflight -i plateP -t GGGAATTCGAGUGCAGCTG)"; rc=$?
[[ $rc != 0 ]] && ! grep -q '^OK|checked' <<< "$out" && pass "trim sequence the gate refuses: no verdict, no OK line" || fail "refused trim: exit $rc, $(codes "$out")"
grep -q 'nas-trim-hits nas1 plateP GGGAATTCGAGU$' "${T}/ssh.log" && pass "(the gate, not this script, did the refusing)" || fail "trim call not made: $(tail -2 "${T}/ssh.log")"
: > "${T}/ssh.log"
out="$(preflight -i plateP -t 'GGGAATTC GAGGTGCAG')"; rc=$?
[[ $rc != 0 ]] && ! grep -q '^OK|checked' <<< "$out" && pass "trim sequence with a space: no verdict, no OK line" || fail "spaced trim: exit $rc, $(codes "$out")"
grep -q 'nas-trim-hits' "${T}/ssh.log" && fail "a trim sequence with a space was sent" || pass "(it was not sent at all)"
out="$(TBRC_PREFLIGHT_GZIP_SECONDS=0 preflight -i plateP)"; rc=$?
[[ $rc != 0 ]] && ! grep -q '^OK|checked' <<< "$out" && pass "time budget the gate refuses: no verdict, no OK line" || fail "refused budget: exit $rc, $(codes "$out")"
# the same value, in a name the gate takes: nothing is refused
out="$(TBRC_PREFLIGHT_GZIP_SECONDS=30 TBRC_PREFLIGHT_GZIP_FILE_SECONDS=12 preflight -i plateP)"
grep -q '^OK|checked' <<< "$out" && grep -q 'nas-integrity nas1 plateP 30 12$' "${T}/ssh.log" && pass "budgets 30 and 12 pass through unchanged" || fail "budgets: $(codes "$out")"

# The loop that follows a detached run is shared by both ways of launching one.
echo "== following a run: gate and shell commands share one loop"
# (timeout: a loop that never ends is a failure here, not a hang)
loop() { timeout 20 bash -c 'source "$1"; set +e; shift; "$@"; echo "rc=$?"' _ "${T}/repo/frontend/scripts/common.sh" "$@" 2>&1; }
check_follow() {  # label, then the two ways of starting a run named alice-<kind>
  local label="$1" kind out
  for kind in chatty failing vanished; do
    : > "${T}/ssh.log"
    if [[ "$label" == gate ]]; then out="$(loop dispatch_via_gate login.invalid alice "$kind")"
    else out="$(loop dispatch_detached login.invalid /cluster/TBRC /usr/bin snakemake "alice-$kind" "data/runs/alice/$kind/sample.json" "data/runs/alice/$kind/")"; fi
    case "$kind" in
      chatty) [[ "$(grep -c '^step ' <<< "$out")" == 3 && "$(grep -c '^step one$' <<< "$out")" == 1 && "$out" == *"rc=0" && "$out" != *WRONG_SEEN* ]] \
                && pass "$label: new lines relayed once each, from where the last poll stopped, then exit 0" || fail "$label chatty: $out" ;;
      failing) [[ "$out" == *boom* && "$out" == *"rc=3" ]] && pass "$label: the run's own exit status is returned" || fail "$label failing: $out" ;;
      vanished) [[ "$out" == *"no longer running and left no exit status"* && "$out" == *"rc=70" ]] \
                && pass "$label: a driver that vanished ends the wait after three polls" || fail "$label vanished: $out" ;;
    esac
  done
}
check_follow shell
check_follow gate
: > "${T}/ssh.log"; loop dispatch_via_gate login.invalid alice chatty > /dev/null
[[ "$(gate_calls | grep poll | sed 's/.* //' | paste -sd' ' -)" == "0 2" ]] && pass "gate: the second poll asks from line 3" || fail "polls: $(gate_calls)"

# ---- who turns the gate on, and off ---------------------------------------------
# cluster.use_gate must be exactly true. TBRC_USE_GATE=0 in the environment turns it off
# whatever the config says (the tbrc command line goes through the user's own, ungated key);
# no other value of the variable turns anything on or off.
echo
echo "== gate switch: use_gate in the config, TBRC_USE_GATE in the environment"
gate_state() {  # USE_GATE_JSON ("-" = key absent, "broken" = unparseable file), then NAME=value words: prints on or off
  local lit="$1" cfg="${T}/cfg_state.json"; shift
  case "$lit" in
    -) jq 'del(.cluster.use_gate)' "${T}/repo_gate/frontend/server_config.json" > "$cfg" ;;
    broken) printf '{ "cluster": { "use_gate": tru' > "$cfg" ;;
    *) jq --argjson v "$lit" '.cluster.use_gate = $v' "${T}/repo_gate/frontend/server_config.json" > "$cfg" ;;
  esac
  env "$@" bash -c 'source "$1"; CONFIG_FILE="$2"; if gate_enabled; then echo on; else echo off; fi' _ "${G}/common.sh" "$cfg"
}
expect_state() {  # WANT, USE_GATE_JSON, env words
  local want="$1" lit="$2"; shift 2
  [[ "$(gate_state "$lit" "$@")" == "$want" ]] && pass "use_gate ${lit}${*:+, $*}: ${want}" || fail "use_gate ${lit}${*:+, $*}: wanted ${want}, got $(gate_state "$lit" "$@")"
}
expect_state on true
expect_state off false
expect_state off -
for lit in null '"true"' 1 '"yes"' '[true]' '{}' broken; do expect_state off "$lit"; done
expect_state off true TBRC_USE_GATE=0
for v in 1 "" 00 false " 0" "0 " yes; do expect_state on true "TBRC_USE_GATE=$v"; done
expect_state off false TBRC_USE_GATE=1
expect_state off - TBRC_USE_GATE=1

echo "== TBRC_USE_GATE=0: a config with use_gate true runs the commands it always ran"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
out="$(TBRC_USE_GATE=0 bash "${G}/goServer.bash" -s nas1 -u bob -f plateD -n plateD 2>&1)"; rc=$?
[[ $rc == 0 ]] && pass "exit 0" || fail "exit $rc: $(tail -3 <<< "$out")"
grep -q "ssh tester@dtn.invalid mkdir -p /cluster/TBRC/data/runs/bob && rsync" "${T}/ssh.log" && pass "DTN pulls the raw folder" || fail "DTN pull: $(head -3 "${T}/ssh.log")"
grep -q "${T}/runs/bob/plateD/ tester@dtn.invalid:/cluster/TBRC/data/runs/bob/plateD/" "${T}/rsync.log" && pass "run folder shipped to the absolute path" || fail "ship: $(cat "${T}/rsync.log")"
grep -Eq "bob-plateD[[:space:]]+data/runs/bob/plateD/sample.json[[:space:]]+data/runs/bob/plateD/" "${T}/ssh.log" && pass "cluster told data/runs/bob/plateD" || fail "dispatch: $(grep -o 'run_isolated.sh [^>]*' "${T}/ssh.log")"
grep -Eq ' (nas-dir|nas-size|free-space|pull-raw|run-bytes|launch|poll) ' "${T}/ssh.log" && fail "a gate action was asked for" || pass "no gate action asked for"
: > "${T}/ssh.log"; : > "${T}/rsync.log"
TBRC_USE_GATE=0 bash "${G}/goUpload.bash" -u alice -r plateB \
  -d "${T}/runs/alice/plateB/plateB_S1_R1_001.fastq.gz" -e "${T}/runs/alice/plateB/plateB_S1_R2_001.fastq.gz" \
  -n plateB_S1_R1_001.fastq.gz -m plateB_S1_R2_001.fastq.gz > /dev/null 2>&1
grep -q "${T}/runs/alice tester@dtn.invalid:/cluster/TBRC/data/runs/" "${T}/rsync.log" && grep -q "run_isolated.sh" "${T}/ssh.log" \
  && pass "upload: absolute push and composed launch" || fail "upload: $(cat "${T}/rsync.log") / $(head -2 "${T}/ssh.log")"
: > "${T}/rsync.log"
TBRC_USE_GATE=0 TBRC_BC_DIR="${T}/bc" bash "${G}/updatebarcodes.bash" > /dev/null 2>&1
[[ "$(cat "${T}/rsync.log")" == "rsync -av --ignore-existing --exclude=bc_names.txt ${T}/bc/ tester@dtn.invalid:/cluster/TBRC/data/barcodes/" ]] \
  && pass "barcodes: absolute target" || fail "barcodes: $(cat "${T}/rsync.log")"
out="$(TBRC_USE_GATE=0 bash "${G}/goServer.bash" -s nas1 -u bob -f plateC -n plateC -l 2>&1)"; rc=$?
[[ $rc == 0 ]] && grep -q '^READS_OK$' <<< "$out" && pass "run on this host from a storage server is allowed again" || fail "local: exit $rc: $(tail -2 <<< "$out")"
[[ "$(env TBRC_USE_GATE=0 bash -c 'source "$1"; cluster_path data/runs/' _ "${G}/common.sh")" == "/cluster/TBRC/data/runs/" ]] \
  && pass "cluster_path gives the configured root" || fail "cluster_path"

# tbrc run --server hands the work to goServer.bash through the lab's own ssh key, which the
# gate does not cover and must not: it has to ask for the commands, not for gate actions.
echo "== tbrc run --server turns the gate off for goServer.bash"
TBRC_PY="$(cd "${SRC}/../.." && pwd)/tbrc"
if [[ -f "${TBRC_PY}" ]] && command -v python3 > /dev/null; then
  mkdir -p "${T}/cli/frontend/scripts" "${T}/cli/backend"
  cp "${TBRC_PY}" "${T}/cli/tbrc"
  cp "${T}/repo_gate/frontend/server_config.json" "${T}/cli/frontend/server_config.json"
  cat > "${T}/cli/frontend/scripts/goServer.bash" <<'SH'
#!/bin/bash
echo "STUB_ARGS $*"
echo "STUB_USE_GATE ${TBRC_USE_GATE-unset}"
SH
  for inherited in unset 1; do
    if [[ "$inherited" == unset ]]; then env_cmd=(env -u TBRC_USE_GATE); else env_cmd=(env TBRC_USE_GATE=1); fi
    out="$("${env_cmd[@]}" python3 "${T}/cli/tbrc" run --user bob --server nas1 --folder "plateK${inherited}" --run-name "plateK${inherited}" \
            --email a@b.org --primer-set B2 --no-trim-primers 2>&1)"; rc=$?
    [[ $rc == 0 ]] && grep -q '^STUB_ARGS -u bob -s nas1 -f plateK' <<< "$out" && pass "tbrc ran goServer.bash (TBRC_USE_GATE ${inherited} in its environment)" || fail "tbrc: exit $rc: $(tail -3 <<< "$out")"
    grep -q '^STUB_USE_GATE 0$' <<< "$out" && pass "goServer.bash got TBRC_USE_GATE=0 (TBRC_USE_GATE ${inherited} in tbrc's environment)" || fail "goServer.bash got: $(grep STUB_USE_GATE <<< "$out")"
  done
else
  echo "  SKIP  tbrc or python3 not found"
fi

# ---- a gate that refuses the status check for good --------------------------------
# A refused poll (status 3) or a gate being reinstalled (4) will not get better by
# asking again, but a network blip (255) will. Only in gate mode, consecutive refusals
# end the wait: status 70, and a message saying the job may still be running.
echo
echo "== gate mode: a gate that keeps refusing the status check"
refusals() { grep -c "poll alice $1 " "${T}/ssh.log"; }
for kind in refused refusedfour; do
  : > "${T}/ssh.log"
  out="$(TBRC_GATE_POLL_REFUSED_LIMIT=4 loop dispatch_via_gate login.invalid alice "$kind")"
  [[ "$out" == *"rc=70" ]] && pass "$kind: the wait ends with status 70" || fail "$kind: $out"
  [[ "$(refusals "$kind")" == 4 ]] && pass "$kind: after exactly 4 refused polls" || fail "$kind: $(refusals "$kind") polls"
  grep -q 'refused the status check for alice-'"$kind" <<< "$out" && grep -q 'may still be running on the cluster' <<< "$out" \
    && pass "$kind: says the gate refused and the job may still be running" || fail "$kind message: $out"
done
: > "${T}/ssh.log"
out="$(TBRC_GATE_POLL_REFUSED_LIMIT=3 loop dispatch_via_gate login.invalid alice flaky)"
[[ "$out" == *"rc=0" && "$out" == *"flaky line"* && "$(refusals flaky)" == 6 ]] \
  && pass "an answer in between resets the count (4 refusals in all, never 3 in a row)" || fail "flaky: $out"
: > "${T}/ssh.log"
out="$(TBRC_GATE_POLL_REFUSED_LIMIT=3 loop dispatch_via_gate login.invalid alice netdown)"
[[ "$out" == *"rc=0" && "$(refusals netdown)" == 7 && "$(grep -c 'Could not reach' <<< "$out")" == 6 ]] \
  && pass "status 255 keeps retrying past the limit" || fail "netdown: $out"
: > "${T}/ssh.log"
out="$(TBRC_GATE_POLL_REFUSED_LIMIT=3 loop dispatch_via_gate login.invalid alice slowgate)"
[[ "$out" == *"rc=70" && "$(refusals slowgate)" == 3 ]] && pass "slowgate: given up on at the limit" || fail "slowgate (gate): $out"
: > "${T}/ssh.log"
out="$(TBRC_GATE_POLL_REFUSED_LIMIT=3 loop dispatch_detached login.invalid /cluster/TBRC /usr/bin snakemake carol-slowgate data/runs/carol/slowgate/sample.json data/runs/carol/slowgate/)"
[[ "$out" == *"rc=0" && "$(grep -c 'Could not reach' <<< "$out")" == 5 ]] \
  && pass "the same five failed polls through the shell commands are retried, never given up on" || fail "slowgate (shell): $out"
# the limit follows the poll interval: about ten minutes of refusals, whatever the interval
limit() { env TBRC_POLL_SECONDS="$1" bash -c 'source "$1"; gate_poll_refused_limit' _ "${G}/common.sh"; }
[[ "$(limit 60) $(limit 30) $(limit 45) $(limit 600) $(limit 1000) $(limit 0)" == "10 20 14 1 1 600" ]] \
  && pass "default limit: 10 polls at 60 s, 20 at 30 s, 14 at 45 s, 1 at 600 s or more (and 600 for no pause)" || fail "limits: $(limit 60) $(limit 30) $(limit 45) $(limit 600) $(limit 1000) $(limit 0)"


# dispatch_detached composes one remote command. Its closing message used to put the run key
# inside double quotes, where $( ) and backticks still run on the cluster; it is quoted with
# shq now, like every other value in that command. The test plays the cluster's shell: it runs
# the command ssh was given, in a throwaway pipeline root, with a key that carries command
# substitutions and a double quote.
echo
echo "== shell launch: the run key in the closing message is quoted"
launch_root="${T}/launch_root"; mkdir -p "${launch_root}"
printf '#!/bin/sh\nexit 0\n' > "${launch_root}/run_isolated.sh"; chmod +x "${launch_root}/run_isolated.sh"
printf '#!/bin/sh\nexit 0\n' > "${T}/bin/conda"; chmod +x "${T}/bin/conda"
for hostile_key in 'k$(touch pwned_a)x`touch pwned_b`' 'k"; touch pwned_a; echo "y'; do
  : > "${T}/ssh.log"
  loop dispatch_detached login.invalid "${launch_root}" "${T}/bin:/usr/bin:/bin" snakemake "${hostile_key}" t.json i/ > /dev/null
  remote_cmd="$(head -1 "${T}/ssh.log")"; remote_cmd="${remote_cmd#ssh login.invalid }"
  launch_out="$(bash -c "${remote_cmd}" 2>&1)"
  [[ "${launch_out}" == "Run ${hostile_key} started on login.invalid, attempt "* ]] \
    && pass "the closing message shows the key exactly as given" || fail "closing message: ${launch_out}"
  [[ ! -e "${launch_root}/pwned_a" && ! -e "${launch_root}/pwned_b" ]] \
    && pass "nothing in the key ran" || fail "a command in the key ran: $(ls "${launch_root}" | grep -c pwned) file(s)"
  rm -f "${launch_root}"/pwned_a "${launch_root}"/pwned_b
done

echo
if (( FAILS )); then echo "${FAILS} check(s) failed"; exit 1; fi
echo "all checks passed"
