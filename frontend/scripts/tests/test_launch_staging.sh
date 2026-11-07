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
# Nothing leaves this machine. The launchers run from a throwaway copy of the
# repository with a stub backend/run_isolated.sh that reports what it was given,
# stub pipeline tools, and stub ssh/rsync: a "storage server" is a local folder,
# and every cluster command is logged instead of sent.
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
# ssh: storage probes succeed, a missing raw folder fails "test -d", cluster
# commands are logged, and a poll for the detached driver reports exit 0.
cat > "${T}/bin/ssh" <<SH
#!/bin/bash
echo "ssh \$*" >> "${T}/ssh.log"
cmd="\${*: -1}"
case "\$cmd" in
  *"/raw/"*"test -d"*|*"test -d"*"/raw/"*) f=\$(printf '%s' "\$cmd" | grep -o '/raw/[A-Za-z0-9._-]*' | head -1); [[ -d "${T}/nas\$f" ]] ;;
  *__TBRC_EXIT__*|*".exit"*) echo "__TBRC_EXIT__=0" ;;
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

echo
if (( FAILS )); then echo "${FAILS} check(s) failed"; exit 1; fi
echo "all checks passed"
