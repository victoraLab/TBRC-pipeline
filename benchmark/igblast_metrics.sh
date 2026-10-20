#!/bin/bash
# Annotate each pipeline's delivered FASTA with IgBlast and report the biological
# metrics: how many delivered sequences are a productive rearrangement, how many
# get a V and J call, and how many yield a CDR3.
#
# This is the evidence the sequence-composition metrics cannot give. A pipeline
# can produce clean-looking FASTA that IgBlast cannot annotate.
set -uo pipefail

BENCH="${BENCH_WORK:?set BENCH_WORK}"
IGROOT="${IG_ROOT:?set IG_ROOT to the extracted igblast directory}"
IMAGE="${BENCH_IMAGE:-tbrc-backend:arm64}"
OUT="${BENCH}/igblast_metrics.tsv"

# Plain lookups rather than associative arrays: macOS ships bash 3.2, which has
# no `declare -A`, and this script is run from the host.
lookup() {
    case "$1" in
        RUN12)               echo "RUN12 igh Ig" ;;
        RUN13)               echo "RUN13 igh Ig" ;;
        RUN14new)            echo "RUN14 igh Ig" ;;
        RUN10)               echo "RUN10 ig_light Ig" ;;
        TRB_RUN1)    echo "TRB_RUN1_sample trb TCR" ;;
        TRB_RUN2)              echo "TRB_RUN2_sample trb TCR" ;;
        TRB_RUN3) echo "TRB_RUN3_sample trb TCR" ;;
        *) echo "" ;;
    esac
}

printf 'dataset\tpipeline\tpanel\tsequences\tv_called\tj_called\tproductive\tcdr3\n' > "${OUT}"

for ds in RUN12 RUN13 RUN14new RUN10 TRB_RUN1 TRB_RUN2 TRB_RUN3; do
  set -- $(lookup "$ds"); s="$1"; panel="$2"; st="$3"
  [ -n "$s" ] || { echo "skip $ds: unknown"; continue; }
  for pl in classic_luka tiago_vsearch; do
    fa="${BENCH}/${ds}/${pl}/results/bench/${ds}/${s}/${s}.imgt.ready.fasta"
    [ -f "${fa}" ] || fa="${BENCH}/${ds}/${pl}/results/bench/${ds}/${s}/${s}.final.fasta"
    if [ ! -f "${fa}" ]; then
      echo "skip ${ds}/${pl}: no FASTA"; continue
    fi
    echo "igblast ${ds} / ${pl} (mouse/${panel}) ..."
    docker run --rm --platform linux/arm64 \
      -v "${IGROOT}":/ig:ro -v "$(dirname "${fa}")":/data \
      -e IGDATA=/ig/ncbi-igblast-1.22.0 \
      -w /data "${IMAGE}" bash -c "
        igblastn \
          -germline_db_V /ig/work/mouse/${panel}/mouse_${panel}_V \
          -germline_db_J /ig/work/mouse/${panel}/mouse_${panel}_J \
          \$( [ -f /ig/work/mouse/${panel}/mouse_${panel}_D.nhr ] && echo '-germline_db_D /ig/work/mouse/${panel}/mouse_${panel}_D' ) \
          -organism mouse -ig_seqtype ${st} \
          -auxiliary_data /ig/ncbi-igblast-1.22.0/optional_file/mouse_gl.aux \
          -query '$(basename "${fa}")' -outfmt 19 -num_threads 4 \
          -out /data/${pl}.airr.tsv" >/dev/null 2>&1

    airr="$(dirname "${fa}")/${pl}.airr.tsv"
    if [ ! -s "${airr}" ]; then echo "  no AIRR output"; continue; fi
    python3 - "${airr}" "${ds}" "${pl}" "${panel}" >> "${OUT}" <<'PY'
import csv, sys
path, ds, pl, panel = sys.argv[1:5]
n = v = j = prod = cdr3 = 0
with open(path) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n += 1
        if (row.get("v_call") or "").strip(): v += 1
        if (row.get("j_call") or "").strip(): j += 1
        if (row.get("productive") or "").strip().upper() in ("T", "TRUE"): prod += 1
        if (row.get("cdr3_aa") or row.get("cdr3") or "").strip(): cdr3 += 1
print("\t".join(map(str, [ds, pl, panel, n, v, j, prod, cdr3])))
PY
    echo "  done"
  done
done
echo
echo "wrote ${OUT}"
