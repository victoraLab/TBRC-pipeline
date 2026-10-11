#!/bin/bash
# Two-pass barcode demultiplexing for the tiago_vsearch pipeline.
#
# Differs from the classic split_fasta.sh in four ways:
#   - works on FASTQ, so quality survives into the clustering stage
#   - mismatch tolerance is derived per barcode file rather than hardcoded
#     (--exact on BC1, default 1 on BC2), which both fixes barcode set B and
#     recovers the reads BC1's exact matching used to discard
#   - wells are filtered on read count, not file size in bytes
#   - every read that falls out is counted, not silently deleted
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

merged_fastq=""; barcode1=""; barcode2=""; destiny=""; sentinel=""
min_reads=10; keep_unmatched="false"

while getopts "f:1:2:d:o:n:u:" opt; do
    case "${opt}" in
        f) merged_fastq="${OPTARG}" ;;
        1) barcode1="${OPTARG}" ;;
        2) barcode2="${OPTARG}" ;;
        d) destiny="${OPTARG}" ;;
        o) sentinel="${OPTARG}" ;;
        n) min_reads="${OPTARG}" ;;
        u) keep_unmatched="${OPTARG}" ;;
        *) echo "Usage: $0 -f merged.fastq -1 BC1 -2 BC2 -d destdir -o sentinel [-n min_reads] [-u keep_unmatched]" >&2; exit 1 ;;
    esac
done

split1="${destiny}split1/"
split2="${destiny}split1/split2/"
qc_demux="${destiny}qc/demux/"
mkdir -p "${split1}" "${split2}" "${qc_demux}"

safety_report="${qc_demux}barcode_safety.tsv"
: > "${safety_report}"
bc1_mm="$(python3 "${SCRIPT_DIR}/bc_safety.py" "${barcode1}" --report "${safety_report}" --label BC1)"
bc2_mm="$(python3 "${SCRIPT_DIR}/bc_safety.py" "${barcode2}" --report "${safety_report}" --label BC2)"
echo "Barcode mismatch tolerance: BC1=${bc1_mm} BC2=${bc2_mm} (derived from barcode separation)"

funnel="${qc_demux}funnel.tsv"
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }

merged_reads=$(count_reads "${merged_fastq}")

# Pass 1: BC1 at the 5' end.
bc1_stats="${qc_demux}BC1.stats.txt"
fastx_barcode_splitter.pl --bcfile "${barcode1}" --bol --mismatches "${bc1_mm}" \
    --prefix "${split1}" --suffix ".fastq" < "${merged_fastq}" \
    | awk -F '\t' 'NF == 3 && $1 != "Barcode" && $1 != "total" { print }' > "${bc1_stats}"

bc1_unmatched=0
for f in "${split1}"*unmatched*; do
    bc1_unmatched=$(( bc1_unmatched + $(count_reads "$f") ))
    [[ "${keep_unmatched}" == "true" ]] || rm -f "$f"
done

# Pass 2: BC2 at the 3' end, within each BC1 bin.
bc_stats="${destiny}BC.stats.txt"
: > "${bc_stats}"
for f in "${split1}"*.fastq; do
    # The first pass writes a file for every barcode, including ones that matched
    # nothing. fastx_barcode_splitter aborts on an empty input ("unknown file
    # format. First character = ''"), so skip those rather than feeding them in.
    if [[ ! -s "$f" ]]; then
        rm -f "$f"
        continue
    fi
    x="${f##*/}"
    fastx_barcode_splitter.pl --bcfile "${barcode2}" --eol --mismatches "${bc2_mm}" \
        --prefix "${split2}${x%.*}" --suffix ".fastq" < "$f" \
        | awk -F '\t' 'NF == 3 && $1 != "Barcode" && $1 != "total" { print }' >> "${bc_stats}"
done

bc2_unmatched=0
for f in "${split2}"*unmatched*; do
    bc2_unmatched=$(( bc2_unmatched + $(count_reads "$f") ))
    [[ "${keep_unmatched}" == "true" ]] || rm -f "$f"
done

# Depth filter on read count, not bytes. A byte threshold rescales with amplicon
# length, which is exactly the kind of silent drift this pipeline is removing.
demuxed_reads=0; kept_reads=0; wells_total=0; wells_kept=0
for f in "${split2}"*.fastq; do
    wells_total=$(( wells_total + 1 ))
    n=$(count_reads "$f")
    demuxed_reads=$(( demuxed_reads + n ))
    if (( n < min_reads )); then
        rm -f "$f"
    else
        wells_kept=$(( wells_kept + 1 ))
        kept_reads=$(( kept_reads + n ))
    fi
done

{
    printf 'stage\treads\twells\n'
    printf 'merged\t%s\t\n' "${merged_reads}"
    printf 'bc1_unmatched\t%s\t\n' "${bc1_unmatched}"
    printf 'bc2_unmatched\t%s\t\n' "${bc2_unmatched}"
    printf 'demultiplexed\t%s\t%s\n' "${demuxed_reads}" "${wells_total}"
    printf 'passing_depth_filter\t%s\t%s\n' "${kept_reads}" "${wells_kept}"
} > "${funnel}"

echo "Demultiplexed ${demuxed_reads}/${merged_reads} reads into ${wells_total} wells; ${wells_kept} wells passed the ${min_reads}-read floor."

if (( wells_kept == 0 )); then
    echo "No well reached the ${min_reads}-read floor. Check that the barcode set matches this library." >&2
    exit 1
fi

touch "${sentinel}"
