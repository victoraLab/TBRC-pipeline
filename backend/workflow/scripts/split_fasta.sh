#!/bin/bash

set -euo pipefail
shopt -s nullglob

#This script is responsible for splitting fasta files into plates based on each barcode.

# Parsing command line arguments
# Defaults for the luka_light preprocessing adapters. -t is the 5' handle that
# sits between the plate barcode and the V region (same sequence the UI exposes
# as "Trim sequence"); -T is the reverse-primer tail that terminates the
# amplicon 3' of J. Both are overridable so this branch cannot drift away from
# the trim sequence the user actually selected.
luka_light_5p="GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"
luka_light_3p="CACTGACACTCCCACGAGCTC"
# Length of the gene-specific portion of the forward primer that follows the 5'
# handle in the light-chain protocol. Those bases are primer-derived, not
# template, so they are cut before collapsing. Measured on real luka_light data
# (RUN10): per-base entropy over the 23 nt after the handle averages
# 0.77 bits versus 1.62 bits from base 24 on, i.e. template diversity starts at
# base 24 exactly. fastx_trimmer -f is 1-based, hence 24.
luka_light_primer_footprint=24

while getopts "f:1:2:o:d:h:m:t:T:" opt; do
    case ${opt} in
        f)  # Input fasta file
            input_fasta="${OPTARG}" ;;
        1)  # Barcode1 file
            barcode1="${OPTARG}" ;;
        2)  # Barcode2 file
            barcode2="${OPTARG}" ;;
        o)  # Output prefix
            output="${OPTARG}" ;;
        d)  # Destination directory
            destiny="${OPTARG}" ;;
        m)  # Logical test for luka_extra step 
            method="${OPTARG}" ;;
        h)  # Number of head lines
            head="${OPTARG}" ;;
        t)  # 5' handle removed by the luka_light branch
            luka_light_5p="${OPTARG}" ;;
        T)  # 3' reverse-primer tail removed by the luka_light branch
            luka_light_3p="${OPTARG}" ;;
        *)  # In case of incorrect usage
            echo "Usage: $0 -f <input_fasta> -1 <barcode1_file> -2 <barcode2_file> -o <output_directory> -d <destination_directory> -m <method> -h <head_lines> [-t <5prime_handle>] [-T <3prime_tail>]"
            exit 1 ;;
    esac
done

# fastx_barcode_splitter tolerates 1 mismatch by default, and on a tie it takes
# whichever barcode appears FIRST in the file (its comparison is a strict "<",
# with no ambiguity check and no warning). That is only safe when every pair of
# barcodes differs in at least 2m+1 positions; correcting m errors needs that
# much separation or a corrupted barcode can sit equidistant from two of them.
#
# Barcode set B violates this: its BC2 has a minimum Hamming distance of 2 across
# 16 pairs -- P01x/P11x and P06x/P12x for every row A-H -- so a single sequencing
# error in one of those barcodes was silently and always resolved to the earlier
# plate. Plates 11 and 12 leaked into plates 1 and 6.
#
# Rather than hardcode a set name, derive the tolerance each file can actually
# support. Sets that are well separated are unaffected.
safe_mismatches() {
    local bcfile="$1"
    awk '
        NF >= 2 { bc[n++] = toupper($2) }
        END {
            if (n < 2) { print 1; exit }
            best = 999
            for (i = 0; i < n; i++)
                for (j = i + 1; j < n; j++) {
                    if (length(bc[i]) != length(bc[j])) continue
                    d = 0
                    for (k = 1; k <= length(bc[i]); k++)
                        if (substr(bc[i], k, 1) != substr(bc[j], k, 1)) d++
                    if (d < best) best = d
                }
            if (best == 999) { print 1; exit }
            # m mismatches are unambiguous only while 2m+1 <= best
            print int((best - 1) / 2)
        }
    ' "${bcfile}"
}

# Never exceed what the file supports, and never loosen what the pipeline does
# today: BC1 stays exact, BC2 keeps its 1-mismatch tolerance unless the set
# cannot support it. Recovering reads by relaxing BC1 is a separate change with
# real effect on results, and belongs with the pipeline candidates.
BC2_SAFE="$(safe_mismatches "${barcode2}")"
BC2_MISMATCHES=1
if [[ "${BC2_SAFE}" -lt 1 ]]; then
    BC2_MISMATCHES=0
    echo "Barcode set has closely spaced BC2 barcodes; demultiplexing BC2 with exact matching to avoid cross-plate misassignment." >&2
fi
echo "BC2 mismatch tolerance: ${BC2_MISMATCHES} (barcode file supports ${BC2_SAFE})"

# Defining directory structures
initial_split_folder="${destiny}split1/"
second_split_folder="${destiny}split1/split2/"
collapsed_folder="${destiny}split1/split2/collapsed/"
head_folder="${destiny}split1/split2/collapsed/head/"
final_output_folder="${destiny}split1/split2/collapsed/head/final/"

# Splitting based on Barcode1
touch ${output}
cat ${input_fasta} | fastx_barcode_splitter.pl --bcfile ${barcode1} --bol --prefix ${output} --suffix .fasta --exact
# Removing unmatched sequences from the initial split
rm ${initial_split_folder}*unmatched*
# Removing small fasta files (less than 4k in size)
find ${initial_split_folder} -name "*fasta" -size -4k -delete

initial_fastas=(${initial_split_folder}*fasta)
if [[ ${#initial_fastas[@]} -eq 0 ]]; then
    echo "No barcode-matched FASTA files remained after the first split and size filtering." >&2
    exit 1
fi

# Splitting based on Barcode2
# fastx_barcode_splitter prints a per-well summary table on stdout. It used to be
# thrown away, which is why plate_stats.R -- which needs exactly this -- never ran.
bc_stats="${destiny}BC.stats.txt"
: > "${bc_stats}"

mkdir -p ${second_split_folder}
for f in ${initial_split_folder}*fasta; do
    x="${f##*/}"
    cat $f | fastx_barcode_splitter.pl --bcfile ${barcode2} --eol --mismatches ${BC2_MISMATCHES} --prefix ${second_split_folder}${x%.*} --suffix .fasta \
        | awk -F '\t' 'NF == 3 && $1 != "Barcode" && $1 != "total" { print }' >> "${bc_stats}"
done

# Cleaning up unmatched sequences and small files
rm ${second_split_folder}*unmatched*
find ${second_split_folder} -name "*fasta" -size -4k -delete
echo ${method}

second_fastas=(${second_split_folder}*fasta)
if [[ ${#second_fastas[@]} -eq 0 ]]; then
    echo "No barcode-matched FASTA files remained after the second split and size filtering." >&2
    exit 1
fi

# If method is "luka_light", apply cutadapt and fastx_trimmer before collapsing sequences
if [ "${method}" == "luka_light" ]; then
    echo "luka_light preprocessing: 5' handle=${luka_light_5p} 3' tail=${luka_light_3p} primer footprint=$((luka_light_primer_footprint - 1)) nt"
    # Intermediates go to a scratch directory, never next to the per-well
    # FASTAs. Writing modified_*/trimmed_* into split2/ meant a re-entered run
    # picked them up in the *fasta glob, and in the deployed January build they
    # were collapsed straight into final.fasta: the delivered RUN10 result
    # carried 2,013 untrimmed duplicate records out of 3,411.
    luka_tmp="$(mktemp -d)"
    trap 'rm -rf "${luka_tmp}"' EXIT
    kept_fastas=()
    for f in "${second_fastas[@]}"; do
        x="${f##*/}"
        # Strip the plate barcode plus the 5' handle. --discard-untrimmed is
        # mandatory here: the fastx_trimmer step below cuts a FIXED number of
        # bases, so it is only correct on reads this call actually put in the
        # right register. Without it, reads carrying no handle keep their
        # barcode and instead lose 23 bases of barcode/handle, entering the
        # collapse step shifted relative to their neighbours -- on real
        # light-chain data that is 26% of reads.
        # -O 15 raises cutadapt's minimum overlap from its default of 3, so a
        # 3-base coincidence at the read start no longer counts as a hit. On
        # real data the floor costs nothing (2007/2007 kappa and 2586/2589
        # heavy sequences still match) and takes spurious hits to zero.
        cutadapt -g "${luka_light_5p}" -e 0.2 -O 15 --discard-untrimmed \
            -o "${luka_tmp}/modified_${x}" "$f"
        # A well with no handle-bearing reads leaves an empty file here, and
        # fastx_trimmer aborts on it with "Premature End-Of-File", taking the
        # whole run down under set -e. Drop the well instead.
        if [[ ! -s "${luka_tmp}/modified_${x}" ]]; then
            echo "luka_light: no reads matched the 5' handle in ${x}; dropping this well." >&2
            rm -f "$f"
            continue
        fi
        # Remove the gene-specific portion of the forward primer (see
        # luka_light_primer_footprint above); those bases are primer, not template.
        fastx_trimmer -f "${luka_light_primer_footprint}" \
            -i "${luka_tmp}/modified_${x}" -o "${luka_tmp}/trimmed_${x}"
        # Strip the reverse-primer tail. --minimum-length 1 is mandatory: when
        # the tail matches at or near position 0 (adapter dimers, short
        # inserts) cutadapt otherwise emits a zero-length FASTA record and
        # fastx_collapser aborts with "found empty nucleotide sequence",
        # killing the run under set -e.
        cutadapt -a "${luka_light_3p}" -e 0.2 -O 10 --minimum-length 1 \
            -o "${luka_tmp}/final_${x}" "${luka_tmp}/trimmed_${x}"
        if [[ ! -s "${luka_tmp}/final_${x}" ]]; then
            echo "luka_light: ${x} was emptied by the 3' tail trim; dropping this well." >&2
            rm -f "$f"
            continue
        fi
        mv "${luka_tmp}/final_${x}" "$f"
        kept_fastas+=("$f")
    done
    rm -rf "${luka_tmp}"
    trap - EXIT

    second_fastas=("${kept_fastas[@]}")
    if [[ ${#second_fastas[@]} -eq 0 ]]; then
        echo "luka_light preprocessing removed every read. The 5' handle (${luka_light_5p}) does not match this library." >&2
        exit 1
    fi
fi

# Collapsing sequences
mkdir -p ${collapsed_folder}
for f in "${second_fastas[@]}"; do
    x="${f##*/}"
    fastx_collapser -i $f -o ${collapsed_folder}${x%.*}.fasta
done

# Per-well collapse depths: one column per well, values are the read counts of
# each collapsed sequence, most abundant first. plate_stats.R turns the top-3
# share of these into a per-well quality score.
collapsed_stats="${destiny}collapsed.stats.txt"
stats_tmp="$(mktemp -d)"
for f in "${collapsed_folder}"*fasta; do
    well="${f##*/}"
    well="${well%.*}"
    {
        printf '%s\n' "${well}"
        # fastx_collapser writes headers as >rank-count
        grep '^>' "$f" 2>/dev/null | sed 's/^>[0-9]*-//' || true
    } > "${stats_tmp}/${well}.col"
done
if compgen -G "${stats_tmp}/*.col" > /dev/null; then
    paste "${stats_tmp}"/*.col > "${collapsed_stats}"
fi
rm -rf "${stats_tmp}"

# Taking the top 'head' lines of each sequence file
mkdir -p ${head_folder}
collapsed_fastas=(${collapsed_folder}*fasta)
if [[ ${#collapsed_fastas[@]} -eq 0 ]]; then
    echo "No collapsed FASTA files were produced after barcode splitting." >&2
    exit 1
fi

for f in "${collapsed_fastas[@]}" ; do
    x="${f##*/}"
    head -n${head} $f > ${head_folder}${x%.*}.fasta
done

# Final processing
mkdir -p ${final_output_folder}
head_fastas=(${head_folder}*fasta)
if [[ ${#head_fastas[@]} -eq 0 ]]; then
    echo "No head FASTA files were available for final assembly." >&2
    exit 1
fi

awk '/>/ { sub(">", ">" FILENAME "_") } 1' "${head_fastas[@]}" > "${final_output_folder}final.fasta"
sed -i 's/>.*\//>/g' "${final_output_folder}final.fasta"
sed -i 's/\.fasta//g' "${final_output_folder}final.fasta"

if [[ ! -s "${final_output_folder}final.fasta" ]]; then
    echo "Final assembled FASTA is missing or empty after genfasta." >&2
    exit 1
fi
