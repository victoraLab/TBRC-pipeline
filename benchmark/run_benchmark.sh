#!/bin/bash
# Run every dataset in datasets.tsv through both assembly pipelines and collect
# comparable metrics. Both pipelines are RE-RUN here rather than compared against
# delivered archives, so neither gets the benefit of different tool versions or
# parameters that happened to be in force on the day.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/.." && pwd)"
RAW="${BENCH_RAW:?set BENCH_RAW to the directory holding <dataset>/*.fastq.gz}"
WORK="${BENCH_WORK:?set BENCH_WORK to a scratch directory}"
IMAGE="${BENCH_IMAGE:-tbrc-backend:arm64}"
CORES="${BENCH_CORES:-8}"
RESULTS="${WORK}/results.jsonl"

mkdir -p "${WORK}"
: > "${RESULTS}"

while IFS=$'\t' read -r name primer preprocess chain source description; do
    [ "${name}" = "name" ] && continue
    [ -d "${RAW}/${name}" ] || { echo "skip ${name}: no raw data"; continue; }

    r1=$(ls "${RAW}/${name}"/*_R1_*.fastq.gz 2>/dev/null | head -1)
    [ -n "${r1}" ] || { echo "skip ${name}: no R1"; continue; }
    sample=$(basename "${r1}" | sed -E 's/\.fastq(\.gz)?$//; s/([_.-])R1([_.-]|$)/\1/; s/([_.-])001([_.-]|$)/\1/; s/([_.-])S[0-9]+([_.-]|$)/\1/; s/[_.-]+$//')

    for pipeline in classic_luka tiago_vsearch; do
        dir="${WORK}/${name}/${pipeline}"
        rm -rf "${dir}"; mkdir -p "${dir}/data/runs/bench/${name}" "${dir}/data/barcodes" "${dir}/logs"
        cp -R "${REPO_ROOT}/backend/." "${dir}/"
        cp -R "${REPO_ROOT}/frontend/bc/${primer}" "${dir}/data/barcodes/"
        cp "${RAW}/${name}"/*.fastq.gz "${dir}/data/runs/bench/${name}/"

        python3 - "${dir}/data/runs/bench/${name}/sample.json" "$sample" "$primer" "$preprocess" "$pipeline" "$name" <<'PY'
import json, io, sys
path, sample, primer, preprocess, pipeline, name = sys.argv[1:7]
handles = {"T": "CCAGGGTTTTCCCAGTCACGAC"}
json.dump({
    "user_id": "bench", "user_email": "", "input_method": "Server",
    "sample_name": sample, "server_name": "storage1",
    "server_path": "", "server_results_path": "", "server_ssh_port": "",
    "server_transfer_mode": "pull",
    "primer_set": primer, "trim_primers": "TRUE",
    "sequence_number": 6, "pipeline": pipeline,
    "luka_light_protocol": preprocess,
    "run_igblast": "FALSE", "run_clonality": "FALSE",
    "archive_format": "zip", "archive_scope": "light",
    "igblast_species": "mouse", "igblast_panel": "ig_all",
    "igblast_bin": "igblastn", "igblast_organism": "mouse",
    "igblast_db_v": None, "igblast_db_d": None, "igblast_db_j": None,
    "igblast_data": None, "igblast_auxiliary_data": None,
    "trim_seq": handles.get(primer, "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"),
    "keep_raw": "FALSE",
    "source": f"bench/{name}", "dest": f"results/bench/{name}",
}, io.open(path, "w"), indent=2)
PY

        echo "running ${name} / ${pipeline} ..."
        start=$(date +%s)
        docker run --rm --platform linux/arm64 -v "${dir}":/work -w /work \
            -e DATA_TABLE="data/runs/bench/${name}/sample.json" \
            -e INPUTPATH="data/runs/bench/${name}/" \
            "${IMAGE}" snakemake --cores "${CORES}" > "${dir}/snakemake.log" 2>&1
        rc=$?
        elapsed=$(( $(date +%s) - start ))

        rundir="${dir}/results/bench/${name}/${sample}"
        if [ "${rc}" -eq 0 ] && [ -f "${rundir}/${sample}.final.fasta" ]; then
            python3 "${HERE}/metrics.py" --run-dir "${rundir}" --sample "${sample}" \
                --dataset "${name}" --pipeline "${pipeline}" --runtime "${elapsed}" >> "${RESULTS}"
            echo "  ok in ${elapsed}s"
        else
            python3 -c "
import json,sys
print(json.dumps({'dataset':'${name}','pipeline':'${pipeline}','error':'run failed (rc=${rc})','runtime_sec':${elapsed}}))" >> "${RESULTS}"
            echo "  FAILED after ${elapsed}s (rc=${rc}) -- see ${dir}/snakemake.log"
        fi
    done
done < "${HERE}/datasets.tsv"

echo
echo "results: ${RESULTS}  ($(wc -l < "${RESULTS}" | tr -d ' ') rows)"
