#!/usr/bin/env python3
"""Compute comparable metrics from a finished TBRC run directory.

Everything here is measured from the delivered FASTA and the run's own QC
output, so the same numbers can be computed for either assembly pipeline
without either one being given special treatment.
"""
import argparse
import collections
import json
import os
import re
import statistics
import sys

ADAPTER = "AGATCGGAAGAGC"          # Illumina read-through
DIMER_TAIL = "CACTGACACTCCCACGAGCTC"  # reverse primer, seen early = zero insert
HEADER = re.compile(r"^>(?P<sample>.*?)(?P<well>[A-Za-z0-9]+)_(?P<rank>\d+)-(?P<count>\d+)$")


def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line, []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse(header, sample):
    body = header[1:]
    if sample and body.startswith(sample):
        body = body[len(sample):]
    m = re.match(r"^(?P<well>.+)_(?P<rank>\d+)-(?P<count>\d+)$", body)
    if not m:
        return None
    return m.group("well"), int(m.group("rank")), int(m.group("count"))


def hamming_le(a, b, k):
    if abs(len(a) - len(b)) > k:
        return False
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--runtime", type=float, default=float("nan"))
    args = ap.parse_args()

    fasta = os.path.join(args.run_dir, f"{args.sample}.final.fasta")
    if not os.path.exists(fasta):
        print(json.dumps({"dataset": args.dataset, "pipeline": args.pipeline, "error": "no final.fasta"}))
        return 1

    wells = collections.defaultdict(list)
    seqs = {}
    for header, seq in read_fasta(fasta):
        p = parse(header, args.sample)
        if not p:
            continue
        well, rank, count = p
        wells[well].append((rank, count, seq))
        seqs[header] = seq

    if not wells:
        print(json.dumps({"dataset": args.dataset, "pipeline": args.pipeline, "error": "no parsable records"}))
        return 1

    records = sum(len(v) for v in wells.values())
    adapter = sum(1 for s in seqs.values() if ADAPTER in s)
    dimer = sum(1 for s in seqs.values() if DIMER_TAIL in s[:130])

    top1 = []
    ones = 0
    near = pairs = 0
    lengths = []
    reads_backing = 0
    for v in wells.values():
        v.sort()
        c1 = max(c for _, c, _ in v)
        top1.append(c1)
        if c1 == 1:
            ones += 1
        reads_backing += sum(c for _, c, _ in v)
        lengths.extend(len(s) for _, _, s in v)
        if len(v) > 1:
            pairs += 1
            if hamming_le(v[0][2], v[1][2], 2):
                near += 1

    lengths.sort()

    # the pipeline's own funnel, when it wrote one
    funnel = {}
    fpath = os.path.join(args.run_dir, "qc", "demux", "funnel.tsv")
    if os.path.exists(fpath):
        for line in open(fpath):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "stage":
                try:
                    funnel[parts[0]] = int(parts[1])
                except ValueError:
                    pass

    out = {
        "dataset": args.dataset,
        "pipeline": args.pipeline,
        "runtime_sec": round(args.runtime, 1),
        "records": records,
        "wells": len(wells),
        "adapter_pct": round(100 * adapter / records, 2),
        "dimer_pct": round(100 * dimer / records, 2),
        "singleton_well_pct": round(100 * ones / len(wells), 2),
        "rank2_near_rank1_pct": round(100 * near / pairs, 2) if pairs else 0.0,
        "support_median": statistics.median(top1),
        "support_mean": round(statistics.mean(top1), 1),
        "reads_backing_records": reads_backing,
        "length_median": lengths[len(lengths) // 2] if lengths else 0,
        "length_p90": lengths[int(0.9 * len(lengths))] if lengths else 0,
        "funnel": funnel,
    }
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
