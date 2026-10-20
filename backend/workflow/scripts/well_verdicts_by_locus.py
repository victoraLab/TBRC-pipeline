#!/usr/bin/env python3
"""Call each well's verdict per locus, from the annotated delivery.

rank_wells.py judges a well on all its clusters together, so a well holding one
cell's heavy chain and its kappa chain -- the normal case on a mixed plate -- is
called multiclonal: 194 of CMP15_all's 331 wells, against 53 on the
light-filtered run (handoff 2026-09-17, section 4). Only the annotation knows
which locus a record is, so this runs after it: the delivered records of each
well are grouped by locus (the prefix of the first V call) and ranked within it,
with the same thresholds rank_wells.py uses.

It sees what was delivered, not every cluster: a well keeps the top
sequence_number/2 records across both loci, so a second cluster of one locus
can fall below the cut when the other locus is present. A "monoclonal" here is
therefore an upper bound; a "multiclonal" is not.

Usage: well_verdicts_by_locus.py --fasta imgt.ready.fasta --airr igblast.tsv --sample S --out TSV
         [--min-dominant-reads 3] [--multiclonal-ratio 0.2]
"""
import argparse
import collections
import csv
import re
import sys

csv.field_size_limit(sys.maxsize)
HEADER = re.compile(r"^(?P<well>.+)_(?P<rank>\d+)-(?P<count>\d+)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--airr", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-dominant-reads", type=int, default=3)
    ap.add_argument("--multiclonal-ratio", type=float, default=0.2)
    args = ap.parse_args()

    locus_of = {}
    with open(args.airr, newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            v_call = (row.get("v_call") or "").strip()
            locus_of[(row.get("sequence_id") or "").strip()] = (
                v_call.split(",")[0][:3].upper() if v_call else "none")

    counts = collections.defaultdict(list)    # (well, locus) -> read counts
    with open(args.fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].split()[0]
            token = name[len(args.sample):] if name.startswith(args.sample) else name
            m = HEADER.match(token)
            if not m:
                continue
            counts[(m.group("well"), locus_of.get(name, "none"))].append(int(m.group("count")))

    with open(args.out, "w") as out:
        out.write("well\tlocus\trecords\trank1_count\trank2_ratio\tverdict\n")
        for (well, locus), sizes in sorted(counts.items()):
            sizes.sort(reverse=True)
            rank1 = sizes[0]
            ratio = sizes[1] / rank1 if len(sizes) > 1 and rank1 else 0.0
            if locus == "none":
                verdict = "unannotated"
            elif rank1 < args.min_dominant_reads:
                verdict = "low_depth"
            elif ratio > args.multiclonal_ratio:
                verdict = "multiclonal"
            else:
                verdict = "monoclonal"
            out.write("{}\t{}\t{}\t{}\t{:.4f}\t{}\n".format(well, locus, len(sizes), rank1, ratio, verdict))

    by_locus = collections.Counter(locus for (_, locus) in counts)
    print("Per-locus well verdicts: {} (well, locus) pairs; loci {}.".format(
        len(counts), dict(sorted(by_locus.items()))))


if __name__ == "__main__":
    main()
