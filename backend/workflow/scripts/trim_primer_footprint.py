#!/usr/bin/env python3
"""Cut each record's 5' primer footprint once its locus is known.

luka_light cut a fixed 23 nt after the 5' handle from every read before
assembly. That is the F_106 kappa primer on a kappa read, but the heavy forward
primer is the handle itself, so on a heavy read the 23 nt were template: every
heavy chain on a mixed luka_light plate lost FR1 nt 24-46 (CMP15_kappa: 424
records). The standard protocol cut nothing and left the kappa primer in the
delivered sequence, where it reads as 19% "mutation" in FR1. The locus decides
how much is primer, and only the annotation knows the locus, so the cut is made
here, after IgBLAST, from the record's own V call. final.script.sh then
annotates the cut FASTA again so the delivered table describes the delivered
sequences.

The input records start where the handle ended. A record whose locus has no
footprint, or that has no V call, is written unchanged.

Usage: trim_primer_footprint.py --fasta F --airr T --footprint IGK=23[,IGL=n] --out O
Prints FOOTPRINT_CUT=<records cut> last, for the caller.
"""
import argparse
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--airr", required=True)
    ap.add_argument("--footprint", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    footprint = {}
    for item in args.footprint.split(","):
        locus, _, n = item.strip().partition("=")
        if locus and n.strip().isdigit() and int(n) > 0:
            footprint[locus.strip().upper()] = int(n)

    locus_of = {}
    with open(args.airr, newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            v_call = (row.get("v_call") or "").strip()
            if v_call:
                locus_of[(row.get("sequence_id") or "").strip()] = v_call.split(",")[0][:3].upper()

    cut = collections.Counter()
    kept = collections.Counter()
    with open(args.out, "w") as out:
        for name, seq in read_fasta(args.fasta):
            locus = locus_of.get(name.split()[0], "")
            n = footprint.get(locus, 0)
            if n and len(seq) > n:
                seq = seq[n:]
                cut[locus] += 1
            else:
                kept[locus or "no V call"] += 1
            out.write(">{}\n{}\n".format(name, seq))

    def fmt(counter):
        return ", ".join("{} {}".format(k, v) for k, v in sorted(counter.items())) or "none"
    print("Primer footprint ({}): cut from {}; left whole: {}.".format(
        args.footprint, fmt(cut), fmt(kept)))
    print("FOOTPRINT_CUT={}".format(sum(cut.values())))


if __name__ == "__main__":
    main()
