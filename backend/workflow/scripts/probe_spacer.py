#!/usr/bin/env python3
"""Verify the library's spacer length before anything is cut.

The classic pipeline hardcodes `fastx_trimmer -f 3`, i.e. it assumes exactly two
bases precede BC1. If a library is built differently that silently cuts into the
barcode and the run degrades without an error. This reads the head of R1, counts
exact BC1 matches at each candidate offset, and fails if the modal offset is not
the configured one.
"""
import argparse
import collections
import gzip
import sys


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_barcodes(path):
    codes = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                codes.append(parts[1].upper())
    return codes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--expect", type=int, required=True)
    ap.add_argument("--max-offset", type=int, default=6)
    ap.add_argument("--limit", type=int, default=200000)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    codes = set(read_barcodes(args.barcodes))
    if not codes:
        print("No barcodes to probe with; skipping spacer check.", file=sys.stderr)
        return 0
    width = len(next(iter(codes)))

    hits = collections.Counter()
    seen = 0
    with open_maybe_gzip(args.reads) as handle:
        for i, line in enumerate(handle):
            if i % 4 != 1:
                continue
            seq = line.strip().upper()
            seen += 1
            for off in range(args.max_offset + 1):
                if seq[off:off + width] in codes:
                    hits[off] += 1
                    break
            if seen >= args.limit:
                break

    with open(args.report, "w") as out:
        out.write("offset\tmatches\tfraction_of_reads\n")
        for off in range(args.max_offset + 1):
            out.write(f"{off}\t{hits.get(off, 0)}\t{hits.get(off, 0) / seen if seen else 0:.6f}\n")

    if not hits:
        print(
            f"Spacer probe found no barcode match in the first {seen} reads. The barcode set "
            f"probably does not belong to this library.",
            file=sys.stderr,
        )
        return 1

    modal = hits.most_common(1)[0][0]
    share = hits[modal] / seen if seen else 0
    print(f"Spacer probe: modal offset {modal} ({share:.2%} of {seen} reads), configured {args.expect}.")

    if modal != args.expect:
        print(
            f"Spacer mismatch: reads put BC1 at offset {modal}, but the pipeline is configured to "
            f"cut {args.expect} bases. Cutting anyway would remove part of the barcode. Set "
            f"prep.spacer_nt to {modal} in backend/config/pipelines.yaml, or check the barcode set.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
