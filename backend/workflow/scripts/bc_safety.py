#!/usr/bin/env python3
"""Report how many mismatches a barcode file can tolerate.

fastx_barcode_splitter breaks a tie by taking whichever barcode comes first in
the file -- its comparison is a strict "<", with no ambiguity check and no
warning. Correcting m errors is therefore only unambiguous while every pair of
barcodes differs in at least 2m+1 positions.

Barcode set B violates this: its BC2 has a minimum Hamming distance of 2 over 16
pairs (P01x/P11x and P06x/P12x for every row), so one sequencing error always
resolved to the earlier plate.
"""
import argparse
import itertools
import sys


def read_barcodes(path):
    codes = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                codes.append((parts[0], parts[1].upper()))
    return codes


def min_distance(codes):
    """Minimum pairwise Hamming distance, and the pairs achieving it."""
    best = None
    closest = []
    for (name_a, a), (name_b, b) in itertools.combinations(codes, 2):
        if len(a) != len(b):
            continue
        d = sum(x != y for x, y in zip(a, b))
        if best is None or d < best:
            best, closest = d, [(name_a, name_b)]
        elif d == best:
            closest.append((name_a, name_b))
    return best, closest


def allowed_mismatches(distance):
    if distance is None:
        return 1
    return max(0, (distance - 1) // 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("barcode_file")
    ap.add_argument("--report", help="write a TSV describing the decision")
    ap.add_argument("--label", default="", help="name for this barcode file in the report")
    args = ap.parse_args()

    codes = read_barcodes(args.barcode_file)
    if len(codes) < 2:
        # Nothing to collide with.
        print(1)
        return 0

    distance, closest = min_distance(codes)
    mismatches = allowed_mismatches(distance)

    if args.report:
        with open(args.report, "a") as out:
            if out.tell() == 0:
                out.write("barcode_file\tlabel\tn_barcodes\tmin_hamming\tallowed_mismatches\tcolliding_pairs\n")
            pairs = ";".join(f"{a}~{b}" for a, b in closest) if distance is not None and distance < 3 else ""
            out.write(
                f"{args.barcode_file}\t{args.label}\t{len(codes)}\t{distance}\t{mismatches}\t{pairs}\n"
            )

    if distance is not None and distance < 3:
        print(
            f"{args.label or args.barcode_file}: minimum Hamming distance {distance} over "
            f"{len(closest)} pair(s); demultiplexing with {mismatches} mismatches to avoid "
            f"misassignment.",
            file=sys.stderr,
        )

    print(mismatches)
    return 0


if __name__ == "__main__":
    sys.exit(main())
