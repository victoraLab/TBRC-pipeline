#!/usr/bin/env python3
"""Rank each well's clustered sequences and emit the legacy assembly outputs.

Replaces `fastx_collapser | head -n` with a deterministic, total ordering and an
explicit per-well verdict.

Two behaviours worth stating plainly:

*   Ranking is by cluster size descending, then sequence ascending. The classic
    pipeline left ties to whatever order fastx_collapser happened to produce, so
    two runs of the same data disagreed on ~2% of records.

*   Wells whose dominant sequence has little read support are FLAGGED, not
    dropped. They are the wells the diagnosis identified as primer-dimer
    artefacts, but a thin well is not necessarily a fake one and the decision
    belongs to whoever reads the plate. The verdict is in
    qc/assembly/well_verdicts.tsv; set wells.drop_low_depth to omit them instead.
"""
import argparse
import collections
import glob
import os
import re
import sys

SIZE_RE = re.compile(r"(?:^|;)size=(\d+)")
# Illumina read-through. Should be gone after the prep stage; if it survives to
# here the well is worth looking at.
ADAPTER = "AGATCGGAAGAGC"


def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def size_of(header):
    m = SIZE_RE.search(header)
    return int(m.group(1)) if m else 1


def well_token(path):
    return os.path.basename(path).split(".")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wells-dir", required=True, help="directory of <well>.top.fa files")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--final-fasta", required=True)
    ap.add_argument("--collapsed-stats", required=True)
    ap.add_argument("--verdicts", required=True)
    ap.add_argument("--funnel", required=True)
    ap.add_argument("--topn", type=int, required=True)
    ap.add_argument("--min-dominant-reads", type=int, default=3)
    ap.add_argument("--drop-low-depth", action="store_true")
    ap.add_argument("--multiclonal-ratio", type=float, default=0.2)
    # A floor on cluster size, for pipelines that keep every cluster in a well
    # rather than the top few. Default 1 keeps the existing behaviour exactly.
    ap.add_argument("--min-cluster-reads", type=int, default=1)
    args = ap.parse_args()

    top_files = sorted(glob.glob(os.path.join(args.wells_dir, "*.top.fa")))
    if not top_files:
        print(f"No clustered wells found in {args.wells_dir}.", file=sys.stderr)
        return 1

    for path in (args.final_fasta, args.collapsed_stats, args.verdicts, args.funnel):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    records = []            # (well, rank, count, seq)
    verdicts = []
    per_well_counts = {}
    delivered = 0
    reads_backing = 0
    flagged = collections.Counter()

    for path in top_files:
        well = well_token(path)
        seqs = [(size_of(h), s) for h, s in read_fasta(path) if s]
        if not seqs:
            continue

        # Total, deterministic order: size descending, then sequence ascending.
        seqs.sort(key=lambda t: (-t[0], t[1]))
        kept = [t for t in seqs[: args.topn] if t[0] >= args.min_cluster_reads]
        if not kept:
            # Nothing in the well cleared the floor. Record it rather than
            # letting the well vanish from the verdict table.
            verdicts.append((well, sum(c for c, _ in seqs), len(seqs), seqs[0][0],
                             0.0, "dropped:below_min_cluster_reads"))
            flagged["below_min_cluster_reads"] += 1
            continue

        rank1 = kept[0][0]
        rank2 = kept[1][0] if len(kept) > 1 else 0
        well_reads = sum(c for c, _ in seqs)
        ratio = (rank2 / rank1) if rank1 else 0.0

        if rank1 < args.min_dominant_reads:
            verdict = "low_depth"
        elif args.multiclonal_ratio >= 1.0:
            # A ratio of 1.0 cannot be exceeded, which is how a caller says the
            # distinction does not apply: a MIBS well holds ~100 cells and is
            # supposed to be polyclonal. Calling it "monoclonal" would be a
            # confident wrong answer written into the delivered table.
            verdict = "multiwell"
        elif ratio > args.multiclonal_ratio:
            verdict = "multiclonal"
        else:
            verdict = "monoclonal"

        if any(ADAPTER in s for _, s in kept):
            verdict += "+adapter"

        flagged[verdict] += 1
        low = verdict.startswith("low_depth")

        if low and args.drop_low_depth:
            verdicts.append((well, well_reads, len(seqs), rank1, ratio, "dropped:" + verdict))
            continue

        verdicts.append((well, well_reads, len(seqs), rank1, ratio, verdict))
        per_well_counts[well] = [c for c, _ in kept]
        for rank, (count, seq) in enumerate(kept, start=1):
            records.append((well, rank, count, seq))
            delivered += 1
            reads_backing += count

    if not records:
        print("Every well was dropped; nothing to deliver.", file=sys.stderr)
        return 1

    # Legacy header contract: ><sample><BC1><BC2>_<rank>-<count>
    with open(args.final_fasta, "w") as out:
        for well, rank, count, seq in records:
            out.write(f">{args.sample}{well}_{rank}-{count}\n{seq}\n")

    # plate_stats.R expects one column per well, counts descending down the rows.
    with open(args.collapsed_stats, "w") as out:
        wells = sorted(per_well_counts)
        out.write("\t".join(f"{args.sample}{w}" for w in wells) + "\n")
        depth = max(len(v) for v in per_well_counts.values())
        for i in range(depth):
            row = [str(per_well_counts[w][i]) if i < len(per_well_counts[w]) else "" for w in wells]
            out.write("\t".join(row) + "\n")

    with open(args.verdicts, "w") as out:
        out.write("well\treads_in_well\tn_clusters\trank1_count\trank2_ratio\tverdict\n")
        for well, reads, nclust, rank1, ratio, verdict in sorted(verdicts):
            out.write(f"{well}\t{reads}\t{nclust}\t{rank1}\t{ratio:.4f}\t{verdict}\n")

    with open(args.funnel, "a") as out:
        if out.tell() == 0:
            out.write("stage\treads\twells\n")
        out.write(f"delivered_records\t{reads_backing}\t{len(per_well_counts)}\n")

    summary = ", ".join(f"{v}={n}" for v, n in sorted(flagged.items()))
    print(f"Delivered {delivered} records from {len(per_well_counts)} wells "
          f"({reads_backing} reads backing them). Verdicts: {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
