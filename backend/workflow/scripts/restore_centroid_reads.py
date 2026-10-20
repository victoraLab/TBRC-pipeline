#!/usr/bin/env python3
"""Put a whole read back in place of each cluster's primer-free centroid.

A well is clustered on copies of its reads with the 5' primer footprint removed,
because the degenerate kappa primer split identical clones into separate
clusters: on CMP15, 153 of 186 "multiclonal" kappa wells were one clone, and 752
of the 770 differences between their clusters sat inside the 23-nt primer
(docs/cmp15-mixed-plate-study-2026-09-17.md, finding 4).

What is delivered is still a real read, untrimmed: among the reads whose
primer-free copy IS the centroid, the most common whole sequence (ties to the
alphabetically first, so reruns agree). Those reads share the centroid's
template exactly, so the delivered template is the one clustering chose; the
vote only settles the masked stretch, where one read's sequencing error would
otherwise be delivered as a mutation on a heavy chain (FR1 nt 24-46 are
template there). Voting across the whole cluster instead let primer
degeneracy decide: the true template's reads are split across primer
variants, so an error copy on a single variant could outvote it (CMP15 kappa,
C07P14: 2 substitutions in a 912-read well). Membership comes from the
dereplication's --uc table; the label, and so the cluster size, is kept.

Usage: restore_centroid_reads.py <top.fa> <well.fastq> <uniques.uc> <out.fa>
"""
import collections
import sys


def strip_size(label):
    return label.split(";size=")[0]


def read_fastq(path):
    reads = {}
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            # vsearch labels a sequence with its header up to the first space.
            reads[header[1:].split()[0]] = seq
    return reads


def parent_map(uc_path):
    """query -> the seed it joined (a seed maps to itself)."""
    parent = {}
    with open(uc_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] not in ("S", "H"):
                continue
            query = strip_size(fields[8])
            parent[query] = query if fields[0] == "S" else strip_size(fields[9])
    return parent


def read_fasta(path):
    label, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if label is not None:
                    yield label, "".join(chunks)
                label, chunks = line[1:], []
            else:
                chunks.append(line.strip())
    if label is not None:
        yield label, "".join(chunks)


def main(top, fastq, uniques_uc, out):
    reads = read_fastq(fastq)
    members = collections.defaultdict(list)
    for read, unique in parent_map(uniques_uc).items():
        members[unique].append(read)

    missing = 0
    with open(out, "w") as dest:
        for label, clustered in read_fasta(top):
            counts = collections.Counter(reads[r] for r in members.get(strip_size(label), ()) if r in reads)
            if counts:
                best = max(counts.values())
                seq = min(s for s, n in counts.items() if n == best)
            else:
                missing += 1
                seq = clustered
            dest.write(">{}\n{}\n".format(label, seq))
    if missing:
        print("restore_centroid_reads: {} centroid(s) in {} had no reads to restore from; "
              "kept the clustered sequence".format(missing, top), file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:5])
