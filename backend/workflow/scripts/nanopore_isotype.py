#!/usr/bin/env python3
"""Isotype tail and IgA locus calling for per-well nanopore consensus.

Two calls, and they are not the same kind of question.

The isotype TAIL is a presence test: does the consensus contain the IgA or the
IgG2b constant-region motif. A consensus either reached the tail or it did not.

The IgA LOCUS -- wild-type against locus-swapped -- turns on a SINGLE
discriminating base inside the WT motif. One base decided by majority over noisy
reads is exactly the situation where nanopore error becomes a wrong answer, so it
is not read off the consensus. The well's reads are mapped back and the base is
counted directly: the call needs coverage of at least 20 reads AND at least 90%
agreement, and anything short of that is reported as MIXED or
INSUFFICIENT_COVERAGE rather than guessed.

The motifs and both thresholds are taken unchanged from the assay's own
implementation (build_platewell_consensus.py), because they encode which allele
is which.
"""
import argparse, csv, os, shutil, subprocess, sys, tempfile
from collections import Counter, defaultdict

IGA_TAIL = "ACTGTTCGAGGCCCGTTTGGCAGCAAAGAGGTCCCCCAGTAC"
IGG2B_TAIL = "AAGGTAAAGTGGATCTTCTCCTCAGTGGTGGAGCTGAAGCAGAAGATCTCCCCTGACTACAGAAACATGATTGGGCAGGGAGCC"
WT_IGA = "TACTCCTGCATGGTGGGCCACGA"
LS_IGA = "TACTCCTGCATGTTGGGCCACGA"

# Position within WT_IGA that distinguishes wild-type from locus-swapped.
ALLELE_INDEX = 12
MIN_LOCUS_COVERAGE = 20
MIN_ALLELE_FRACTION = 0.90


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            elif line:
                buf.append(line.upper())
    if name is not None:
        yield name, "".join(buf)


def parse_header(h):
    out = {}
    for part in h.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def tail_call(sequence):
    iga, igg2b = IGA_TAIL in sequence, IGG2B_TAIL in sequence
    if iga and igg2b:
        return iga, igg2b, "BOTH_tail"
    if iga:
        return iga, igg2b, "IGA_tail"
    if igg2b:
        return iga, igg2b, "IGG2B_tail"
    return iga, igg2b, "no_tail"


def find_motif(consensus, motif, max_mismatch=2):
    """Best near-match position for the motif, ignoring the allele position.

    The allele base is the thing being measured, so including it in the match
    would bias which allele is found."""
    best = None
    for i in range(max(0, len(consensus) - len(motif) + 1)):
        window = consensus[i:i + len(motif)]
        d = sum(a != b for j, (a, b) in enumerate(zip(window, motif)) if j != ALLELE_INDEX)
        if best is None or d < best[0]:
            best = (d, i)
    if best is None or best[0] > max_mismatch:
        return None
    return best[1]


def allele_pileup(reads, consensus, pos, minimap2="minimap2"):
    """Count bases the well's own reads put at one consensus position."""
    d = tempfile.mkdtemp()
    try:
        ref, qry = os.path.join(d, "c.fa"), os.path.join(d, "r.fa")
        with open(ref, "w") as fh:
            fh.write(f">cons\n{consensus}\n")
        with open(qry, "w") as fh:
            for i, s in enumerate(reads):
                fh.write(f">r{i}\n{s}\n")
        p = subprocess.run([minimap2, "-a", "-x", "map-ont", "-t", "1",
                            "--secondary=no", ref, qry],
                           capture_output=True, text=True)
        if p.returncode != 0:
            return Counter()
        counts = Counter()
        import re as _re
        for line in p.stdout.splitlines():
            if line.startswith("@"):
                continue
            f = line.split("\t")
            if len(f) < 10 or f[2] == "*":
                continue
            try:
                ref_pos = int(f[3]) - 1
            except ValueError:
                continue
            cig, seq = f[5], f[9]
            qi = 0
            for n, op in _re.findall(r"(\d+)([MIDNSHP=X])", cig):
                n = int(n)
                if op in "SH":
                    if op == "S":
                        qi += n
                elif op in "M=X":
                    if ref_pos <= pos < ref_pos + n:
                        base = seq[qi + (pos - ref_pos)]
                        if base in "ACGT":
                            counts[base] += 1
                        break
                    ref_pos += n; qi += n
                elif op == "I":
                    qi += n
                elif op in "DN":
                    if ref_pos <= pos < ref_pos + n:
                        break
                    ref_pos += n
        return counts
    finally:
        shutil.rmtree(d, ignore_errors=True)


def locus_call(counts):
    coverage = sum(counts[b] for b in "ACGT")
    g, t = counts.get("G", 0), counts.get("T", 0)
    gf = g / coverage if coverage else 0.0
    tf = t / coverage if coverage else 0.0
    if coverage < MIN_LOCUS_COVERAGE:
        call = "INSUFFICIENT_COVERAGE"
    elif gf >= MIN_ALLELE_FRACTION:
        call = "WT_IGA"
    elif tf >= MIN_ALLELE_FRACTION:
        call = "LS_IGA"
    else:
        call = "MIXED"
    return coverage, g, t, gf, tf, call


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("clean_fasta")
    ap.add_argument("long_consensus")
    ap.add_argument("output_tsv")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    longs = {}
    for name, seq in read_fasta(args.long_consensus):
        longs[name.split("|")[0]] = seq
    if not longs:
        with open(args.output_tsv, "w") as fh:
            fh.write("platewell\tiga_tail\tigg2b_tail\ttail_call\tiga_locus\t"
                     "locus_coverage\twt_g_reads\tls_t_reads\twt_fraction\tls_fraction\n")
        sys.stderr.write("[isotype] no consensus sequences\n")
        return

    wells = defaultdict(list)
    for name, seq in read_fasta(args.clean_fasta):
        meta = parse_header(name)
        key = meta.get("platewell") or meta.get("well")
        if key in longs:
            wells[key].append(seq)

    from concurrent.futures import ThreadPoolExecutor

    def call(well):
        cons = longs[well]
        iga, igg2b, tail = tail_call(cons)
        pos = find_motif(cons, WT_IGA)
        if pos is None:
            return (well, iga, igg2b, tail, "NOT_REACHED_OR_UNRESOLVED", 0, 0, 0, 0.0, 0.0)
        counts = allele_pileup(wells.get(well, []), cons, pos + ALLELE_INDEX, args.minimap2)
        cov, g, t, gf, tf, locus = locus_call(counts)
        return (well, iga, igg2b, tail, locus, cov, g, t, gf, tf)

    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for r in ex.map(call, sorted(longs)):
            rows.append(r)

    with open(args.output_tsv, "w") as fh:
        fh.write("platewell\tiga_tail\tigg2b_tail\ttail_call\tiga_locus\t"
                 "locus_coverage\twt_g_reads\tls_t_reads\twt_fraction\tls_fraction\n")
        for r in rows:
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(*r))
    tally = Counter(r[3] for r in rows)
    loci = Counter(r[4] for r in rows)
    sys.stderr.write("[isotype] tails: {}\n".format(dict(tally)))
    sys.stderr.write("[isotype] loci : {}\n".format(dict(loci)))


if __name__ == "__main__":
    main()
