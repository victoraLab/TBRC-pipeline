#!/usr/bin/env python3
"""Build the second per-well consensus: VDJ only, quality-filtered.

This assay needs two sequences per well, for two different questions:

  long  the whole amplicon including the constant-region tail, used for
        tail-swap and IgA locus-swap analysis. Built by nanopore_consensus.py.
  VDJ   V through J only, called more strictly, used for IgBLAST functionality
        and V-gene substitution analysis.

They are not the same sequence cropped two ways, and the reason matters. In
mutation analysis a miscalled base becomes a reported somatic hypermutation, so
the VDJ consensus is called with a higher majority threshold (0.70 against 0.60)
and drops reads below Q12 -- deliberately trading a few Ns for not inventing
mutations.

The reference implementation crops each READ to VDJ using that read's own
IgBLAST coordinates, which requires annotating every read: about 1.04M queries
and roughly 11.6 hours on a full run. This gets the coordinates from the long
consensus instead -- one IgBLAST query per well -- and then maps the well's reads
back onto it to recover each read's VDJ span. minimap2 does that in milliseconds
per well.

One genuine consequence of cropping to VDJ, worth knowing: the exact-collapse
that is a no-op on full-length reads becomes meaningful here. At Q24 an
error-free 1.6kb read is a 0.2% event, but an error-free 400 nt VDJ is roughly
20%, so identical VDJ sequences do recur and the read counts become real
weights rather than all being 1.
"""

import argparse, os, shutil, subprocess, sys, tempfile
from collections import Counter, defaultdict


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


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def vdj_bounds(airr_path):
    """V start and J end per well, from the long consensus's own annotation."""
    import csv
    csv.field_size_limit(10 ** 9)
    out = {}
    if not os.path.exists(airr_path):
        return out
    with open(airr_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            well = r.get("sequence_id", "").split("|")[0]
            try:
                start = int(float(r.get("v_sequence_start") or 0))
                end = int(float(r.get("j_sequence_end") or 0))
            except (TypeError, ValueError):
                continue
            if well and start > 0 and end > start:
                out[well] = (start, end)
    return out


def crop_reads_to_consensus(reads, consensus, start, end, minimap2="minimap2"):
    """Map a well's reads onto its own consensus and take each read's VDJ span.

    Using the consensus as the reference means the coordinates come from one
    annotation rather than one per read, and a read only has to overlap the
    region rather than span the whole amplicon."""
    d = tempfile.mkdtemp()
    try:
        ref, qry = os.path.join(d, "c.fa"), os.path.join(d, "r.fa")
        with open(ref, "w") as fh:
            fh.write(f">cons\n{consensus}\n")
        with open(qry, "w") as fh:
            for i, s in enumerate(reads):
                fh.write(f">r{i}\n{s}\n")
        p = run([minimap2, "-x", "map-ont", "-t", "1", "--secondary=no", ref, qry])
        if p.returncode != 0 or not p.stdout.strip():
            return []
        cropped = []
        for line in p.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            try:
                qi = int(f[0][1:]); qs, qe = int(f[2]), int(f[3])
                ts, te = int(f[7]), int(f[8])
            except (ValueError, IndexError):
                continue
            # Keep only reads whose alignment actually covers the VDJ region.
            if ts > start or te < end:
                continue
            span = te - ts
            if span <= 0:
                continue
            # Project consensus coordinates onto the read, proportionally. This
            # is approximate, and deliberately so: the consensus is then called
            # by majority over many reads, which absorbs a few bases of slop at
            # the edges far better than it would absorb a wrong read.
            rs = qs + int((start - ts) / span * (qe - qs))
            re_ = qs + int((end - ts) / span * (qe - qs))
            seq = reads[qi][max(0, rs):min(len(reads[qi]), re_)]
            if len(seq) > 50:
                cropped.append(seq)
        return cropped
    finally:
        shutil.rmtree(d, ignore_errors=True)


def call_consensus(seqs, weights, min_depth, min_fraction, spoa="abpoa"):
    """Majority call over an alignment of the well's VDJ reads.

    Uses abPOA to align and then calls columns explicitly rather than taking
    abPOA's own consensus, because this step needs the stricter majority rule
    and the ability to emit N where support is thin. A base that is not clearly
    supported must become N: in mutation analysis an N is a gap in the evidence,
    while a wrong base is a somatic hypermutation that never happened."""
    if not seqs:
        return "", "no_reads"
    d = tempfile.mkdtemp()
    try:
        fa = os.path.join(d, "v.fa")
        with open(fa, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">r{i}\n{s}\n")
        # -r1 emits the row-column MSA rather than the consensus.
        p = run([spoa, "-r", "1", fa])
        if p.returncode != 0 or not p.stdout.strip():
            return "", "align_failed"
        aln, cur, names = [], [], []
        for line in p.stdout.splitlines():
            if line.startswith(">"):
                if cur:
                    aln.append("".join(cur))
                cur = []
                names.append(line[1:].strip())
            else:
                cur.append(line.strip())
        if cur:
            aln.append("".join(cur))
        aln = [a for a in aln if a]
        if not aln:
            return "", "align_empty"
        width = len(aln[0])
        w = [weights[i] if i < len(weights) else 1 for i in range(len(aln))]
        # Only count a sequence at columns it actually spans, so a read that
        # stops early does not vote "gap" across the remainder.
        spans = []
        for a in aln:
            occ = [i for i, b in enumerate(a) if b != "-"]
            spans.append((occ[0], occ[-1]) if occ else (0, -1))
        out, ns = [], 0
        for col in range(width):
            bases, gaps = Counter(), 0
            for n, a in enumerate(aln):
                if spans[n][0] <= col <= spans[n][1]:
                    b = a[col]
                    if b == "-":
                        gaps += w[n]
                    elif b in "ACGT":
                        bases[b] += w[n]
            depth = sum(bases.values())
            if depth < min_depth or gaps > depth:
                continue
            base, support = bases.most_common(1)[0]
            if support / depth >= min_fraction:
                out.append(base)
            else:
                out.append("N"); ns += 1
        seq = "".join(out)
        if not seq:
            return "", "insufficient_depth"
        return seq, ("pass" if ns == 0 else "contains_low_confidence_N")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("clean_fasta", help="demultiplexed reads")
    ap.add_argument("long_consensus", help="per-well long consensus FASTA")
    ap.add_argument("long_airr", help="AIRR annotation of the long consensus")
    ap.add_argument("output_prefix")
    ap.add_argument("--min-mean-q", type=float, default=12.0,
                    help="stricter than the long consensus: a low-quality read "
                         "contributes a false mutation, not just noise")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--max-unique", type=int, default=300)
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-base-fraction", type=float, default=0.70)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--abpoa", default="abpoa")
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    bounds = vdj_bounds(args.long_airr)
    sys.stderr.write(f"[vdj] VDJ coordinates for {len(bounds)} wells\n")
    if not bounds:
        sys.stderr.write("[vdj] no annotated wells; nothing to do\n")
        open(args.output_prefix + ".vdj_per_well.fasta", "w").close()
        open(args.output_prefix + ".vdj_per_well.tsv", "w").write(
            "platewell\treads_used\tunique_used\tvdj_length\tn_count\tstatus\n")
        return

    longs = {}
    for name, seq in read_fasta(args.long_consensus):
        longs[name.split("|")[0]] = seq

    wells = defaultdict(list)
    for name, seq in read_fasta(args.clean_fasta):
        meta = parse_header(name)
        key = meta.get("platewell") or meta.get("well")
        if not key or key not in bounds:
            continue
        q = meta.get("meanq")
        if q is not None:
            try:
                if float(q) < args.min_mean_q:
                    continue
            except ValueError:
                pass
        wells[key].append(seq)

    from concurrent.futures import ThreadPoolExecutor

    def build(item):
        well, reads = item
        if len(reads) < args.min_reads or well not in longs:
            return well, "", 0, 0, 0, "below_min_reads"
        start, end = bounds[well]
        cropped = crop_reads_to_consensus(reads, longs[well], start, end, args.minimap2)
        if len(cropped) < args.min_reads:
            return well, "", len(reads), 0, 0, "too_few_spanning_reads"
        # Cropping to VDJ makes exact collapse meaningful, so counts become real
        # weights here rather than being uniformly 1 as they are full-length.
        counts = Counter(cropped)
        ranked = sorted(counts, key=lambda s: (counts[s], len(s)), reverse=True)[:args.max_unique]
        weights = [counts[s] for s in ranked]
        seq, status = call_consensus(ranked, weights, args.min_depth,
                                     args.min_base_fraction, args.abpoa)
        return well, seq, len(cropped), len(ranked), seq.count("N"), status

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for r in ex.map(build, sorted(wells.items())):
            results.append(r)

    fa = args.output_prefix + ".vdj_per_well.fasta"
    tsv = args.output_prefix + ".vdj_per_well.tsv"
    n = 0
    with open(fa, "w") as fh, open(tsv, "w") as th:
        th.write("platewell\treads_used\tunique_used\tvdj_length\tn_count\tstatus\n")
        for well, seq, used, uniq, ncount, status in results:
            th.write(f"{well}\t{used}\t{uniq}\t{len(seq)}\t{ncount}\t{status}\n")
            if seq:
                fh.write(f">{well}|reads={used}|unique={uniq}|N={ncount}|region=VDJ\n{seq}\n")
                n += 1
    sys.stderr.write(f"[vdj] {n} VDJ consensus sequences written\n")


if __name__ == "__main__":
    main()
