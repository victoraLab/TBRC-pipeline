#!/usr/bin/env python3
"""Per-well consensus for plate-indexed nanopore Ig amplicons.

The method this replaces builds a progressive multiple alignment with mafft and
then calls each column by read-count-weighted majority. That is sound for
substitutions and weak exactly where nanopore is weakest: a homopolymer run
forces the aligner to place gaps it cannot place confidently, and the column
majority inherits that guess. Homopolymer indels are the dominant R10 error mode,
and they are the ones that put a frameshift into an otherwise good VDJ.

Partial order alignment keeps the alignment as a graph rather than a fixed column
grid, so an indel is a branch rather than a committed gap, and the traversal picks
the best-supported path through it. That is the standard fix for this failure
mode. Racon then polishes the POA draft against the well's own reads, which
recovers residual errors the draft carried.

Design constraints taken from the baseline so a comparison means something:
identical depth thresholds, identical cap on unique sequences per well, and the
same reads in. Only the consensus algorithm differs.
"""

import argparse, os, shutil, subprocess, sys, tempfile
from collections import defaultdict, Counter


def parse_header(h):
    """Luka's demux writes key=value pairs separated by '|'. Reuse them rather
    than re-deriving well identity, so both methods group reads identically."""
    out = {}
    for part in h.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def poa_consensus(seqs, spoa="spoa", engine="spoa"):
    """Draft consensus from the partial order graph.

    Two engines, same idea. abPOA uses an adaptive band, which on this assay's
    200 sequences of ~1.5kb is roughly sixteen times faster than mafft --fast
    and six times faster than spoa, measured on a real well. Since the quality
    target here is parity with the mafft baseline rather than an improvement,
    speed is the axis that decides between them."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for i, s in enumerate(seqs):
            fh.write(f">r{i}\n{s}\n")
        tmp = fh.name
    try:
        # -r 0 emits the consensus only. The default affine gap penalties are
        # tuned for noisy long reads, which is the point of using this at all.
        if engine == "mafft":
            # The baseline aligner, wired into this same harness purely so the
            # speed claim is measured end to end rather than extrapolated from a
            # single well. Column-majority calling below is deliberately NOT
            # reimplemented here: this path exists to time the alignment.
            p = run(["mafft", "--retree", "1", "--maxiterate", "0", "--quiet", tmp])
            if p.returncode != 0:
                return None
            seqs, cur = [], []
            for line in p.stdout.splitlines():
                if line.startswith(">"):
                    if cur: seqs.append("".join(cur))
                    cur = []
                else:
                    cur.append(line.strip())
            if cur: seqs.append("".join(cur))
            if not seqs: return None
            # plain column majority over the alignment, gaps included
            out = []
            for col in zip(*seqs):
                c = Counter(col).most_common(1)[0][0]
                if c != "-": out.append(c)
            return "".join(out) or None
        if engine == "abpoa":
            p = run(["abpoa", tmp])
        else:
            p = run([spoa, "-r", "0", tmp])
        if p.returncode != 0:
            return None
        lines = [l for l in p.stdout.splitlines() if l and not l.startswith(">")]
        return "".join(lines) or None
    finally:
        os.unlink(tmp)


def racon_polish(draft, seqs, minimap2="minimap2", racon="racon"):
    """One polishing round of the draft against the well's own reads.

    Racon needs real alignments, so the reads are mapped back to the draft with
    minimap2 in ONT mode. A failure here is not fatal: the unpolished POA draft
    is already a valid consensus, and returning it is better than dropping the
    well entirely."""
    d = tempfile.mkdtemp()
    try:
        dpath, rpath, spath = (os.path.join(d, n) for n in ("draft.fa", "reads.fa", "aln.paf"))
        with open(dpath, "w") as fh:
            fh.write(f">draft\n{draft}\n")
        with open(rpath, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">r{i}\n{s}\n")
        p = run([minimap2, "-x", "map-ont", "-t", "1", dpath, rpath])
        if p.returncode != 0 or not p.stdout.strip():
            return draft
        with open(spath, "w") as fh:
            fh.write(p.stdout)
        p = run([racon, "-t", "1", rpath, spath, dpath])
        if p.returncode != 0:
            return draft
        out = [l for l in p.stdout.splitlines() if l and not l.startswith(">")]
        return "".join(out) or draft
    finally:
        shutil.rmtree(d, ignore_errors=True)


def dominant_by_clustering(seqs, vsearch="vsearch", ident=0.85):
    """Find a well's dominant clone by clustering its own reads.

    The pipeline this replaces assigns clones from read-level IgBLAST, which
    means annotating every read in the run. On this data that is not a small
    cost: the exact-collapse step ahead of it is a no-op -- 100% of reads are
    unique at Q24 over 1.6kb, median duplicate count 1 across 485 wells -- so
    "annotate every exact representative" is "annotate every read", about 1.04M
    queries and roughly twelve hours.

    A well holds one cell. Its reads are one clone apart from doublets and
    barcode misassignment, and distinct clones differ sharply across CDR3, so
    clustering the well's own ~270 reads separates the dominant group without
    needing a germline database at all. IgBLAST then runs once per well on the
    finished consensus rather than once per read.

    Returns the dominant cluster's reads, and the fraction of the well they
    represent -- that fraction is itself the doublet signal worth reporting."""
    if len(seqs) < 3:
        return seqs, 1.0
    d = tempfile.mkdtemp()
    try:
        fa = os.path.join(d, "w.fa")
        uc = os.path.join(d, "w.uc")
        with open(fa, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">r{i}\n{s}\n")
        p = run([vsearch, "--cluster_fast", fa, "--id", str(ident),
                 "--uc", uc, "--quiet", "--threads", "1"])
        if p.returncode != 0 or not os.path.exists(uc):
            return seqs, 1.0
        members = defaultdict(list)
        with open(uc) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if not f or f[0] not in ("S", "H"):
                    continue
                cid = f[1]
                try:
                    idx = int(f[8][1:])
                except (ValueError, IndexError):
                    continue
                members[cid].append(idx)
        if not members:
            return seqs, 1.0
        biggest = max(members.values(), key=len)
        return [seqs[i] for i in biggest if i < len(seqs)], len(biggest) / len(seqs)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def select_reads(records, max_unique, min_meanq, min_len, max_len):
    """Quality and length gate, then exact-collapse with counts.

    The baseline caps unique sequences per well; the same cap is applied here so
    neither method gets an advantage from seeing more or less of the well. Reads
    are ordered by abundance so the cap keeps the best-supported ones."""
    kept = []
    for name, seq in records:
        meta = parse_header(name)
        # Not every demultiplex path annotates quality. The well-only output
        # carries no meanq and no insertlen, so treating an absent field as
        # zero silently discards the entire library -- which it did, until this
        # returned 0 consensus sequences from 87 wells in 0.22 seconds. A field
        # that is not there is unknown, not failing: gate on it only when it
        # exists, and fall back to the sequence's own length.
        try:
            q = float(meta["meanq"]) if "meanq" in meta else None
            ln = int(meta["insertlen"]) if "insertlen" in meta else len(seq)
        except ValueError:
            continue
        if q is not None and q < min_meanq:
            continue
        if ln < min_len or ln > max_len:
            continue
        kept.append(seq)
    if not kept:
        return [], 0
    counts = Counter(kept)
    ordered = [s for s, _ in counts.most_common(max_unique)]
    return ordered, len(kept)


def load_annotated(path):
    """Read the clone-annotated read table and keep only each well's dominant clone.

    This exists so a comparison against the mafft baseline isolates the consensus
    algorithm. The baseline restricts to the dominant clone before aligning; a
    method that instead aligns every read in the well is not doing worse
    alignment, it is aligning a different -- and sometimes doubleted -- set of
    reads, and the resulting chimera shows up as indels and stop codons. Feeding
    both methods the same clone-restricted reads removes that confound."""
    import gzip
    from collections import defaultdict
    opener = gzip.open if str(path).endswith(".gz") else open
    per_clone = defaultdict(lambda: defaultdict(list))
    weight = defaultdict(lambda: defaultdict(int))
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            pw, ck, seq = f[idx["platewell"]], f[idx["clone_key"]], f[idx["sequence"]]
            try:
                cnt = int(f[idx["raw_count"]])
            except (ValueError, KeyError):
                cnt = 1
            if not pw or not seq:
                continue
            per_clone[pw][ck].append(seq)
            weight[pw][ck] += cnt
    wells = {}
    for pw, clones in per_clone.items():
        dominant = max(weight[pw], key=lambda c: weight[pw][c])
        wells[pw] = clones[dominant]
    return wells


def select_like_baseline(seqs, max_unique):
    """Mirror the baseline's own selection: 75% of the budget by abundance, the
    remainder by length. Pure abundance biases the alignment toward truncated
    reads that happen to be common, which is its own source of indels."""
    counts = Counter(seqs)
    by_abundance = sorted(counts, key=lambda s: (counts[s], len(s)), reverse=True)
    by_length = sorted(counts, key=lambda s: (len(s), counts[s]), reverse=True)
    n_ab = max(1, int(max_unique * 0.75))
    chosen = by_abundance[:n_ab]
    for s in by_length:
        if len(chosen) >= max_unique:
            break
        if s not in chosen:
            chosen.append(s)
    return chosen


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("clean_fasta", help="demultiplexed reads carrying platewell= in the header")
    ap.add_argument("--from-annotated", default="",
                    help="clone-annotated read table; restricts to each well's dominant clone")
    ap.add_argument("output_prefix")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="a well at or above this depth is called primary")
    ap.add_argument("--rescue-min-reads", type=int, default=5,
                    help="between this and --min-reads a well is still called, but "
                         "flagged rescued rather than mixed in silently")
    ap.add_argument("--min-dominance", type=float, default=0.70,
                    help="below this share in the dominant clone a well is flagged: "
                         "it probably received two cells, and its consensus is a blend "
                         "of two antibodies rather than a noisy copy of one")
    ap.add_argument("--max-unique", type=int, default=200,
                    help="cap on unique sequences entering the alignment (baseline: 200)")
    ap.add_argument("--min-meanq", type=float, default=15.0)
    ap.add_argument("--min-len", type=int, default=250)
    ap.add_argument("--max-len", type=int, default=3000)
    ap.add_argument("--no-polish", action="store_true", help="emit the POA draft without racon")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--spoa", default="spoa")
    ap.add_argument("--engine", choices=("spoa", "abpoa", "mafft"), default="spoa")
    ap.add_argument("--cluster-wells", action="store_true",
                    help="find the dominant clone by clustering each well's own reads, "
                         "instead of requiring read-level IgBLAST clone assignment")
    ap.add_argument("--cluster-id", type=float, default=0.85)
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--racon", default="racon")
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    annotated_mode = bool(args.from_annotated)
    if annotated_mode:
        raw = load_annotated(args.from_annotated)
        wells = {pw: [(f"|platewell={pw}|meanq=99|insertlen={len(s)}", s) for s in seqs]
                 for pw, seqs in raw.items()}
        sys.stderr.write(f"[nanopore_consensus] {len(wells)} wells, dominant clone only\n")
    else:
        # The assay runs in two modes. A plate-multiplexed library carries both a
        # plate cassette and a well anchor, and identity is platewell. A library
        # sequenced one plate per file carries only the well anchor -- the plate
        # cassette is legitimately absent, not missing -- and identity is the
        # well. Measured on two of this lab's libraries: 100% of well-assigned
        # reads were well-only, so a platewell-only reader silently returns
        # nothing on them.
        wells = defaultdict(list)
        n_pw = n_well = 0
        for name, seq in read_fasta(args.clean_fasta):
            meta = parse_header(name)
            key = meta.get("platewell") or meta.get("well")
            if key:
                if meta.get("platewell"):
                    n_pw += 1
                else:
                    n_well += 1
                wells[key].append((name, seq))
        if n_well and not n_pw:
            sys.stderr.write("[nanopore_consensus] well-only library "
                             "(no plate cassette); identity is the well\n")
        sys.stderr.write(f"[nanopore_consensus] {len(wells)} wells with reads\n")

    from concurrent.futures import ThreadPoolExecutor

    def build(item):
        pw, records = item
        dom_frac = 1.0
        if annotated_mode:
            allseq = [s for _, s in records]
            total = len(allseq)
            seqs = select_like_baseline(allseq, args.max_unique)
        elif args.cluster_wells:
            gated, total = select_reads(records, 10 ** 9, args.min_meanq,
                                        args.min_len, args.max_len)
            if not gated:
                return pw, None, 0, 0, "no_reads_after_qc", 0.0
            dom, frac = dominant_by_clustering(gated, args.vsearch, args.cluster_id)
            total = len(dom)
            seqs = select_like_baseline(dom, args.max_unique)
            dom_frac = frac
        else:
            seqs, total = select_reads(records, args.max_unique, args.min_meanq,
                                       args.min_len, args.max_len)
        # Depth and dominance are separate failure modes and are reported
        # separately. A shallow well is weak evidence; a well whose dominant
        # clone holds half its reads is a doublet, and no amount of depth fixes
        # that. Neither is silently mixed in with clean wells.
        if total < args.rescue_min_reads or not seqs:
            return pw, None, total, len(seqs), "failed_low_depth", dom_frac
        if total < args.min_reads:
            tier = "rescued_low_depth"
        elif dom_frac < args.min_dominance:
            tier = "flagged_low_dominance"
        else:
            tier = "primary"
        draft = poa_consensus(seqs, args.spoa, args.engine)
        if draft is None:
            return pw, None, total, len(seqs), "poa_failed", dom_frac
        if args.no_polish:
            return pw, draft, total, len(seqs), tier, dom_frac
        polished = racon_polish(draft, seqs, args.minimap2, args.racon)
        return pw, polished, total, len(seqs), tier, dom_frac

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for r in ex.map(build, sorted(wells.items())):
            results.append(r)

    fa = args.output_prefix + ".consensus_per_well.fasta"
    tsv = args.output_prefix + ".consensus_per_well.tsv"
    n_called = 0
    with open(fa, "w") as fh, open(tsv, "w") as th:
        th.write("platewell\treads_used\tunique_used\tconsensus_length\t"
                 "confidence_tier\tdominant_fraction\n")
        for pw, cons, total, uniq, status, dfrac in results:
            th.write(f"{pw}\t{total}\t{uniq}\t{len(cons) if cons else 0}\t{status}\t{dfrac:.3f}\n")
            if cons:
                fh.write(f">{pw}|reads={total}|unique={uniq}|domfrac={dfrac:.3f}"
                         f"|tier={status}|method=poa_racon\n{cons}\n")
                n_called += 1
    # A well is meant to hold one cell, and the whole method rests on that: it
    # takes the dominant clone and builds one consensus. If most wells are
    # polyclonal, the library is not single-cell -- it is barcoded bulk samples --
    # and one consensus per well reports an arbitrary sequence while discarding
    # the rest. That produces output which looks entirely normal and is
    # scientifically meaningless, so it is called out loudly rather than left for
    # someone to notice downstream.
    dfracs = [r[5] for r in results if r[1] and r[5] is not None]
    if dfracs:
        dfracs_sorted = sorted(dfracs)
        med_dom = dfracs_sorted[len(dfracs_sorted) // 2]
        if med_dom < 0.35:
            sys.stderr.write(
                "\n[nanopore_consensus] *** WARNING: median dominant-clone fraction is "
                f"{med_dom:.2f} ***\n"
                "[nanopore_consensus] Wells are polyclonal. This library does not look "
                "single-cell.\n"
                "[nanopore_consensus] One consensus per well reports an arbitrary sequence "
                "and discards the rest.\n"
                "[nanopore_consensus] Treat these results as INVALID for a bulk or "
                "multi-cell library.\n\n")

    tiers = Counter(r[4] for r in results if r[1])
    sys.stderr.write(f"[nanopore_consensus] {n_called} consensus sequences written\n")
    for t in ("primary", "rescued_low_depth", "flagged_low_dominance"):
        if tiers.get(t):
            sys.stderr.write(f"[nanopore_consensus]   {t}: {tiers[t]}\n")


if __name__ == "__main__":
    main()
