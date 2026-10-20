#!/usr/bin/env python3
"""Library structure QC for plate-indexed nanopore Ig/TCR amplicons.

The demultiplexer answers "which well is this read from"; it drops everything
else under a handful of counters, and one of them -- plate_assigned_short_insert
-- hides the commonest way a library fails. A read can carry a well, a plate and
both primers and still be dropped as short, because it has 25-49 nt between the
primers: the forward primer joined straight to the reverse primer, a primer
dimer. A library dominated by dimers looks like a run with low yield, but it is
a wet-lab problem, not a sequencing or pipeline one. This script says which it
is.

Every read is put in exactly one class, by its structure:

  full_length      one 5' handle, one reverse-primer cassette downstream, and an
                   insert of at least --min-insert nt (what the pipeline uses)
  short_product    the same, insert between --dimer-max and --min-insert: a
                   template-primed product too short to be a full V(D)J --
                   mispriming inside V (kappa: ~160 nt, primed in IGKV4 FR3)
  primer_dimer     the same, insert under --dimer-max: no template between the
                   primers
  concatemer       two or more 5' handles or reverse cassettes in one read,
                   either strand: amplicons joined in PCR or at ligation
  mixed_strands    one handle and one cassette, but on opposite strands
  no_3p_cassette   a handle and nothing recognisable downstream (truncated read,
                   or a reverse primer not in the barcode set's list)
  no_5p_handle     a cassette and no handle
  unrecognised     neither

The motifs and tolerances are the demultiplexer's own (imported from it), so a
read the demultiplexer would call clean is full_length here. The handle is
matched exactly, as there; the reverse primer within --primer-mm substitutions,
found with pigeonhole seeds so the whole read can be searched, not only its
3' end -- a concatemer's second cassette is in the middle.
"""

import argparse, collections, importlib.util
from multiprocessing import Pool
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "demux_nanopore_wells", Path(__file__).with_name("demux_nanopore_wells.py"))
DMX = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(DMX)
rc, hd = DMX.rc, DMX.hd

CLASSES = ("full_length", "short_product", "primer_dimer", "concatemer", "mixed_strands",
           "no_3p_cassette", "no_5p_handle", "unrecognised")
PRODUCTS = ("full_length", "short_product", "primer_dimer")   # both ends found, one amplicon
BIN = 10          # insert-length histogram bin, nt
MAX_BIN = 1000    # everything longer lands in the last bin

G = {}            # per-worker state, set by init()


def seeds_for(primer, mm):
    """mm+1 non-overlapping pieces: a match with <= mm substitutions contains one exactly."""
    k = len(primer) // (mm + 1)
    return [(i * k, primer[i * k:(i + 1) * k]) for i in range(mm + 1)]


def find_exact(s, motif):
    out, i = [], s.find(motif)
    while i >= 0:
        out.append(i); i = s.find(motif, i + 1)
    return out


def find_cassettes(s):
    """Start positions of any reverse primer (in read orientation) within primer_mm.
    Hits a few nt apart are one cassette matched by two primer variants."""
    hits = {}
    for pr, seeds in G["primers"]:
        for off, seed in seeds:
            for i in find_exact(s, seed):
                st = i - off
                if st < 0 or st + len(pr) > len(s) or st in hits: continue
                if hd(s[st:st + len(pr)], pr) <= G["primer_mm"]: hits[st] = len(pr)
    merged = []
    for st in sorted(hits):
        if merged and st - merged[-1][0] < 10: continue
        merged.append((st, hits[st]))
    return merged


def call_well(t, fp):
    best = sorted((hd(t[fp - len(bc):fp], bc), w) for w, bc in G["wells"].items() if fp >= len(bc))
    if best and best[0][0] <= G["well_mm"] and (len(best) == 1 or best[1][0] > best[0][0]): return best[0][1]
    return None


def call_plate(t, end):
    obs = t[end:end + 10]
    if len(obs) < 10: return None
    best = sorted((hd(obs, b), p) for p, b in G["plates"].items())
    if best[0][0] <= G["plate_mm"] and (len(best) == 1 or best[1][0] > best[0][0]): return best[0][1]
    return None


def classify(seq):
    """-> (class, insert_length or None, platewell or None)"""
    r = rc(seq)
    h5 = [("F", p) for p in find_exact(seq, DMX.PARTIAL_FR1)] + [("R", p) for p in find_exact(r, DMX.PARTIAL_FR1)]
    c3 = [("F",) + c for c in find_cassettes(seq)] + [("R",) + c for c in find_cassettes(r)]
    if len(h5) >= 2 or len(c3) >= 2: return "concatemer", None, None
    if len(h5) == 1 and len(c3) == 1:
        (o5, fp), (o3, pp, plen) = h5[0], c3[0]
        trim5 = fp + len(DMX.FULL_FR1)
        if o5 != o3 or pp < fp: return "mixed_strands", None, None
        t = seq if o5 == "F" else r
        insert = max(0, pp - trim5)
        cls = ("primer_dimer" if insert < G["dimer_max"] else
               "short_product" if insert < G["min_insert"] else "full_length")
        well, plate = call_well(t, fp), call_plate(t, pp + plen)
        return cls, insert, (plate + well if well and plate else None)
    if len(h5) == 1: return "no_3p_cassette", None, None
    if len(c3) == 1: return "no_5p_handle", None, None
    return "unrecognised", None, None


def init(cfg):
    G.update(cfg)
    G["primers"] = [(rc(p), seeds_for(rc(p), cfg["primer_mm"])) for p in cfg["primers_fwd"]]


def work(batch):
    by = collections.defaultdict(collections.Counter)      # (source, class) -> stats
    qh = collections.defaultdict(collections.Counter)      # (source, class) -> int Q histogram
    lh = collections.defaultdict(collections.Counter)      # (source, class) -> read-length histogram (10 nt)
    ins, wells = collections.Counter(), collections.defaultdict(collections.Counter)
    for source, seq, qual in batch:
        cls, insert, pw = classify(seq)
        by[(source, cls)]["reads"] += 1
        qh[(source, cls)][int(DMX.meanq(qual))] += 1
        lh[(source, cls)][len(seq) // 10 * 10] += 1
        if insert is not None: ins[min(insert // BIN * BIN, MAX_BIN)] += 1
        if pw: wells[pw][cls] += 1
    return by, qh, lh, ins, wells


def batches(paths, n=4000):
    b = []
    for source, h, seq, qual in DMX.fastq_records(paths):
        b.append((source, seq, qual))
        if len(b) == n: yield b; b = []
    if b: yield b


def hist_median(h):
    tot = sum(h.values())
    if not tot: return ""
    run = 0
    for k in sorted(h):
        run += h[k]
        if run * 2 >= tot: return k
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", nargs="+")
    ap.add_argument("--well-anchors", required=True)
    ap.add_argument("--plate-barcodes", required=True)
    ap.add_argument("--prefix", required=True, help="writes <prefix>.library_qc.tsv, .library_qc_inserts.tsv, .library_qc_wells.tsv")
    ap.add_argument("--inner-reverse-primer", action="append", default=None)
    ap.add_argument("--primer-mm", type=int, default=2)
    ap.add_argument("--well-barcode-mm", type=int, default=2)
    ap.add_argument("--plate-barcode-mm", type=int, default=1)
    ap.add_argument("--min-insert", type=int, default=250)
    ap.add_argument("--dimer-max", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()

    anchors = DMX.table(a.well_anchors)
    cfg = dict(primers_fwd=[p.upper() for p in (a.inner_reverse_primer or [DMX.INNER_REVERSE_PRIMER])],
               primer_mm=a.primer_mm, well_mm=a.well_barcode_mm, plate_mm=a.plate_barcode_mm,
               min_insert=a.min_insert, dimer_max=a.dimer_max,
               wells={w: s[:-len(DMX.PARTIAL_FR1)] for w, s in anchors.items()},
               plates={p: rc(b) for p, b in DMX.table(a.plate_barcodes).items()})

    by = collections.defaultdict(collections.Counter)
    qh = collections.defaultdict(collections.Counter)
    lh = collections.defaultdict(collections.Counter)
    ins, wells = collections.Counter(), collections.defaultdict(collections.Counter)
    with Pool(max(1, a.jobs), initializer=init, initargs=(cfg,)) as pool:
        for b, q, l, i, w in pool.imap_unordered(work, batches(a.fastq), chunksize=1):
            for k, v in b.items(): by[k].update(v)
            for k, v in q.items(): qh[k].update(v)
            for k, v in l.items(): lh[k].update(v)
            ins.update(i)
            for k, v in w.items(): wells[k].update(v)

    sources = sorted({s for s, _ in by})
    with open(a.prefix + ".library_qc.tsv", "w") as out:
        # fraction_of_products: share among reads with both ends found, i.e. what
        # the PCR made. It is the number to read first: in a noisy run half the
        # reads cannot be classified at all, so a dimer share of all reads
        # understates, several-fold, how much of the product was dimer.
        out.write("scope\tclass\treads\tfraction\tfraction_of_products\tmedian_read_length\tmedian_mean_q\n")
        for scope in ["all"] + sources:
            keys = sources if scope == "all" else [scope]
            tot = sum(by[(s, c)]["reads"] for s in keys for c in CLASSES) or 1
            prod = sum(by[(s, c)]["reads"] for s in keys for c in PRODUCTS) or 1
            for c in CLASSES:
                n = sum(by[(s, c)]["reads"] for s in keys)
                q = collections.Counter(); l = collections.Counter()
                for s in keys: q.update(qh[(s, c)]); l.update(lh[(s, c)])
                fp = f"{n / prod:.4f}" if c in PRODUCTS else ""
                out.write(f"{scope}\t{c}\t{n}\t{n / tot:.4f}\t{fp}\t{hist_median(l)}\t{hist_median(q)}\n")
    with open(a.prefix + ".library_qc_inserts.tsv", "w") as out:
        # The thresholds travel with the histogram so the report draws the ones
        # this run used, not whatever its own defaults are.
        out.write("insert_bin_start\treads\tdimer_max\tmin_insert\n")
        for k in sorted(ins): out.write(f"{k}\t{ins[k]}\t{a.dimer_max}\t{a.min_insert}\n")
    with open(a.prefix + ".library_qc_wells.tsv", "w") as out:
        out.write("platewell\tfull_length\tshort_product\tprimer_dimer\tdimer_fraction\n")
        for pw in sorted(wells):
            f, s, d = (wells[pw][c] for c in ("full_length", "short_product", "primer_dimer"))
            out.write(f"{pw}\t{f}\t{s}\t{d}\t{d / max(f + s + d, 1):.3f}\n")
    tot = sum(v["reads"] for v in by.values()) or 1
    print("library QC: " + ", ".join(f"{c} {sum(by[(s, c)]['reads'] for s in sources) / tot:.1%}" for c in CLASSES))


if __name__ == "__main__":
    main()
