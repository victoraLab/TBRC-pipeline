#!/usr/bin/env python3
# Plate + well demultiplexing for nanopore Ig/TCR amplicons.
#
# VENDORED, essentially unchanged, from the Victora lab's
# Nanopore_single_cell_IGH_TS_LS pipeline
# (bin/demultiplex_nanopore_plate_longreads.py), written by Luka. It is kept
# here so a TBRC run does not depend on a script living in someone's Dropbox,
# and so the version a result was produced with is recorded alongside it.
#
# It is vendored rather than rewritten deliberately. The barcode tolerances,
# the terminal search window, the orientation handling and the both-end trimming
# encode real knowledge about this assay, and it demultiplexes correctly: on the
# libraries tested it assigned 50-72% of raw reads to a well, and where a plate
# cassette was absent it correctly reported the reads as well-only rather than
# discarding them. Replacing working, assay-specific code with a reimplementation
# would risk a regression for no benefit.
#
# Upstream is the authority. If it changes, re-vendor rather than diverge.
#
# ONE local change, and it is a parameter, not a rewrite: the inner reverse
# primer that the plate cassette sits behind was a constant (the IgA primer the
# upstream libraries used), so a library amplified with any other constant-region
# primer had every plate cassette missed -- on IgG heavy and kappa light MinION
# libraries no well-assigned read got a plate, although the plate barcode and the
# adapter are right there behind the IgG and kappa primers.
# --inner-reverse-primer (repeatable) defaults to the upstream constant, so a run
# that does not pass it demultiplexes exactly as before. Carry this forward when
# re-vendoring.
"""Orient, Well+Plate demultiplex, and trim both ends of long VDJ reads."""

import argparse, gzip, math, re
from collections import Counter, defaultdict
from pathlib import Path

PARTIAL_FR1 = "GGGAATTCGAGGTGCAGCTGCAGGA"
FULL_FR1 = "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"
A1_CORE = "TACACGACGCTCTTCCGATCT"
INNER_REVERSE_PRIMER = "CCACCATCTCCATGACTGCCA"
ADAPTER = "CTGCTGAACCGCTCTTCCGATCT"


def rc(s): return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
def hd(a, b): return sum(x != y for x, y in zip(a, b))
def meanq(q): return sum(ord(x)-33 for x in q) / len(q) if q else 0.0


def contains_hd(s, motif, mm):
    return any(hd(s[i:i+len(motif)], motif) <= mm for i in range(max(0, len(s)-len(motif)+1)))


def motif_positions(s, motif, mm, start=0):
    for i in range(start, len(s)-len(motif)+1):
        d = hd(s[i:i+len(motif)], motif)
        if d <= mm: yield i, d


def openfq(p): return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)
def sid(h): return re.sub(r"[^A-Za-z0-9_.:-]", "_", h[1:].strip().split()[0])
def source_id(path):
    name = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(suffix): return name[:-len(suffix)]
    return name


def fastq_records(paths):
    seen = Counter()
    for path in paths:
        base = source_id(path); seen[base] += 1
        source = base if seen[base] == 1 else f"{base}_{seen[base]}"
        with openfq(path) as fin:
            while True:
                h = fin.readline()
                if not h: break
                seq = fin.readline().strip().upper(); plus = fin.readline(); qual = fin.readline().strip()
                if not plus or len(seq) != len(qual): raise RuntimeError(f"Malformed FASTQ near {h.strip()} in {path}")
                yield source, h, seq, qual


def table(path):
    out = {}
    with open(path) as h:
        for line in h:
            if line.strip():
                k, v = line.split()[:2]; out[k] = v.upper()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq"); ap.add_argument("well_anchors"); ap.add_argument("plate_barcodes"); ap.add_argument("prefix")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--extra-fastq", action="append", default=[], help="Additional FASTQ or FASTQ.gz from the same library")
    ap.add_argument("--well-barcode-mm", type=int, default=2)
    ap.add_argument("--plate-barcode-mm", type=int, default=1)
    ap.add_argument("--primer-mm", type=int, default=2)
    ap.add_argument("--terminal-window", type=int, default=300)
    ap.add_argument("--min-insert-length", type=int, default=250)
    ap.add_argument("--inner-reverse-primer", action="append", default=None,
                    help="Constant-region primer the plate cassette follows, 5'->3' as synthesised; repeatable. "
                         "Default: the upstream IgA primer " + INNER_REVERSE_PRIMER)
    args = ap.parse_args()

    anchors, plates_forward = table(args.well_anchors), table(args.plate_barcodes)
    if len(anchors) != 96 or len(plates_forward) != 16: raise SystemExit("Expected 96 well anchors and 16 plate barcodes")
    well_barcodes = {w:a[:-len(PARTIAL_FR1)] for w,a in anchors.items()}
    plates_read = {p:rc(b) for p,b in plates_forward.items()}
    primers_read = [rc(p.upper()) for p in (args.inner_reverse_primer or [INNER_REVERSE_PRIMER])]
    adapter_read = rc(ADAPTER)

    prefix = Path(args.prefix)
    clean_path = Path(str(prefix)+".platewell_clean.fasta")
    partial_path = Path(str(prefix)+".well_only_partial.fasta")
    qc_path = Path(str(prefix)+".platewell_read_qc.tsv")
    summary_path = Path(str(prefix)+".platewell_summary.tsv")
    counts, by_source, by_cell, by_well = Counter(), defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
    inputs = [args.fastq] + args.extra_fastq
    for path in inputs:
        if not Path(path).is_file(): raise SystemExit(f"Missing FASTQ: {path}")

    def bump(metric, source, amount=1):
        counts[metric] += amount
        by_source[source][metric] += amount

    with open(clean_path,"w") as clean, open(partial_path,"w") as partial, open(qc_path,"w") as qc:
        qc.write("source_run\tread_id\tplate\twell\tplatewell\torientation\twell_barcode_mm\tplate_barcode_mm\tprimer_mm\tread_length\tinsert_length\tmean_q\ta1\tfull_fr1\tadapter_verified\tstatus\n")
        for source,h,seq,qual in fastq_records(inputs):
            bump("raw_reads",source); rid=sid(h)

            wm=[]
            for ori,s,q in (("F",seq,qual),("R",rc(seq),qual[::-1])):
                pos=s.find(PARTIAL_FR1)
                while pos>=0:
                    for well,bc in well_barcodes.items():
                        if pos>=len(bc):
                            d=hd(s[pos-len(bc):pos],bc)
                            if d<=args.well_barcode_mm: wm.append((d,well,ori,pos,s,q))
                    pos=s.find(PARTIAL_FR1,pos+1)
            if not wm:
                bump("no_well_anchor",source); continue
            dmin=min(x[0] for x in wm); best=[x for x in wm if x[0]==dmin]
            if len({x[1] for x in best})!=1:
                bump("ambiguous_well",source); continue
            wmm,well,ori,fp,s,q=best[0]
            bump("well_assigned",source); by_well[well]["well_assigned"]+=1
            full_obs=s[fp:fp+len(FULL_FR1)]
            full_fr1=len(full_obs)==len(FULL_FR1) and hd(full_obs,FULL_FR1)<=2
            trim5=fp+len(FULL_FR1)
            if trim5>len(s): bump("incomplete_fr1",source); continue
            a1=contains_hd(s[max(0,fp-180):fp],A1_CORE,2)

            pm=[]
            search_start=max(trim5,len(s)-args.terminal_window)
            for primer_read in primers_read:
                for pp,pmm in motif_positions(s,primer_read,args.primer_mm,search_start):
                    bstart=pp+len(primer_read); bend=bstart+10
                    if bend>len(s): continue
                    observed=s[bstart:bend]
                    calls=sorted((hd(observed,b),p) for p,b in plates_read.items())
                    if calls[0][0]<=args.plate_barcode_mm:
                        tied=[p for d,p in calls if d==calls[0][0]]
                        if len(tied)==1: pm.append((pmm,calls[0][0],tied[0],pp,bend,s))
            if not pm:
                bump("well_only_missing_plate_cassette",source); by_well[well]["partial"]+=1
                insert=s[trim5:]
                partial.write(f">{args.sample}|source_run={source}|well={well}|read={rid}|ori={ori}|wmm={wmm}|status=well_only|n=1\n{insert}\n")
                continue
            pm.sort(key=lambda x:(x[0]+x[1],-x[3])); score=pm[0][0]+pm[0][1]
            tied={(x[2],x[3]) for x in pm if x[0]+x[1]==score}
            if len(tied)!=1:
                bump("ambiguous_plate",source); by_well[well]["ambiguous_plate"]+=1; continue
            pmm,bmm,plate,pp,bend,_=pm[0]; cell=plate+well
            insert=s[trim5:pp]
            if len(insert)<args.min_insert_length:
                bump("plate_assigned_short_insert",source); by_cell[cell]["short"]+=1; continue
            after_barcode=s[bend:]
            adapter_verified=contains_hd(after_barcode,adapter_read,3)
            bump("platewell_clean_reads",source); bump(f"plate_{plate}",source)
            by_cell[cell]["clean"]+=1; by_cell[cell]["a1"]+=int(a1); by_cell[cell]["adapter_verified"]+=int(adapter_verified)
            header=(f"{args.sample}|source_run={source}|plate={plate}|well={well}|platewell={cell}|read={rid}|ori={ori}|"
                    f"wmm={wmm}|pmm={bmm}|primermm={pmm}|len={len(seq)}|insertlen={len(insert)}|"
                    f"meanq={meanq(q):.2f}|a1={int(a1)}|fullfr1={int(full_fr1)}|adapter={int(adapter_verified)}|strict=1|n=1")
            clean.write(f">{header}\n{insert}\n")
            qc.write(f"{source}\t{rid}\t{plate}\t{well}\t{cell}\t{ori}\t{wmm}\t{bmm}\t{pmm}\t{len(seq)}\t{len(insert)}\t{meanq(q):.2f}\t{int(a1)}\t{int(full_fr1)}\t{int(adapter_verified)}\tplatewell_clean\n")

    with open(summary_path,"w") as out:
        out.write("section\tid\tmetric\tvalue\n")
        for k,v in sorted(counts.items()): out.write(f"global\tALL\t{k}\t{v}\n")
        for source in sorted(by_source):
            for k,v in sorted(by_source[source].items()): out.write(f"source\t{source}\t{k}\t{v}\n")
        for w in sorted(anchors):
            for k in ("well_assigned","partial","ambiguous_plate"): out.write(f"well\t{w}\t{k}\t{by_well[w][k]}\n")
        for p in sorted(plates_forward):
            for w in sorted(anchors):
                c=p+w
                if by_cell[c]["clean"]:
                    for k in ("clean","a1","adapter_verified","short"): out.write(f"platewell\t{c}\t{k}\t{by_cell[c][k]}\n")
    print(f"Clean PlateWell FASTA: {clean_path}\nWell-only partial FASTA: {partial_path}\nQC: {qc_path}\nSummary: {summary_path}")

if __name__=="__main__": main()
