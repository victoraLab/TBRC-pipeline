#!/usr/bin/env python3
"""Run the vendored nanopore demultiplexer on several input FASTQs at once.

demux_nanopore_wells.py reads every input in one process, one record at a time,
and does its barcode and primer search in pure Python: a full MinION library
(40 FASTQ chunks, 8.1 GB, 14 M reads) kept one core busy for 56 min while 31 of
the workstation's 32 threads idled -- 56 of the run's 60 min. Every read is
demultiplexed on its own, so nothing is lost by splitting the work: this runs the
vendored script once per input file, --jobs at a time, and merges the per-file
outputs back in input order. On that library, with --jobs 8, it took 6.6 min, and
the merged FASTAs, read table and summary were byte-identical to the serial run.

The vendored script is driven, not modified, apart from one reporting flag
(--summary-all-platewells, see its header), so it can still be re-vendored from
upstream. Takes the vendored script's own arguments plus --jobs; with one input
or --jobs 1 it simply runs the vendored script as before.
"""

import argparse, collections, importlib.util, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

VENDORED = Path(__file__).with_name("demux_nanopore_wells.py")
_spec = importlib.util.spec_from_file_location("demux_nanopore_wells", VENDORED)
DMX = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(DMX)

# The summary's row layout, as demux_nanopore_wells.py writes it.
WELL_METRICS = ("well_assigned", "partial", "ambiguous_plate")
CELL_METRICS = ("clean", "a1", "adapter_verified", "short")
SUFFIXES = (".platewell_clean.fasta", ".well_only_partial.fasta",
            ".platewell_read_qc.tsv", ".platewell_summary.tsv")


def vendored_args(a, fastq, extra, prefix):
    """The vendored script's command line for these inputs; options left unset
    keep the vendored defaults."""
    cmd = [sys.executable, str(VENDORED), fastq, a.well_anchors, a.plate_barcodes, prefix]
    for flag, val in (("--sample", a.sample), ("--well-barcode-mm", a.well_barcode_mm),
                      ("--plate-barcode-mm", a.plate_barcode_mm), ("--primer-mm", a.primer_mm),
                      ("--terminal-window", a.terminal_window),
                      ("--min-insert-length", a.min_insert_length)):
        if val is not None: cmd += [flag, str(val)]
    for p in a.inner_reverse_primer or []: cmd += ["--inner-reverse-primer", p]
    for p in extra: cmd += ["--extra-fastq", p]
    return cmd


def read_summary(path):
    rows = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    with open(path) as h:
        if h.readline().rstrip("\n") != "section\tid\tmetric\tvalue":
            raise SystemExit(f"unexpected summary header in {path}")
        for line in h:
            section, ident, metric, value = line.rstrip("\n").split("\t")
            rows[section][ident][metric] += int(value)
    return rows


def merge(a, parts, inputs):
    out = a.prefix
    for suffix in SUFFIXES[:2]:
        with open(out + suffix, "wb") as dst:
            for p in parts:
                with open(p + suffix, "rb") as src: shutil.copyfileobj(src, dst, 1 << 20)
    with open(out + SUFFIXES[2], "wb") as dst:
        header = None
        for p in parts:
            with open(p + SUFFIXES[2], "rb") as src:
                first = src.readline()
                if header is None: header = first; dst.write(first)
                elif first != header: raise SystemExit(f"read-table header differs in {p}")
                shutil.copyfileobj(src, dst, 1 << 20)

    tot = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for p in parts:
        for section, ids in read_summary(p + SUFFIXES[3]).items():
            for ident, metrics in ids.items(): tot[section][ident].update(metrics)
    anchors, plates = DMX.table(a.well_anchors), DMX.table(a.plate_barcodes)
    with open(out + SUFFIXES[3], "w") as h:
        h.write("section\tid\tmetric\tvalue\n")
        for k, v in sorted(tot["global"]["ALL"].items()): h.write(f"global\tALL\t{k}\t{v}\n")
        for source in sorted(tot["source"]):
            for k, v in sorted(tot["source"][source].items()): h.write(f"source\t{source}\t{k}\t{v}\n")
        for w in sorted(anchors):
            for k in WELL_METRICS: h.write(f"well\t{w}\t{k}\t{tot['well'][w][k]}\n")
        for p in sorted(plates):
            for w in sorted(anchors):
                c = p + w
                if tot["platewell"][c]["clean"]:
                    for k in CELL_METRICS: h.write(f"platewell\t{c}\t{k}\t{tot['platewell'][c][k]}\n")
    for suffix, label in zip(SUFFIXES, ("Clean PlateWell FASTA", "Well-only partial FASTA", "QC", "Summary")):
        print(f"{label}: {out}{suffix}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq"); ap.add_argument("well_anchors"); ap.add_argument("plate_barcodes"); ap.add_argument("prefix")
    ap.add_argument("--extra-fastq", action="append", default=[])
    ap.add_argument("--sample"); ap.add_argument("--well-barcode-mm", type=int); ap.add_argument("--plate-barcode-mm", type=int)
    ap.add_argument("--primer-mm", type=int); ap.add_argument("--terminal-window", type=int)
    ap.add_argument("--min-insert-length", type=int); ap.add_argument("--inner-reverse-primer", action="append")
    ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()

    inputs = [a.fastq] + a.extra_fastq
    # The vendored script names each read's source after its file, numbering
    # repeats of one name in input order; split one file per process, two files
    # with one name would both be the first. Rare (inputs come from one folder),
    # so such a run simply stays serial.
    names = [DMX.source_id(p) for p in inputs]
    if a.jobs <= 1 or len(inputs) == 1 or len(set(names)) != len(names):
        sys.exit(subprocess.call(vendored_args(a, a.fastq, a.extra_fastq, a.prefix)))
    for p in inputs:
        if not Path(p).is_file(): raise SystemExit(f"Missing FASTQ: {p}")

    work = Path(a.prefix + ".parts")
    if work.exists(): shutil.rmtree(work)
    work.mkdir(parents=True)
    parts = [str(work / f"part{i:04d}") for i in range(len(inputs))]

    def one(i):
        cmd = vendored_args(a, inputs[i], [], parts[i]) + ["--summary-all-platewells"]
        with open(parts[i] + ".log", "w") as log:
            return subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)

    t0 = time.time()
    # Largest files first, so the last process started is a short one.
    order = sorted(range(len(inputs)), key=lambda i: -Path(inputs[i]).stat().st_size)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        codes = dict(zip(order, pool.map(one, order)))
    failed = [i for i in range(len(inputs)) if codes[i] != 0]
    if failed:
        for i in failed:
            sys.stderr.write(f"demux failed on {inputs[i]}:\n" + Path(parts[i] + ".log").read_text())
        raise SystemExit(1)
    merge(a, parts, inputs)
    shutil.rmtree(work)
    print(f"[parallel] {len(inputs)} inputs, {a.jobs} jobs, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
