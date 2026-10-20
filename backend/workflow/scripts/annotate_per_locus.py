#!/usr/bin/env python3
"""Annotate a plate whose loci need different germline sets, into one table.

IgBLAST takes one V, one D and one J database per run, and the receptor check
could only pick one species for a whole plate. Neither fits a plate whose chains
come from different germline sets. On CMP15 the kappa chains are mouse and the
heavy chains are a human IGHV1-69 rearranged to mouse D and J: annotated as
mouse, the heavy V matched at 75%; annotated as human, the mouse JH matched at
89% and FR4 looked 13% mutated. Both went out as "completed".

detect_receptor.py now writes a plan, locus -> V species / D-J species, when a
plate needs one. This follows it:

  1. Annotate every sequence as the declared species, as a plain run would, and
     as every other V species in the plan. Each sequence's locus comes from
     whichever of those V calls scores best. From the declared table alone, 34
     of CMP15's heavy chains read as another locus: against mouse germlines a
     human V aligns at 75% and can lose to nothing at all.
  2. A locus whose V and J species agree takes its rows from that species'
     table. A locus that needs both is run again: V database and FR/CDR
     definitions (-organism) from the V species, D and J databases and J frames
     (the auxiliary file) from the D-J species. On CMP15 that gives the heavy
     chains 98% median J identity and FR4 WGQGTLVTVSA (mouse JH3).
  3. Write one table, one row per query in input order like a plain IgBLAST
     run, with a germline_set column naming the set each row was called against.

Exit status is run_igblast.sh's: 0, 64 for "not set up on this host", anything
else for a failure.

Usage: annotate_per_locus.py --fasta F --out-dir D --sample S --plan IGH:human/mouse,IGK:mouse/mouse
         --declared-species mouse --panel ig_all --igblast-script run_igblast.sh
         [--organism mouse] [--aux FILE_OR_DIR]
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

csv.field_size_limit(sys.maxsize)


def read_fasta(path):
    name, chunks = None, []
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                # IgBLAST's sequence_id is the header up to the first space.
                name, chunks = (line[1:].split() or [""])[0], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse_plan(text):
    plan = {}
    for item in (text or "").split(","):
        item = item.strip()
        if not item:
            continue
        locus, _, sets = item.partition(":")
        v_sp, _, j_sp = sets.partition("/")
        if not (locus and v_sp and j_sp):
            raise SystemExit("annotate_per_locus: cannot read plan entry {!r}".format(item))
        plan[locus.strip().upper()] = (v_sp.strip().lower(), j_sp.strip().lower())
    return plan


def read_table(path):
    with open(path, newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {(row.get("sequence_id") or "").strip(): row for row in reader}
        return reader.fieldnames or [], rows


def database(work_root, species, panel, segment):
    """The panel's database for one segment, or the class-wide panel's.

    Same two spellings run_igblast.sh accepts (<sp>_<panel>_<seg>, and
    <sp>_<panel without ig_>_<seg> for hand-built panels)."""
    fallback = "tcr_all" if panel in ("tra", "trb", "tcr_all") else "ig_all"
    for name in dict.fromkeys((panel, fallback)):
        for stem in ("{}_{}_{}".format(species, name, segment),
                     "{}_{}_{}".format(species, name[3:] if name.startswith("ig_") else name, segment)):
            path = os.path.join(work_root, species, name, stem)
            if os.path.exists(path + ".nhr"):
                return path
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--declared-species", required=True)
    ap.add_argument("--panel", default="ig_all")
    ap.add_argument("--organism", default="")
    ap.add_argument("--aux", default="")
    ap.add_argument("--igblast-script", required=True)
    args = ap.parse_args()

    declared = args.declared_species.strip().lower()
    plan = parse_plan(args.plan)
    script = os.path.abspath(args.igblast_script)
    work_root = os.path.join(os.path.dirname(script), "..", "..", "igblast", "work")

    def aux_for(species):
        # A file pinned for the declared species would give another species
        # wrong J frames; a directory lets run_igblast.sh pick by species.
        if species == declared or os.path.isdir(args.aux or ""):
            return args.aux or ""
        return os.path.dirname(args.aux) if os.path.isfile(args.aux or "") else ""

    def run(fasta, label, species, organism, db_v="", db_d="", db_j=""):
        out = os.path.join(workdir, label)
        cmd = [script, fasta, out, label, species, args.panel, db_v, db_d, db_j,
               aux_for(species), organism]
        proc = subprocess.run(cmd, stdout=sys.stderr)
        if proc.returncode != 0:
            sys.exit(proc.returncode)
        return read_table(os.path.join(out, "{}.igblast.tsv".format(label)))

    def v_score(row):
        try:
            return float(row.get("v_score") or 0.0)
        except ValueError:
            return 0.0

    records = list(read_fasta(args.fasta))
    workdir = tempfile.mkdtemp(prefix="tbrc_per_locus_")
    try:
        columns, base_rows = run(args.fasta, "declared", declared, args.organism or declared)
        species_rows = {declared: base_rows}
        for v_sp in sorted({v for v, _ in plan.values()} - {declared}):
            species_rows[v_sp] = run(args.fasta, v_sp, v_sp, v_sp)[1]
        locus_of = {}
        for name, _ in records:
            # Best V across the species, ties to the declared one.
            best = max(((v_score(rows.get(name, {})), sp == declared, rows.get(name, {}))
                        for sp, rows in species_rows.items()), key=lambda t: t[:2])
            v_call = (best[2].get("v_call") or "").strip()
            locus_of[name] = v_call.split("*")[0][:3].upper() if v_call else ""

        own_rows, labels = {}, {}
        for locus, (v_sp, j_sp) in sorted(plan.items()):
            label = v_sp if v_sp == j_sp else "{} V; {} D, J".format(v_sp, j_sp)
            names = [n for n, _ in records if locus_of.get(n) == locus]
            if not names or (v_sp, j_sp) == (declared, declared):
                continue
            if v_sp == j_sp:
                rows = species_rows[v_sp]
            else:
                subset = os.path.join(workdir, "{}.fasta".format(locus))
                with open(subset, "w") as handle:
                    for name, seq in records:
                        if locus_of.get(name) == locus:
                            handle.write(">{}\n{}\n".format(name, seq))
                db = {seg: database(work_root, sp, args.panel, seg)
                      for seg, sp in (("V", v_sp), ("D", j_sp), ("J", j_sp))}
                if not db["V"] or not db["J"]:
                    print("annotate_per_locus: no {} V or {} J database built for panel {}; "
                          "{} keeps the {} annotation.".format(v_sp, j_sp, args.panel, locus, declared),
                          file=sys.stderr)
                    continue
                # -organism is the V species: it selects the FR/CDR definitions
                # of the V. The J frames come from the auxiliary file, which
                # run_igblast.sh picks by its species argument, here the D-J one.
                _, rows = run(subset, locus, j_sp, v_sp, db["V"], db["D"], db["J"])
            for name in names:
                if name in rows:
                    own_rows[name] = rows[name]
                    labels[name] = label
            print("annotate_per_locus: {} {} sequences annotated as {}.".format(
                len(names), locus, label), file=sys.stderr)

        os.makedirs(args.out_dir, exist_ok=True)
        out_path = os.path.join(args.out_dir, "{}.igblast.tsv".format(args.sample))
        fieldnames = [c for c in columns if c != "germline_set"] + ["germline_set"]
        with open(out_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for name, _ in records:
                row = own_rows.get(name) or base_rows.get(name)
                if row is None:
                    continue
                row = dict(row)
                row["germline_set"] = labels.get(name, declared)
                writer.writerow(row)
        print("IGBlast output: {}".format(out_path))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
