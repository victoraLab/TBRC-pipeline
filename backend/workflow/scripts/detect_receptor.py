#!/usr/bin/env python3
"""Decide what a library is -- receptor class, loci and species -- and say so.

The panel and the species a run is annotated with are free-text choices made on
the submission form, and until this existed nothing checked either against the
sequences. Both failures look like success:

  * A T-cell library submitted as `igh` does not fail. IgBLAST aligns the reads
    against an immunoglobulin germline database, finds almost nothing, writes a
    table that is syntactically valid and biologically empty, and clonality
    then reports no clones. Both steps exit zero. The run is delivered.

  * A mixed heavy-plus-kappa plate submitted as `igh` does not fail either: every
    kappa sequence is forced onto the best-scoring IGHV gene. On CMP15 that was
    623 of 946 records, all "annotated".

  * A plate whose heavy chains are human, submitted as mouse, aligns at 75%
    identity to mouse germlines and is delivered with 206 calls to one gene.
    The report says completed.

This probes a small deterministic subsample against the broad germline panels
(ig_all and tcr_all) of every species that has them built, and lets the
alignment scores decide. The question "which germline database do these
sequences actually come from" is answered by the same tool that will do the
real annotation, so the probe predicts the outcome rather than approximating it.

Receptor-class outcomes (the `verdict` field):
  agree         the declared panel matches the sequences. Nothing changes.
  widened       the declared single-locus Ig panel misses a locus that is
                present in quantity; ig_all is used instead.
  corrected     the declared class is wrong and the evidence is one-sided.
  mixed         both receptor classes are present in quantity.
  contradicted  the panel is wrong but the mode says report, not act.
  narrow_panel  the panel misses a locus but the mode says report, not act.
Species outcomes (the `species_verdict` field):
  agree         the declared species wins decisively.
  corrected     another species wins decisively; the run is annotated with it.
  contradicted  another species wins decisively but the mode says report only.
  mixed         no species is decisive, typically one locus per species (a
                knock-in). The declared species is kept and the per-locus
                composition says which locus is which.
  unchecked     no other species has its panels built.
Anything the evidence cannot support is `inconclusive`, and inconclusive never
changes the panel or the species. A library that fails to align to anything is
a different problem -- a failed prep, a non-receptor amplicon -- and silently
switching would only hide it.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

# What each selectable panel is looking for. `all_receptors` deliberately has no
# class: it is the answer to a mixed library, never a thing to be corrected away
# from.
PANEL_CLASS = {
    "igh": "Ig",
    "ig_light": "Ig",
    "ig_all": "Ig",
    "ig_b6": "Ig",
    "tra": "TCR",
    "trb": "TCR",
    "tcr_all": "TCR",
    "all_receptors": None,
}

# Which loci a single-locus Ig panel actually carries. A panel not listed here
# is taken to cover its whole class.
PANEL_LOCI = {
    "igh": {"IGH"},
    "ig_light": {"IGK", "IGL"},
    "tra": {"TRA"},
    "trb": {"TRB"},
}

LOCUS_CLASS = {
    "IGH": "Ig", "IGK": "Ig", "IGL": "Ig",
    "TRA": "TCR", "TRB": "TCR", "TRD": "TCR", "TRG": "TCR",
}

# Species whose panels the probe will look for besides the declared one. Only
# those with a built ig_all or tcr_all panel take part; the rest are ignored.
DEFAULT_SPECIES = ("mouse", "human")


def read_fasta(path):
    name, chunks = None, []
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def subsample(records, limit):
    """Evenly spaced, not random: the same input must give the same verdict.

    Spacing across the file rather than taking the head also matters, because
    the final FASTA is ordered by well and the head of it is one corner of one
    plate."""
    total = len(records)
    if total <= limit:
        return records
    step = total / float(limit)
    return [records[int(i * step)] for i in range(limit)]


def run_probe(igblast_script, fasta, workdir, label, species, panel, aux, organism,
              igblast_bin, igblast_data, env):
    out_dir = os.path.join(workdir, label)
    os.makedirs(out_dir, exist_ok=True)
    # The auxiliary file is species-specific. A directory lets the script pick
    # the right one for whichever species is being probed; a file pinned for
    # the declared species would be wrong for every other, so it is only passed
    # through when the probe is for the declared species itself.
    cmd = [igblast_script, fasta, out_dir, label, species, panel, "", "", "",
           aux or "", organism]
    probe_env = dict(env)
    if igblast_data:
        probe_env["IGBLAST_DATA"] = igblast_data
    if igblast_bin:
        probe_env["IGBLAST_BIN"] = igblast_bin
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=probe_env)
    except OSError as exc:
        return None, "probe could not be launched: {}".format(exc)
    result = os.path.join(out_dir, "{}.igblast.tsv".format(label))
    if proc.returncode != 0 or not os.path.exists(result):
        detail = (proc.stderr or "").strip().splitlines()
        return None, detail[-1] if detail else "probe exited {}".format(proc.returncode)
    return result, ""


def parse_airr(path, min_score, max_evalue):
    """sequence_id -> (locus, v_score), keeping only hits worth believing.

    IgBLAST reports a best V hit for almost any input. Measured on this data:
    a real V region scores 337-418 against its own germline set with E-values
    around 1e-115, while 400 nt of random DNA still draws a V call -- at a score
    of 22-33 and an E-value above 1. Without a floor, noise votes. With one, a
    library that is neither Ig nor TCR reports as inconclusive, which is the
    honest answer and leaves the panel alone."""
    hits = {}
    with open(path, newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seq_id = (row.get("sequence_id") or "").strip()
            v_call = (row.get("v_call") or "").strip()
            if not seq_id or not v_call:
                continue
            try:
                score = float(row.get("v_score") or 0.0)
            except ValueError:
                score = 0.0
            if score < min_score:
                continue
            raw_support = (row.get("v_support") or "").strip()
            if raw_support:
                try:
                    if float(raw_support) > max_evalue:
                        continue
                except ValueError:
                    pass
            # The V call names the locus. IgBLAST's own locus column is derived
            # from more than the V hit and disagrees with it on a third of the
            # rows when a single-locus panel is forced onto a mixed plate
            # (CMP15_kappa: 516 of 1,425), so it is not used here.
            locus = v_call.split("*")[0][:3].upper()
            locus = locus if locus in LOCUS_CLASS else ""
            hits[seq_id] = (locus, score)
    return hits


def refine_panel(cls, loci, declared):
    """Pick the narrowest panel the evidence supports.

    If the declared panel is already in the right class, keep it: the user knows
    whether they amplified alpha or beta, and this check is not entitled to
    second-guess a choice that is not wrong. Whether it IS wrong -- a
    single-locus panel on a plate that carries another locus -- is judged
    separately, by uncovered_share."""
    if PANEL_CLASS.get(declared) == cls:
        return declared
    total = sum(loci.values())
    if not total:
        return "ig_all" if cls == "Ig" else "tcr_all"
    share = {k: v / float(total) for k, v in loci.items()}
    if cls == "Ig":
        if share.get("IGH", 0) >= 0.90:
            return "igh"
        if share.get("IGK", 0) + share.get("IGL", 0) >= 0.90:
            return "ig_light"
        return "ig_all"
    if share.get("TRA", 0) >= 0.90:
        return "tra"
    if share.get("TRB", 0) >= 0.90:
        return "trb"
    return "tcr_all"


def uncovered_share(declared, loci):
    """Share of aligned sequences whose locus the declared panel does not carry."""
    covered = PANEL_LOCI.get(declared)
    if not covered:
        return 0.0
    total = sum(loci.values())
    if not total:
        return 0.0
    return sum(v for k, v in loci.items() if k not in covered) / float(total)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--igblast-script", required=True)
    ap.add_argument("--species", default="mouse")
    ap.add_argument("--organism", default="")
    ap.add_argument("--declared-panel", default="")
    ap.add_argument("--aux", default="")
    ap.add_argument("--igblast-bin", default="igblastn")
    ap.add_argument("--igblast-data", default="")
    ap.add_argument("--limit", type=int, default=400)
    ap.add_argument("--report", required=True)
    ap.add_argument("--composition", default="")
    ap.add_argument("--mode", default="auto", choices=["auto", "warn", "off"])
    # Other species to probe. Empty means the built-in list minus the declared
    # species; a species without built panels is skipped without comment.
    ap.add_argument("--species-candidates", default="")
    # A class or species needs this share of the sequences that aligned to
    # anything before it is called the library's class or species.
    ap.add_argument("--decisive-fraction", type=float, default=0.80)
    # A single-locus Ig or TCR panel is widened to its class panel when this
    # share of aligned sequences belongs to a locus it does not carry. On CMP15,
    # submitted as igh, 62% of the plate was kappa; the threshold is set well
    # below that so a plate that is one-fifth light chain is not silently
    # annotated as heavy only.
    ap.add_argument("--widen-fraction", type=float, default=0.10)
    # And this share of the subsample has to align at all for any verdict.
    ap.add_argument("--min-assigned-fraction", type=float, default=0.20)
    # What counts as an alignment. Both are far below anything a real V region
    # produces and far above anything random DNA reaches; see parse_airr.
    ap.add_argument("--min-v-score", type=float, default=100.0)
    ap.add_argument("--max-evalue", type=float, default=1e-20)
    args = ap.parse_args()

    declared = (args.declared_panel or "").strip().lower()
    declared_class = PANEL_CLASS.get(declared)
    declared_species = (args.species or "mouse").strip().lower()
    organism = args.organism or declared_species

    fields = {
        "declared_panel": declared or "(none)",
        "declared_class": declared_class or "unset",
        "mode": args.mode,
        "sequences_tested": 0,
        "aligned": 0,
        "ig_sequences": 0,
        "tcr_sequences": 0,
        "ig_fraction": "",
        "tcr_fraction": "",
        "dominant_locus": "",
        "observed_class": "",
        "recommended_panel": declared,
        "applied_panel": declared,
        "verdict": "inconclusive",
        "detail": "",
        "species_declared": declared_species,
        "species_probed": declared_species,
        "species_observed": "",
        "species_fraction": "",
        "species_verdict": "unchecked",
        "species_detail": "",
        "applied_species": declared_species,
    }
    loci_counts = {}
    locus_species = {}   # (locus, species) -> count

    def emit(code=0):
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        with open(args.report, "w") as handle:
            handle.write("metric\tvalue\n")
            for key, value in fields.items():
                handle.write("{}\t{}\n".format(key, value))
        if args.composition:
            os.makedirs(os.path.dirname(args.composition) or ".", exist_ok=True)
            total = fields["sequences_tested"] or 1
            with open(args.composition, "w") as handle:
                handle.write("locus\tspecies\tsequences\tfraction\n")
                for (locus, species), count in sorted(locus_species.items(),
                                                      key=lambda kv: (-kv[1], kv[0])):
                    handle.write("{}\t{}\t{}\t{:.4f}\n".format(locus, species, count, count / float(total)))
                unassigned = fields["sequences_tested"] - fields["aligned"]
                if unassigned > 0:
                    handle.write("none\t\t{}\t{:.4f}\n".format(unassigned, unassigned / float(total)))
        # The decisions go to stdout on their own lines so the caller can
        # capture them without parsing the table.
        print("APPLIED_PANEL={}".format(fields["applied_panel"]))
        print("RECEPTOR_VERDICT={}".format(fields["verdict"]))
        print("APPLIED_SPECIES={}".format(fields["applied_species"]))
        print("SPECIES_VERDICT={}".format(fields["species_verdict"]))
        sys.exit(code)

    if args.mode == "off":
        fields["verdict"] = "skipped"
        fields["species_verdict"] = "skipped"
        fields["detail"] = "receptor check disabled for this run"
        emit()

    if not os.path.exists(args.fasta):
        fields["detail"] = "no assembled FASTA to test"
        emit()

    records = list(read_fasta(args.fasta))
    records = [(n, s) for n, s in records if len(s) >= 100]
    if not records:
        fields["detail"] = "no sequence long enough to align"
        emit()

    sample = subsample(records, args.limit)
    fields["sequences_tested"] = len(sample)

    if args.species_candidates.strip():
        candidates = [s.strip().lower() for s in args.species_candidates.split(",") if s.strip()]
    else:
        candidates = list(DEFAULT_SPECIES)
    species_list = [declared_species] + [s for s in candidates if s != declared_species]

    workdir = tempfile.mkdtemp(prefix="tbrc_receptor_")
    try:
        probe_fa = os.path.join(workdir, "probe.fasta")
        with open(probe_fa, "w") as handle:
            for name, seq in sample:
                handle.write(">{}\n{}\n".format(name, seq))

        env = os.environ.copy()
        # (species, class) -> hits. The declared species keeps the submitted
        # organism and auxiliary file; every other species is probed under its
        # own name so the script resolves its own auxiliary data.
        tables = {}
        errors = {}
        for species in species_list:
            own = species == declared_species
            for panel, cls in (("ig_all", "Ig"), ("tcr_all", "TCR")):
                tsv, err = run_probe(args.igblast_script, probe_fa, workdir,
                                     "{}_{}".format(species, panel), species, panel,
                                     args.aux if own else (args.aux if os.path.isdir(args.aux or "") else ""),
                                     organism if own else species,
                                     args.igblast_bin, args.igblast_data, env)
                if tsv is None:
                    errors[(species, cls)] = err
                    continue
                tables[(species, cls)] = parse_airr(tsv, args.min_v_score, args.max_evalue)

        probed_species = sorted({sp for sp, _ in tables})
        fields["species_probed"] = ",".join(probed_species) or declared_species

        if not any(sp == declared_species for sp, _ in tables):
            fields["detail"] = "neither germline panel could be probed for {} ({})".format(
                declared_species,
                errors.get((declared_species, "Ig")) or errors.get((declared_species, "TCR")) or "unknown")
            emit()

        # One vote per sequence, to whichever probe aligned it best. A
        # per-sequence contest rather than summary rates, because a mixed
        # library -- two classes, or two species -- has to be visible as a
        # mixture. Ties go to the declared species, then to Ig.
        species_counts = {}
        for name, _ in sample:
            best = None
            for (species, cls), hits in tables.items():
                hit = hits.get(name)
                if hit is None:
                    continue
                key = (hit[1], species == declared_species, cls == "Ig")
                if best is None or key > best[0]:
                    best = (key, species, cls, hit[0])
            if best is None:
                continue
            _, species, cls, locus = best
            locus = locus or ("IG?" if cls == "Ig" else "TR?")
            loci_counts[locus] = loci_counts.get(locus, 0) + 1
            species_counts[species] = species_counts.get(species, 0) + 1
            locus_species[(locus, species)] = locus_species.get((locus, species), 0) + 1
            fields["aligned"] += 1
            fields["ig_sequences" if cls == "Ig" else "tcr_sequences"] += 1
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    tested = fields["sequences_tested"]
    aligned = fields["aligned"]
    if aligned < max(1, int(args.min_assigned_fraction * tested)):
        fields["detail"] = (
            "only {}/{} sequences produced a credible germline alignment (V score "
            ">= {:g}, E <= {:g}); too few to judge the receptor type. Check the "
            "species, and whether this library is a receptor amplicon at all.".format(
                aligned, tested, args.min_v_score, args.max_evalue))
        emit()

    # ---------------------------------------------------------------- species
    def locus_species_text():
        parts = []
        by_locus = {}
        for (locus, species), count in locus_species.items():
            by_locus.setdefault(locus, {})[species] = count
        for locus in sorted(by_locus, key=lambda k: -sum(by_locus[k].values())):
            counts = by_locus[locus]
            total = sum(counts.values())
            top = max(counts.items(), key=lambda kv: kv[1])
            parts.append("{}: {}% {} ({} sequences)".format(
                locus, int(round(100.0 * top[1] / total)), top[0], total))
        return "; ".join(parts)

    if len(probed_species) < 2:
        fields["species_verdict"] = "unchecked"
        fields["species_observed"] = declared_species
        fields["species_detail"] = (
            "only the {} germline panels are built on this server, so the species "
            "could not be checked against another.".format(declared_species))
    else:
        declared_share = species_counts.get(declared_species, 0) / float(aligned)
        top_species, top_count = max(species_counts.items(), key=lambda kv: kv[1])
        top_share = top_count / float(aligned)
        fields["species_observed"] = top_species
        fields["species_fraction"] = "{:.4f}".format(declared_share)
        if declared_share >= args.decisive_fraction:
            fields["species_verdict"] = "agree"
            fields["species_detail"] = "{}% of the {} sequences that aligned score best against {} germlines.".format(
                int(round(100 * declared_share)), aligned, declared_species)
        elif top_share >= args.decisive_fraction:
            fields["species_verdict"] = "corrected" if args.mode == "auto" else "contradicted"
            if args.mode == "auto":
                fields["applied_species"] = top_species
            fields["species_detail"] = (
                "the run was submitted as {}, but {}% of the {} sequences that aligned "
                "score best against {} germlines. {}".format(
                    declared_species, int(round(100 * top_share)), aligned, top_species,
                    "Annotation and clonality ran against {} germlines instead.".format(top_species)
                    if args.mode == "auto" else
                    "The species was left as submitted, so V identity will be far below "
                    "what the sequences deserve and the gene calls will be wrong."))
        else:
            fields["species_verdict"] = "mixed"
            fields["species_detail"] = (
                "no single species accounts for {}% of the sequences ({}). The declared "
                "species, {}, was kept: IgBLAST annotates against one germline set, so "
                "the loci from the other species carry low V identity and unreliable "
                "gene calls in this archive and need their own annotation.".format(
                    int(round(100 * args.decisive_fraction)), locus_species_text(), declared_species))

    # ---------------------------------------------------------------- class
    ig_frac = fields["ig_sequences"] / float(aligned)
    tcr_frac = fields["tcr_sequences"] / float(aligned)
    fields["ig_fraction"] = "{:.4f}".format(ig_frac)
    fields["tcr_fraction"] = "{:.4f}".format(tcr_frac)
    if loci_counts:
        fields["dominant_locus"] = max(loci_counts.items(), key=lambda kv: kv[1])[0]

    if max(ig_frac, tcr_frac) < args.decisive_fraction:
        fields["observed_class"] = "mixed"
        fields["verdict"] = "mixed"
        fields["recommended_panel"] = "all_receptors"
        fields["detail"] = (
            "{:.0f}% immunoglobulin and {:.0f}% T-cell receptor among the {} "
            "sequences that aligned. Neither class is dominant, so a single-class "
            "panel would discard the other half.".format(
                100 * ig_frac, 100 * tcr_frac, aligned))
        if args.mode == "auto" and declared != "all_receptors":
            fields["applied_panel"] = "all_receptors"
        emit()

    observed = "Ig" if ig_frac >= tcr_frac else "TCR"
    fields["observed_class"] = observed
    recommended = refine_panel(observed, loci_counts, declared)
    fields["recommended_panel"] = recommended

    if declared_class is None:
        # all_receptors, or nothing declared. Both already cover the data.
        fields["verdict"] = "agree"
        fields["detail"] = (
            "{}% of aligned sequences are {}; the declared panel covers it.".format(
                int(round(100 * max(ig_frac, tcr_frac))), observed))
        emit()

    if declared_class == observed:
        missing = uncovered_share(declared, loci_counts)
        if missing >= args.widen_fraction:
            wide = "ig_all" if observed == "Ig" else "tcr_all"
            fields["recommended_panel"] = wide
            fields["verdict"] = "widened" if args.mode == "auto" else "narrow_panel"
            if args.mode == "auto":
                fields["applied_panel"] = wide
            fields["detail"] = (
                "the run was submitted as {}, but {}% of the {} sequences that aligned "
                "belong to a locus that panel does not carry ({}). {}".format(
                    declared, int(round(100 * missing)), aligned,
                    ", ".join("{} {}".format(v, k) for k, v in sorted(loci_counts.items(), key=lambda kv: -kv[1])),
                    "Annotation and clonality ran with {} instead, so every locus on the plate "
                    "is called against its own germlines.".format(wide)
                    if args.mode == "auto" else
                    "The panel was left as submitted, so those sequences are forced onto "
                    "the wrong locus's genes."))
            emit()
        fields["verdict"] = "agree"
        fields["detail"] = (
            "{}% of the {} sequences that aligned are {}, matching the "
            "declared panel.".format(int(round(100 * max(ig_frac, tcr_frac))),
                                     aligned, observed))
        emit()

    fields["verdict"] = "corrected" if args.mode == "auto" else "contradicted"
    if args.mode == "auto":
        fields["applied_panel"] = recommended
    fields["detail"] = (
        "the run was submitted as {} ({}), but {}% of the {} sequences that "
        "aligned are {} -- mostly {}. {}".format(
            declared, declared_class, int(round(100 * max(ig_frac, tcr_frac))),
            aligned, observed, fields["dominant_locus"] or "unclear",
            "Annotation and clonality ran with {} instead.".format(recommended)
            if args.mode == "auto" else
            "The panel was left as submitted, so annotation will be largely empty."))
    emit()


if __name__ == "__main__":
    main()
