#!/usr/bin/env python3
"""V/J gene concordance between the two pipelines on the wells they share.

The productive-rate metric asks whether a delivered sequence is callable. This
asks the harder question: when both pipelines deliver the same well, do they
agree on WHICH receptor it holds? That is what matters if an old dataset is
reprocessed and the gene call moves.

Compared at gene level on the rank-1 (dominant) sequence of each shared well.
IGHV3-6CHI1-64 is collapsed to IGHV3-6 first -- the same gene under two names,
one used by the custom reference and one by IMGT.
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

# One gene, two names. Not collapsing it inflates apparent disagreement.
ALIASES = {"IGHV3-6CHI1-64": "IGHV3-6"}


def gene_of(call):
    """Top gene from an AIRR v_call/j_call, allele and ambiguity stripped."""
    if not call:
        return ""
    top = call.split(",")[0].strip()
    top = top.split("*")[0]
    return ALIASES.get(top, top)


def load(airr, sample):
    """rank-1 record per well: (v gene, j gene, cdr3, read count)."""
    out = {}
    if not os.path.exists(airr):
        return out
    with open(airr) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row.get("sequence_id", "")
            m = re.match(rf"^{re.escape(sample)}(.+)_(\d+)-(\d+)$", sid)
            if not m or m.group(2) != "1":
                continue
            out[m.group(1)] = (
                gene_of(row.get("v_call")),
                gene_of(row.get("j_call")),
                (row.get("cdr3_aa") or "").strip(),
                int(m.group(3)),
            )
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bench", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    SPEC = [("CMP12", "CMP12", "IgH"), ("CMP13", "CMP13", "IgH"),
            ("CMP14new", "CMP14", "IgH"), ("CMP10", "CMP10", "IgK"),
            ("MS_biotin_202604", "Miseq_1_MS", "TRB"), ("ABGRES", "AB221104", "TRB"),
            ("parabiosis_trial_rp", "RP230620", "TRB")]

    results = []
    for ds, samp, typ in SPEC:
        base = f"{args.bench}/{ds}"
        C = load(f"{base}/classic_luka/results/bench/{ds}/{samp}/classic_luka.airr.tsv", samp)
        V = load(f"{base}/tiago_vsearch/results/bench/{ds}/{samp}/tiago_vsearch.airr.tsv", samp)
        shared = sorted(set(C) & set(V))
        if not shared:
            continue

        vsame = jsame = both = cdr3same = 0
        vcomp = jcomp = cdr3comp = 0
        disagree = collections.Counter()
        for w in shared:
            cv, cj, cc, _ = C[w]
            vv, vj, vc, _ = V[w]
            if cv and vv:
                vcomp += 1
                if cv == vv:
                    vsame += 1
                else:
                    disagree[f"{cv} -> {vv}"] += 1
            if cj and vj:
                jcomp += 1
                if cj == vj:
                    jsame += 1
            if cv and vv and cj and vj and cv == vv and cj == vj:
                both += 1
            if cc and vc:
                cdr3comp += 1
                if cc == vc:
                    cdr3same += 1

        results.append({
            "dataset": ds, "type": typ,
            "classic_wells": len(C), "vsearch_wells": len(V), "shared": len(shared),
            "v_comparable": vcomp, "v_agree": vsame,
            "v_pct": round(100 * vsame / vcomp, 2) if vcomp else None,
            "j_comparable": jcomp, "j_agree": jsame,
            "j_pct": round(100 * jsame / jcomp, 2) if jcomp else None,
            "vj_agree": both,
            "vj_pct": round(100 * both / vcomp, 2) if vcomp else None,
            "cdr3_comparable": cdr3comp, "cdr3_agree": cdr3same,
            "cdr3_pct": round(100 * cdr3same / cdr3comp, 2) if cdr3comp else None,
            "top_v_disagreements": disagree.most_common(5),
        })

    json.dump(results, open(args.out, "w"), indent=1)

    print(f"  {'dataset':<21}{'type':<6}{'shared':>8}{'V agree':>10}{'J agree':>10}{'V+J':>8}{'CDR3':>9}")
    print("  " + "-" * 72)
    for r in results:
        print(f"  {r['dataset']:<21}{r['type']:<6}{r['shared']:>8,}"
              f"{r['v_pct']:>9.1f}%{r['j_pct']:>9.1f}%{r['vj_pct']:>7.1f}%{r['cdr3_pct']:>8.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
