#!/usr/bin/env python3
"""Clone assignment and incidence counting for MIBS / sulkseq libraries.

MIBS sorts ~100 cells into each of 16-32 wells and asks a different question
from every other pipeline here. Clone size is not read depth: it is
*incidence*, the number of wells a clone turns up in. A PCR artefact or a
chimera appears in one well; a clone expanded enough to be sampled repeatedly
appears in several. That is the method's error control, and it is why the
published filters are "at least 6 reads in some well" and "present in at least
2 wells of the same sample" rather than any depth threshold.

Two consequences shape this script:

  Clones are assigned ACROSS wells within a sample -- the opposite direction
  from every other clonality step in TBRC, which groups within a run.

  Every cluster in a well is a candidate, not just the dominant one. A well
  holding 100 cells has no dominant clone to find.

The clone rule is the paper's: same V gene, same J gene, same junction length,
junction Hamming distance <= 4. Validated against the 586,821 annotated
sequences the MIBS repository ships with their published clone calls; rerun
that comparison at any time with `--validate-against <splited_database dir>`.
Of every pair of sequences this rule places in one clone, the
published calls agree on 100.00%; of the pairs the published calls join, this
recovers 99.78%. After the published filters the two agree on the number that
matters: identical maximum wells per clone (32) and an identical wells-per-clone
distribution at the median, 90th and 99th percentile.

The one detail that is easy to get wrong: IgBLAST reports ties as a
comma-separated call, and 12% of MIBS V calls are ambiguous. Matching on the
first call alone splits clones the published analysis keeps together -- 75 of
them in one sample. Genes must be compared as sets.
"""

import argparse
import collections
import csv
import os
import re
import sys

RECEPTOR_PREFIXES = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
# >{sample}{plate}{row}{col}_{rank}-{reads}, written by rank_wells.py. The well
# token is whatever sits between the sample name and the underscore.
HEADER_RE = re.compile(r"^(?P<well>.+?)_(?P<rank>\d+)-(?P<reads>\d+)$")


def gene_set(call):
    """Every gene named in a call, as a set."""
    out = set()
    for part in (call or "").split(","):
        for token in part.split():
            upper = token.upper()
            if upper.startswith(RECEPTOR_PREFIXES):
                out.add(token.split("*")[0])
                break
    return out


# A well's identity is plate, row and column, but the token the demultiplexer
# builds concatenates the two barcode names in whichever order the barcode set
# puts them: set B names columns in BC1 and plate+row in BC2, giving "01P05B",
# while set T is the other way round and gives "P05B01". Both mean the same
# well. Everything downstream -- the plate map, the published PLATE_ID, the
# plate figures -- expects plate, row, column, so normalise once here.
# The plate letter is matched in either case for the same reason the row is:
# these patterns read a plate map typed by hand as well as a well name built by
# the demultiplexer, and p01a01 is what somebody types. PLATE_ONLY was already
# case-insensitive, so a bare "p1" was understood and "p01a01" was not.
WELL_CANONICAL = re.compile(r"^(?P<plate>[Pp]\d{1,2})(?P<row>[A-Ha-h])(?P<col>\d{1,2})$")
WELL_COL_FIRST = re.compile(r"^(?P<col>\d{1,2})(?P<plate>[Pp]\d{1,2})(?P<row>[A-Ha-h])$")
# Barcode set B2 names its wells row+column and its plates separately, giving
# A01P01 -- the plate token last. Same well, third spelling.
WELL_PLATE_LAST = re.compile(r"^(?P<row>[A-Ha-h])(?P<col>\d{1,2})(?P<plate>[Pp]\d{1,2})$")
# A library run on a single plate carries no plate barcode at all.
WELL_NO_PLATE = re.compile(r"^(?P<row>[A-Ha-h])(?P<col>\d{1,2})$")
WELL_NO_PLATE_COL_FIRST = re.compile(r"^(?P<col>\d{1,2})(?P<row>[A-Ha-h])$")
PLATE_ONLY = re.compile(r"^P\d{1,2}$", re.I)


def canonical_plate(token):
    """P05, whether it arrived as P5 or P05.

    The demultiplexer always pads, because every shipped barcode set names its
    plates P01 upwards. A plate map is typed by hand and P1 is what people
    write, and the unpadded key matched nothing: an entire plate then fell out
    of the design and its wells were reported as a phantom '(unmapped wells)'
    sample sitting beside the real ones."""
    text = (token or "").strip().upper()
    if not PLATE_ONLY.match(text):
        return text
    return "P{:02d}".format(int(text[1:]))


def canonical_well(token):
    """plate + row + column, e.g. P05B04. Unrecognised tokens pass through."""
    for pattern in (WELL_CANONICAL, WELL_COL_FIRST, WELL_PLATE_LAST):
        match = pattern.match(token)
        if match:
            return "{}{}{:02d}".format(canonical_plate(match.group("plate")),
                                       match.group("row").upper(),
                                       int(match.group("col")))
    for pattern in (WELL_NO_PLATE, WELL_NO_PLATE_COL_FIRST):
        match = pattern.match(token)
        if match:
            return "{}{:02d}".format(match.group("row").upper(), int(match.group("col")))
    return token


def parse_header(seq_id, sample):
    """(well, reads) from a delivered header, or (None, None)."""
    match = HEADER_RE.match(seq_id)
    if not match:
        return None, None
    well = match.group("well")
    if sample and well.startswith(sample):
        well = well[len(sample):]
    if not well:
        return None, None
    return canonical_well(well), int(match.group("reads"))


def _merge_by_gene(records, field):
    """Group records transitively by shared gene -- Change-O's other reading.

    A shares a gene with B and B with C puts all three in one group even when A
    and C share nothing. Membership is then decided by the group, and junction
    distance alone links within it."""
    parent = list(range(len(records)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    first = {}
    for index, rec in enumerate(records):
        for gene in rec[field]:
            if gene in first:
                a, b = find(first[gene]), find(index)
                if a != b:
                    parent[a] = b
            else:
                first[gene] = index
    return [find(i) for i in range(len(records))]


def assign_clones_transitive(records, max_distance):
    """Gene grouping first, then single linkage on junction distance alone."""
    v_group = _merge_by_gene(records, "v_genes")
    j_group = _merge_by_gene(records, "j_genes")
    groups = collections.defaultdict(list)
    for index, rec in enumerate(records):
        groups[(v_group[index], j_group[index], len(rec["junction"]))].append(index)

    labels = [None] * len(records)
    next_id = 1
    for members in groups.values():
        unique = {}
        for member in members:
            unique.setdefault(records[member]["junction"], []).append(member)
        keys = list(unique)
        parent = list(range(len(keys)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                distance = 0
                for a, b in zip(keys[i], keys[j]):
                    if a != b:
                        distance += 1
                        if distance > max_distance:
                            break
                if distance <= max_distance:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb

        local = {}
        for i, key in enumerate(keys):
            root = find(i)
            if root not in local:
                local[root] = next_id
                next_id += 1
            for member in unique[key]:
                labels[member] = local[root]
    return labels


def assign_clones(records, max_distance):
    """Single linkage over: shared V gene, shared J gene, equal junction length,
    Hamming <= max_distance.

    Each record is indexed under every (V, J) pair its calls allow, so any two
    records sharing a V and a J meet in some bucket. That is the same relation
    as comparing all pairs, without the quadratic over the whole sample."""
    buckets = collections.defaultdict(list)
    for index, rec in enumerate(records):
        for v_gene in rec["v_genes"]:
            for j_gene in rec["j_genes"]:
                buckets[(v_gene, j_gene, len(rec["junction"]))].append(index)

    parent = list(range(len(records)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_a] = root_b

    for members in buckets.values():
        # Identical junctions are distance zero; collapse before the pairwise
        # pass so it runs over distinct sequences rather than over reads.
        unique = {}
        for member in members:
            unique.setdefault(records[member]["junction"], []).append(member)
        keys = list(unique)
        for key in keys:
            head = unique[key][0]
            for other in unique[key][1:]:
                union(head, other)
        for i in range(len(keys)):
            left = keys[i]
            for j in range(i + 1, len(keys)):
                right = keys[j]
                distance = 0
                for a, b in zip(left, right):
                    if a != b:
                        distance += 1
                        if distance > max_distance:
                            break
                if distance <= max_distance:
                    union(unique[left][0], unique[right][0])

    labels, seen, next_id = [], {}, 1
    for index in range(len(records)):
        root = find(index)
        if root not in seen:
            seen[root] = next_id
            next_id += 1
        labels.append(seen[root])
    return labels


def load_airr(path, sample):
    """Annotated sequences with a junction, one row per delivered sequence."""
    records, skipped = [], collections.Counter()
    with open(path, newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seq_id = (row.get("sequence_id") or "").strip()
            junction = (row.get("junction") or "").strip().upper()
            v_genes = gene_set(row.get("v_call"))
            j_genes = gene_set(row.get("j_call"))
            if not seq_id:
                skipped["no_sequence_id"] += 1
                continue
            if not junction:
                skipped["no_junction"] += 1
                continue
            if not v_genes or not j_genes:
                skipped["no_v_or_j_call"] += 1
                continue
            well, reads = parse_header(seq_id, sample)
            if well is None:
                skipped["unparsable_header"] += 1
                continue
            records.append({
                "sequence_id": seq_id,
                "well": well,
                "reads": reads,
                "junction": junction,
                "junction_aa": (row.get("junction_aa") or "").strip(),
                "v_call": (row.get("v_call") or "").strip(),
                "j_call": (row.get("j_call") or "").strip(),
                "d_call": (row.get("d_call") or "").strip(),
                "v_genes": v_genes,
                "j_genes": j_genes,
                "productive": (row.get("productive") or "").strip().upper() in ("T", "TRUE"),
                "sequence": (row.get("sequence") or "").strip(),
                "v_identity": (row.get("v_identity") or "").strip(),
            })
    return records, skipped


def load_platemap(path):
    """well -> {sample, mouse, tissue, group}. Absent means one sample.

    Three shapes are accepted, because the experimental design is not something
    the pipeline can see and the point is to make supplying it easy rather than
    correct-format-or-nothing.

      short   two columns, no header, like the barcode files:
                  P01<TAB>GF_M1
              A plate on the left, a label on the right. Duplicates are expected
              -- several plates can belong to the same group.

      long    the same two columns, but keyed by well:
                  P01A01<TAB>GF_M1
              For a plate split across groups. A third column, if present, is
              taken as the broader experimental group.

      full    a header row naming the columns (well, sample, mouse, tissue,
              group), for when mouse and tissue matter separately.

    Plate keys expand to every well the run actually recovered on that plate, so
    the short form does not require knowing which wells came back."""
    if not path or not os.path.exists(path):
        return {}

    # Sniff the delimiter from the first line that carries data rather than from
    # the physical first line. A plate map written by hand often opens with a
    # comment naming the experiment, and a comment with no tab in it made a
    # tab-separated file fall through to the whitespace branch, which splits a
    # label like "GF M1" into two fields and drops the second.
    probe = ""
    with open(path, newline="", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.lstrip().startswith("#"):
                probe = line
                break
    if not probe:
        return {}
    delim = "\t" if "\t" in probe else ("," if "," in probe else None)

    rows = []
    with open(path, newline="", errors="replace") as handle:
        if delim is None:
            # Whitespace-separated, which is how some barcode files are written.
            for line in handle:
                if line.strip() and not line.lstrip().startswith("#"):
                    rows.append(line.split())
        else:
            for row in csv.reader(handle, delimiter=delim):
                if row and any(field.strip() for field in row) and \
                        not row[0].lstrip().startswith("#"):
                    rows.append([field.strip() for field in row])
    if not rows:
        return {}

    header = [field.strip().lower() for field in rows[0]]
    # plate_or_well is here because the frontend template used to ship with that
    # header. It was recognised as a header and then matched no key column, so
    # every row was skipped and the map came back empty with no error. The
    # template now says "well"; this keeps the copies already in people's
    # folders working.
    has_header = any(h in header for h in ("well", "plate", "sample", "bio_id",
                                           "bio.id", "group", "plate_well",
                                           "plate_id", "plate_or_well"))

    mapping = {}
    if has_header:
        index = {name: i for i, name in enumerate(header)}
        def get(row, *names):
            for name in names:
                i = index.get(name)
                if i is not None and i < len(row):
                    return row[i].strip()
            return ""
        for row in rows[1:]:
            if not row:
                continue
            key = get(row, "well", "plate_well", "plate_id", "plate_or_well",
                      "plate")
            if not key:
                continue
            sample = get(row, "sample", "bio_id", "bio.id") or get(row, "group") or ""
            mapping[key] = {
                "sample": sample,
                "mouse": get(row, "mouse"),
                "tissue": get(row, "tissue"),
                "group": get(row, "group", "strain") or sample,
            }
    else:
        for row in rows:
            if len(row) < 2 or not row[0]:
                continue
            key = row[0].strip()
            label = row[1].strip()
            group = row[2].strip() if len(row) > 2 and row[2].strip() else label
            mapping[key] = {"sample": label, "mouse": "", "tissue": "", "group": group}

    # Normalise well keys; leave plate keys alone so they can be expanded later.
    resolved = {}
    for key, value in mapping.items():
        resolved[canonical_plate(key) if PLATE_ONLY.match(key)
                 else canonical_well(key)] = value
    return resolved


def platemap_for_well(platemap, well):
    """The entry for a well, falling back to its plate for the short format."""
    if not platemap:
        return {}
    hit = platemap.get(well)
    if hit is not None:
        return hit
    match = WELL_CANONICAL.match(well) or WELL_NO_PLATE.match(well)
    if match:
        try:
            plate = canonical_plate(match.group("plate"))
        except IndexError:
            return {}
        return platemap.get(plate, {})
    return {}


def write_tsv_and_csv(path_tsv, header, rows):
    """Both formats, every time. A TSV is what the pipeline reads; a CSV is what
    opens in a spreadsheet without an import dialog."""
    os.makedirs(os.path.dirname(path_tsv) or ".", exist_ok=True)
    with open(path_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    path_csv = re.sub(r"\.tsv$", ".csv", path_tsv)
    with open(path_csv, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def validate_against(directory, max_distance):
    """Re-derive clones from published annotations and score against their calls.

    Agreement is measured over pairs: of all pairs of sequences, how many does
    each rule place in the same clone. Precision is the fraction of pairs this
    rule joins that the published calls also join -- the number that matters,
    because a false merge inflates a clone size and clone size is the result."""
    import glob
    files = sorted(glob.glob(os.path.join(directory, "*.tab")))
    if not files:
        print("No *.tab files in {}".format(directory), file=sys.stderr)
        return 1
    total_pub = total_mine = total_both = 0
    print("{:<20} {:>7} {:>7} {:>7} {:>8} {:>8}".format(
        "sample", "seqs", "pub", "mine", "recall", "prec"))
    for path in files:
        with open(path, newline="", errors="replace") as handle:
            rows = [r for r in csv.DictReader(handle, delimiter="\t")
                    if (r.get("JUNCTION") or "").strip() and r.get("CLONE")]
        if len(rows) < 2:
            continue
        records = [{"v_genes": gene_set(r["V_CALL"]), "j_genes": gene_set(r["J_CALL"]),
                    "junction": r["JUNCTION"].strip()} for r in rows]
        mine = assign_clones(records, max_distance)
        pub = [r["CLONE"] for r in rows]
        pairs = lambda c: c * (c - 1) // 2
        same_pub = sum(pairs(v) for v in collections.Counter(pub).values())
        same_mine = sum(pairs(v) for v in collections.Counter(mine).values())
        same_both = sum(pairs(v) for v in collections.Counter(zip(pub, mine)).values())
        total_pub += same_pub; total_mine += same_mine; total_both += same_both
        print("{:<20} {:>7} {:>7} {:>7} {:>8.4f} {:>8.4f}".format(
            os.path.basename(path).replace("_clone-pass.tab", "")[:20], len(rows),
            len(set(pub)), len(set(mine)),
            same_both / same_pub if same_pub else 1.0,
            same_both / same_mine if same_mine else 1.0))
    print("\nTOTAL  recall {:.4f}  precision {:.4f}  over {} pairs the published "
          "calls join".format(total_both / total_pub if total_pub else 1.0,
                              total_both / total_mine if total_mine else 1.0,
                              total_pub))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--airr", help="IgBLAST AIRR table for the delivered sequences")
    ap.add_argument("--sample", help="run name, stripped from headers")
    ap.add_argument("--validate-against", default="",
                    help="directory of published *_clone-pass.tab files; score the "
                         "clone rule against their CLONE column and exit")
    ap.add_argument("--platemap", default="", help="well -> sample map; optional")
    ap.add_argument("--outdir", default="")
    ap.add_argument("--max-distance", type=int, default=4,
                    help="junction Hamming distance within a clone (paper: 4)")
    ap.add_argument("--min-reads", type=int, default=6,
                    help="reads a clone must have in a well for that well to count (paper: 6)")
    ap.add_argument("--min-wells", type=int, default=2,
                    help="wells a clone must appear in to be kept (paper: 2)")
    ap.add_argument("--productive-only", action="store_true")
    # Which reading of Change-O's --act set to use. Strict never merges two
    # sequences sharing no gene; transitive lets a chain of ambiguous calls put
    # them in one group and links on junction distance alone.
    #
    # Strict is the default because it cannot invent a clone, and against the
    # published calls it scores precision 1.0000 / recall 0.9978.
    #
    # On sequences assembled here it does split more, because V calls disagree
    # within an identical junction 77.3% of the time against 53.7% in the
    # published annotations. But the effect is modest: reproducing the 2020
    # library gives 30,213 clones strict and 28,490 transitive, against a
    # published 21,132. The choice is worth exposing; it is not the explanation
    # for the difference, which is still open.
    ap.add_argument("--gene-grouping", default="strict",
                    choices=["strict", "transitive"])
    # Where clone assignment happens, which is not the same question as where
    # clones are counted.
    #
    # The published pipeline ran Change-O once over every sample together: the
    # released tables carry only 4,736 distinct CLONE ids for 586,821 sequences,
    # and 97.6% of them appear in more than one sample -- one in all forty. The
    # notebook then scopes them per sample with CLONE_ID = BIO.ID + "_" + CLONE.
    #
    # So a clone is defined globally and only afterwards attributed to a sample.
    # Clustering inside each sample instead is a stricter operation, because a
    # single-linkage chain cannot pass through a sequence in another sample, and
    # it yields finer clones. 'global' reproduces the published method;
    # 'sample' keeps clone assignment independent between samples.
    ap.add_argument("--clone-scope", default="sample",
                    choices=["sample", "global"])
    args = ap.parse_args()

    if args.validate_against:
        return validate_against(args.validate_against, args.max_distance)
    if not args.airr or not args.sample or not args.outdir:
        ap.error("--airr, --sample and --outdir are required unless --validate-against is given")

    records, skipped = load_airr(args.airr, args.sample)
    if not records:
        print("No annotated sequence carried a junction and a parsable well. "
              "Nothing to assign clones to.", file=sys.stderr)
        return 1

    if args.productive_only:
        before = len(records)
        records = [r for r in records if r["productive"]]
        skipped["non_productive"] += before - len(records)

    platemap = load_platemap(args.platemap)
    # A well the plate map does not mention is not part of any sample. Naming it
    # after the run made it look like a real sample in the report -- five stray
    # wells of the 2020 library turned up as a sample called MIBS2020 sitting
    # beside the forty genuine ones. Label them for what they are instead.
    UNMAPPED = "(unmapped wells)"
    unmapped_wells = set()
    for rec in records:
        entry = platemap_for_well(platemap, rec["well"])
        if platemap and not entry:
            unmapped_wells.add(rec["well"])
            rec["bio_id"] = UNMAPPED
            rec["mouse"] = ""
            rec["tissue"] = ""
            # Named as a sample, but given no group. The group panels drop rows
            # with an empty group, which is where a diagnostic bucket belongs:
            # labelling it put "(unmapped wells)" on the x axis of "clone size
            # by group" as though five stray wells were an experimental arm.
            rec["group"] = ""
        else:
            rec["bio_id"] = entry.get("sample") or args.sample
            rec["mouse"] = entry.get("mouse", "")
            rec["tissue"] = entry.get("tissue", "")
            rec["group"] = entry.get("group", "")

    # Clones are assigned within a sample. Two mice sharing a junction are a
    # public clonotype, which is a separate question answered further down --
    # never one clone.
    by_sample = collections.defaultdict(list)
    for rec in records:
        by_sample[rec["bio_id"]].append(rec)

    assigner = (assign_clones_transitive if args.gene_grouping == "transitive"
                else assign_clones)
    if args.clone_scope == "global":
        # One pass over everything, then attribute each sequence's global clone
        # to its own sample -- exactly CLONE_ID = BIO.ID + "_" + CLONE.
        labels = assigner(records, args.max_distance)
        for rec, label in zip(records, labels):
            rec["clone"] = label
            rec["clone_id"] = "{}_{}".format(rec["bio_id"], label)
    else:
        for bio_id, group in sorted(by_sample.items()):
            labels = assigner(group, args.max_distance)
            for rec, label in zip(group, labels):
                rec["clone"] = label
                rec["clone_id"] = "{}_{}".format(bio_id, label)

    # Collapse a clone within a well and sum its reads: the published filter
    # operates on reads per clone per well, not per sequence.
    well_reads = collections.Counter()
    well_repr = {}
    for rec in records:
        key = (rec["clone_id"], rec["well"])
        well_reads[key] += rec["reads"] or 0
        current = well_repr.get(key)
        if current is None or (rec["reads"] or 0) > (current["reads"] or 0):
            well_repr[key] = rec

    passing = {k: v for k, v in well_reads.items() if v >= args.min_reads}
    wells_per_clone = collections.Counter(clone for clone, _ in passing)
    kept_clones = {c for c, n in wells_per_clone.items() if n >= args.min_wells}

    # ---- clone table ------------------------------------------------------
    clone_rows = []
    for clone_id in sorted(wells_per_clone, key=lambda c: (-wells_per_clone[c], c)):
        members = [well_repr[k] for k in passing if k[0] == clone_id]
        if not members:
            continue
        best = max(members, key=lambda r: r["reads"] or 0)
        wells = wells_per_clone[clone_id]
        reads = sum(passing[k] for k in passing if k[0] == clone_id)
        clone_rows.append([
            clone_id, best["bio_id"], best["mouse"], best["tissue"], best["group"],
            wells, reads, wells * reads,
            best["v_call"], best["d_call"], best["j_call"],
            best["junction"], len(best["junction"]), best["junction_aa"],
            "kept" if clone_id in kept_clones else "single_well",
        ])
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.clones.tsv".format(args.sample)),
        ["clone_id", "sample", "mouse", "tissue", "group", "wells", "reads",
         "clone_by_well", "v_call", "d_call", "j_call", "junction",
         "junction_length", "junction_aa", "verdict"],
        clone_rows)

    # ---- clone x well table ----------------------------------------------
    incidence_rows = []
    for (clone_id, well), reads in sorted(well_reads.items()):
        rec = well_repr[(clone_id, well)]
        incidence_rows.append([
            clone_id, well, rec["bio_id"], reads,
            "pass" if reads >= args.min_reads else "below_min_reads",
            "kept" if clone_id in kept_clones else "single_well",
            rec["sequence_id"],
        ])
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.clone_by_well.tsv".format(args.sample)),
        ["clone_id", "well", "sample", "reads", "well_verdict", "clone_verdict",
         "sequence_id"],
        incidence_rows)

    # ---- public clonotypes ------------------------------------------------
    # An exact junction match across mice. Every sample in a MIBS library is
    # sequenced together, so this needs no machinery outside the run -- the
    # thing the plan expected to be an architectural addition is already here.
    by_junction = collections.defaultdict(set)
    junction_clones = collections.defaultdict(set)
    for clone_id in kept_clones:
        members = [well_repr[k] for k in passing if k[0] == clone_id]
        if not members:
            continue
        rec = members[0]
        identity = rec["mouse"] or rec["bio_id"]
        by_junction[rec["junction_aa"] or rec["junction"]].add(identity)
        junction_clones[rec["junction_aa"] or rec["junction"]].add(clone_id)

    public_rows = []
    # The same owner sets, unjoined, because the per-sample table below has to
    # ask which of them a sample is in and splitting the ";"-joined column back
    # apart would break on any label containing a semicolon.
    public_owners = []
    for junction, owners in sorted(by_junction.items(), key=lambda kv: -len(kv[1])):
        if len(owners) < 2:
            continue
        public_rows.append([junction, len(owners), ";".join(sorted(owners)),
                            ";".join(sorted(junction_clones[junction]))])
        public_owners.append(owners)
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.public_clonotypes.tsv".format(args.sample)),
        ["junction_aa", "individuals", "shared_by", "clone_ids"],
        public_rows)

    # ---- per-sample statistics --------------------------------------------
    # Everything above is a number for the run, and a MIBS run is several mice.
    # Incidence only means something within a sample, so a run-level clone count
    # is an average over designs that may differ -- one mouse with sixteen wells
    # and another with thirty-two do not contribute comparably to it. This is
    # the same summary computed per sample, which is the thing a plate map is
    # attached in order to get. Without one there is a single row here, and that
    # is the honest answer rather than a missing file.
    sample_clones = collections.defaultdict(set)
    for clone_id, well in passing:
        sample_clones[well_repr[(clone_id, well)]["bio_id"]].add(clone_id)

    def joined(values):
        return ";".join(sorted({v for v in values if v}))

    stats_rows = []
    for bio_id in sorted(by_sample):
        group = by_sample[bio_id]
        mice = sorted({rec["mouse"] for rec in group if rec["mouse"]})
        # The public table keys on the mouse when a record names one and falls
        # back to the sample, so the count below has to ask the same question
        # the same way. A sample naming several mice is its own identity there.
        identity = mice[0] if len(mice) == 1 else bio_id
        called = sample_clones.get(bio_id, set())
        kept = called & kept_clones
        incidence = [wells_per_clone[c] for c in kept]
        stats_rows.append([
            bio_id,
            joined(rec["group"] for rec in group),
            joined(mice),
            joined(rec["tissue"] for rec in group),
            len({rec["well"] for rec in group}),
            len(group),
            sum(rec["reads"] or 0 for rec in group),
            len(called),
            len(kept),
            # Everything that failed either filter, which is not the same thing
            # as a single-well clone once --min-wells is settable from the form:
            # at --min-wells 5 this counts clones seen in four wells. Named for
            # the filter rather than for the default.
            len(called) - len(kept),
            max(incidence) if incidence else 0,
            round(sum(incidence) / len(incidence), 3) if incidence else 0,
            sum(1 for owners in public_owners if identity in owners),
        ])
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.sample_stats.tsv".format(args.sample)),
        ["sample", "group", "mouse", "tissue", "wells", "sequences", "reads",
         "clones_called", "clones_kept", "clones_below_min_wells",
         "max_wells_per_clone", "mean_wells_per_clone", "public_clonotypes"],
        stats_rows)

    # ---- per-sequence table, in the schema the lab's notebooks read --------
    changeo_rows = []
    for rec in sorted(records, key=lambda r: (r["bio_id"], r["well"], -(r["reads"] or 0))):
        plate = rec["well"][:3] if len(rec["well"]) > 3 else ""
        changeo_rows.append([
            rec["sequence_id"], rec["sequence"], "T" if rec["productive"] else "F",
            rec["v_call"], rec["d_call"], rec["j_call"], rec["junction"],
            len(rec["junction"]), rec["junction_aa"], rec["bio_id"], plate,
            rec["well"][3:] if len(rec["well"]) > 3 else rec["well"], rec["well"],
            rec["reads"], rec["group"], rec["mouse"], rec["tissue"],
            rec.get("clone", ""), rec["v_identity"],
        ])
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.sequences.tsv".format(args.sample)),
        ["SEQUENCE_ID", "SEQUENCE_INPUT", "FUNCTIONAL", "V_CALL", "D_CALL",
         "J_CALL", "JUNCTION", "JUNCTION_LENGTH", "JUNCTION_AA", "BIO.ID",
         "PLATE.NUMBER", "PLATE.WELL", "PLATE_ID", "READ.COUNTS", "STRAIN",
         "MOUSE", "TISSUE", "CLONE", "V_IDENTITY"],
        changeo_rows)

    # ---- the design as the run resolved it --------------------------------
    # A finished run's only record of which design it used is "platemap_used
    # yes|no" in the summary and a line on stderr, so two runs of one library
    # under two different designs are indistinguishable from their outputs.
    # This writes the resolution itself: every well that came back, the sample
    # it was attributed to, and whether the plate map named it or it fell
    # through to the unmapped bucket. Written on every run, mapped or not,
    # because Snakemake fails a rule that does not produce a declared output.
    resolved_rows = []
    seen_wells = set()
    for rec in sorted(records, key=lambda r: r["well"]):
        if rec["well"] in seen_wells:
            continue
        seen_wells.add(rec["well"])
        resolved_rows.append([
            rec["well"], rec["bio_id"], rec["mouse"], rec["tissue"], rec["group"],
            "yes" if platemap and rec["well"] not in unmapped_wells else "no",
        ])
    write_tsv_and_csv(
        os.path.join(args.outdir, "{}.platemap_resolved.tsv".format(args.sample)),
        ["well", "sample", "mouse", "tissue", "group", "in_platemap"],
        resolved_rows)

    # ---- plate map template ----------------------------------------------
    # TBRC is handed a folder of FASTQs and a barcode set. It has no way to know
    # there was an experimental design, which mouse a well came from, or what
    # tissue it was. The paper's private-repertoire analysis needs exactly that,
    # and without it every well collapses into one sample and the comparison
    # disappears.
    #
    # The pipeline cannot invent the design, but it does know which wells exist.
    # So it writes the form: every well actually recovered, in the order they
    # were found, with the columns to fill in. Filling it in and re-running the
    # clone step is cheap -- annotation, the expensive part, does not repeat.
    template = os.path.join(args.outdir, "{}.platemap_template.tsv".format(args.sample))
    if not platemap:
        wells_sorted = sorted({r["well"] for r in records})
        rows = [[w, args.sample, "", "", ""] for w in wells_sorted]
        write_tsv_and_csv(template, ["well", "sample", "mouse", "tissue", "group"], rows)

    # ---- summary ----------------------------------------------------------
    wells_seen = {r["well"] for r in records}
    samples_seen = {r["bio_id"] for r in records}
    kept_wells = [wells_per_clone[c] for c in kept_clones]
    summary = [
        ["sequences_annotated", len(records)],
        ["wells", len(wells_seen)],
        ["samples", len(samples_seen)],
        # The number the report should show. '(unmapped wells)' is a diagnostic
        # bucket, not a mouse, and counting it as a sample over-reports the run
        # by one whenever a single well falls outside the design.
        ["samples_mapped", len(samples_seen - {UNMAPPED})],
        # Over the mapped design only. Counting the unmapped bucket as a sample
        # and its wells as that sample's put a made-up mouse in the denominator
        # and its strays in the numerator, so the average moved in both
        # directions at once.
        ["wells_per_sample",
         round(len(wells_seen - unmapped_wells) /
               max(1, len(samples_seen - {UNMAPPED})), 2)],
        ["clones_called", len(wells_per_clone)],
        ["clones_kept", len(kept_clones)],
        # Historical name, kept because scripts read this file: it counts every
        # clone that failed either filter, which is only "single well" while
        # min_wells is 2.
        ["clones_single_well", len(wells_per_clone) - len(kept_clones)],
        ["max_wells_per_clone", max(kept_wells) if kept_wells else 0],
        ["mean_wells_per_clone", round(sum(kept_wells) / len(kept_wells), 3) if kept_wells else 0],
        ["public_clonotypes", len(public_rows)],
        ["min_reads_per_well", args.min_reads],
        ["min_wells_per_clone", args.min_wells],
        ["junction_max_distance", args.max_distance],
        ["gene_grouping", args.gene_grouping],
        ["clone_scope", args.clone_scope],
        ["platemap_used", "yes" if platemap else "no"],
        ["wells_not_in_platemap", len(unmapped_wells)],
    ]
    for reason, count in sorted(skipped.items()):
        summary.append(["skipped_" + reason, count])
    write_tsv_and_csv(os.path.join(args.outdir, "{}.summary.tsv".format(args.sample)),
                      ["metric", "value"], summary)

    print("{} sequences -> {} clones in {} samples; {} kept (>= {} reads in a well, "
          ">= {} wells), {} dropped as single-well, {} public clonotypes.".format(
              len(records), len(wells_per_clone), len(samples_seen), len(kept_clones),
              args.min_reads, args.min_wells,
              len(wells_per_clone) - len(kept_clones), len(public_rows)))
    if platemap and unmapped_wells:
        print("{} well(s) were recovered that the plate map does not list; they are "
              "reported as '{}' rather than mixed into a sample.".format(
                  len(unmapped_wells), UNMAPPED), file=sys.stderr)
    if not platemap:
        print("No plate map was supplied, so all {} wells were treated as one "
              "sample and no per-sample or per-group comparison is possible. "
              "A template listing every well found has been written to "
              "{}.platemap_template.tsv -- fill in sample, mouse, tissue and "
              "group, save it as platemap.tsv beside the FASTQ files, and rerun. "
              "Clone incidence is only meaningful within a sample.".format(
                  len(wells_seen), args.sample), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
