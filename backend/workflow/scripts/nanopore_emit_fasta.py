#!/usr/bin/env python3
"""Rewrite per-well consensus headers into the shared final.fasta convention.

The QC scripts every pipeline shares parse plate and well out of the FASTA
header with a fixed pattern:

    ^(.*?)(A01P01|P01A01)_(\\d+)-(\\d+)$

so a nanopore header like

    P01A09|reads=16|unique=16|domfrac=0.889|method=poa_racon

stops plate_depth_qc.R with "Failed to parse plate metadata". The information is
all there, just not in the shape the shared code reads.

Rather than loosen a parser that four other scripts depend on, the shared
final.fasta is written in the shared convention and the rich header is kept on
the nanopore-specific consensus file, which is what the archive delivers and
what carries the dominance and method fields.

Well-only libraries have no plate. They get P00 so the identifier still parses,
and the plate map then shows a single plate, which is what a one-plate-per-file
library actually is.
"""
import re, sys

src, dst = sys.argv[1], sys.argv[2]

PLATEWELL = re.compile(r"^(P\d{2})([A-H]\d{2})$")
WELLONLY = re.compile(r"^([A-H]\d{2})$")


def kv(header):
    out = {}
    for part in header.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


written = skipped = 0
with open(src) as fh, open(dst, "w") as out:
    name, buf = None, []
    def flush():
        global written, skipped
        if name is None:
            return
        ident = name.split("|")[0]
        meta = kv(name)
        m = PLATEWELL.match(ident)
        if m:
            tag = ident
        else:
            w = WELLONLY.match(ident)
            if not w:
                skipped += 1
                return
            # No plate cassette in this library; P00 keeps the identifier
            # parseable and the plate map honest about there being one plate.
            tag = "P00" + ident
        reads = meta.get("reads", "0")
        uniq = meta.get("unique", reads)
        out.write(">%s_%s-%s\n%s\n" % (tag, reads, uniq, "".join(buf)))
        written += 1
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            name, buf = line[1:], []
        elif line:
            buf.append(line.strip())
    flush()

sys.stderr.write("[emit] %d sequences written, %d headers unparseable\n"
                 % (written, skipped))
if skipped:
    sys.stderr.write("[emit] unparseable identifiers are a bug, not a data "
                     "property: every consensus is keyed by plate+well or well\n")
