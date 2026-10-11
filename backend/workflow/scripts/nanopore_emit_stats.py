#!/usr/bin/env python3
"""Write collapsed.stats.txt from the per-well consensus table.

plate_stats.R reads this file to draw the plate maps, and it expects one line
per well of "<well identifier><tab><count>". For a nanopore run there is exactly
one consensus per well, so the count that belongs here is the read support
behind that consensus -- which is what the plate map should be showing.

This is a separate script rather than an inline heredoc because Snakemake
formats shell blocks, so any literal brace in embedded Python is read as a
Snakemake placeholder and the rule dies with "tuple index out of range".
"""
import csv, sys

src, dst = sys.argv[1], sys.argv[2]
written = 0
with open(src) as fh, open(dst, "w") as out:
    for row in csv.DictReader(fh, delimiter="\t"):
        length = row.get("consensus_length", "0")
        if length not in ("", "0"):
            out.write("%s\t%s\n" % (row["platewell"], row.get("reads_used", 0)))
            written += 1
sys.stderr.write("[emit] %d wells written to %s\n" % (written, dst))
