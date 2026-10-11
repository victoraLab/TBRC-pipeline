#!/bin/bash
# Regenerate docs/changelog.tsv from the git history.
#
# git is the only honest source for "what changed and when" -- a hand-kept file
# drifts the moment someone forgets to update it. Each commit becomes one row:
# its date, its subject, the finding IDs it mentions, and its first paragraph.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
OUT="${1:-${ROOT}/docs/changelog.tsv}"

cd "${ROOT}"
{
    printf 'date\tcommit\tkind\tfindings\ttitle\tdetail\n'
    git log --reverse --date=short --pretty=format:'%ad%x1f%h%x1f%s%x1f%b%x1e' | \
    python3 -c '
import sys, re
raw = sys.stdin.read()
for rec in raw.split("\x1e"):
    rec = rec.strip("\n")
    if not rec.strip():
        continue
    parts = rec.split("\x1f")
    if len(parts) < 4:
        continue
    date, sha, subject, body = parts[0], parts[1], parts[2], parts[3]

    findings = ",".join(sorted(set(re.findall(r"\bF\d{1,2}\b", subject + " " + body)),
                               key=lambda f: int(f[1:])))

    low = subject.lower()
    FIX = ("fix", "stop", "correct", "close", "keep", "restore", "cover",
           "parse", "make ", "sync", "prevent", "repair", "unblock",
           "fill", "validate")
    FEATURE = ("add", "give", "show", "propose", "offer", "introduce", "expose")
    HOUSE = ("document", "rename", "continue", "derive", "merge", "initial",
             "improve deployment", "unify")

    if any(low.startswith(w) for w in HOUSE):
        kind = "housekeeping"
    elif findings:
        # A commit that cites an audit finding is fixing it, whatever verb it opens with.
        kind = "fix"
    elif any(low.startswith(w) for w in FIX):
        kind = "fix"
    elif any(low.startswith(w) for w in FEATURE):
        kind = "feature"
    else:
        kind = "change"

    # first paragraph of the body, collapsed to one line
    para = ""
    for block in body.split("\n\n"):
        block = block.strip()
        if block and not block.startswith("Co-Authored-By"):
            para = " ".join(block.split())
            break
    para = re.sub(r"\s*Co-Authored-By:.*$", "", para).strip()

    for field in (subject, para):
        if "\t" in field:
            field = field.replace("\t", " ")
    print("\t".join([date, sha, kind, findings,
                     subject.replace("\t", " "), para.replace("\t", " ")]))
'
} > "${OUT}"

echo "wrote ${OUT}: $(( $(wc -l < "${OUT}") - 1 )) entries"
