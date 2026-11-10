#!/usr/bin/env python3
"""Compare two folders of Snakemake dry runs (written by dryrun_capture.sh) command by command.

    python3 backend/tests/compare_dryrun.py BEFORE_DIR AFTER_DIR

Quoting every value in a rule's shell command changes how the command is written, not what it does.
This is the check that says so. Each job's command is split into words the way a shell does
(shlex), so 'x', "x" and x are the same word, and the words of the two runs are compared in order.
An empty word counts: an argument that quoting loses (-t '' becoming -t) is a difference.

A job's inputs, outputs, log files, threads and wildcards are compared as printed too, so a rule
that starts reading or writing a different file is caught as well.

Exit status 0: every pipeline matches. 1: something differs (each difference is printed, with the
exact words on both sides). 2: a folder could not be read.

Plain standard library; runs on Python 3.6 and later.
"""
import difflib
import os
import re
import shlex
import sys

# The line Snakemake prints before each job, e.g. "[Mon Nov  9 04:53:39 2026]".
_STAMP = re.compile(r"^\[[A-Z][a-z]{2} [A-Z][a-z]{2} +\d{1,2} \d\d:\d\d:\d\d \d{4}\]$")
_RULE = re.compile(r"^(?:local)?rule (\S+):$")
# Printed under the rule line, always ending in "resources:"; what follows is the command.
_HEADER = re.compile(r"^    ([a-z_]+): ?(.*)$")
_COMPARED_HEADERS = ("input", "output", "log", "threads", "wildcards", "temporary")
# Snakemake's own note after a job whose output is temp(). It is not part of the command, and it
# lists the files in an order that changes from one run to the next, so the files are compared as a
# sorted list.
_TEMP_NOTE = "Would remove temporary output "


class Unreadable(Exception):
    pass


def parse_dryrun(text, label):
    """[(rule, {header: text}, [words])] in print order."""
    lines = text.splitlines()
    jobs = []
    index = 0
    while index < len(lines):
        match = _RULE.match(lines[index])
        if not match:
            index += 1
            continue
        rule = match.group(1)
        index += 1
        headers = {}
        while index < len(lines):
            header = _HEADER.match(lines[index])
            if not header:
                break
            headers[header.group(1)] = header.group(2)
            index += 1
            if header.group(1) == "resources":
                break
        if "resources" not in headers:
            raise Unreadable("%s: rule %s has no resources line" % (label, rule))
        body = []
        temporary = []
        while index < len(lines) and not _STAMP.match(lines[index]) \
                and not lines[index].startswith("Job stats:"):
            if lines[index].startswith(_TEMP_NOTE):
                temporary.append(lines[index][len(_TEMP_NOTE):])
            else:
                body.append(lines[index])
            index += 1
        if temporary:
            headers["temporary"] = ", ".join(sorted(temporary))
        command = "\n".join(body).strip()
        try:
            words = shlex.split(command, comments=True)
        except ValueError as err:
            raise Unreadable("%s: the command of rule %s cannot be split into words (%s)"
                             % (label, rule, err))
        jobs.append((rule, headers, words))
    return jobs


def describe_words(before, after):
    """The places two word lists differ, each with two words of context."""
    out = []
    matcher = difflib.SequenceMatcher(None, before, after, autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        out.append("      at word %d (after %r): before %r, after %r"
                   % (i1, before[max(0, i1 - 2):i1], before[i1:i2], after[j1:j2]))
    return out


def compare_pipeline(name, before_jobs, after_jobs):
    """The lines describing every difference between two parsed dry runs (none when they match)."""
    problems = []
    before_by_rule = {}
    after_by_rule = {}
    for rule, headers, words in before_jobs:
        before_by_rule.setdefault(rule, []).append((headers, words))
    for rule, headers, words in after_jobs:
        after_by_rule.setdefault(rule, []).append((headers, words))
    for rule in sorted(set(before_by_rule) | set(after_by_rule)):
        old = before_by_rule.get(rule, [])
        new = after_by_rule.get(rule, [])
        if len(old) != len(new):
            problems.append("  %s: rule %s runs %d time(s) before and %d after"
                            % (name, rule, len(old), len(new)))
            continue
        for number, ((old_headers, old_words), (new_headers, new_words)) in enumerate(zip(old, new)):
            where = "%s: rule %s%s" % (name, rule, "" if len(old) == 1 else " #%d" % (number + 1))
            for key in _COMPARED_HEADERS:
                if old_headers.get(key) != new_headers.get(key):
                    problems.append("  %s: %s differs: before %r, after %r"
                                    % (where, key, old_headers.get(key), new_headers.get(key)))
            if old_words != new_words:
                problems.append("  %s: command differs (%d words before, %d after)"
                                % (where, len(old_words), len(new_words)))
                problems.extend(describe_words(old_words, new_words))
    return problems


def read_folder(folder):
    """{pipeline: text} for every .txt file in the folder."""
    if not os.path.isdir(folder):
        raise Unreadable("%s is not a folder" % folder)
    found = {}
    for entry in sorted(os.listdir(folder)):
        if entry.endswith(".txt"):
            with open(os.path.join(folder, entry), encoding="utf-8", errors="replace") as handle:
                found[entry[:-4]] = handle.read()
    if not found:
        raise Unreadable("%s holds no .txt dry runs" % folder)
    return found


def compare_folders(before_dir, after_dir):
    """(report lines, number of pipelines that differ)."""
    before = read_folder(before_dir)
    after = read_folder(after_dir)
    report = []
    differing = 0
    for name in sorted(set(before) | set(after)):
        if name not in before or name not in after:
            report.append("%s: only in %s" % (name, "AFTER" if name in after else "BEFORE"))
            differing += 1
            continue
        old_jobs = parse_dryrun(before[name], "BEFORE/" + name)
        new_jobs = parse_dryrun(after[name], "AFTER/" + name)
        if not old_jobs:
            raise Unreadable("BEFORE/%s holds no jobs; was the dry run captured?" % name)
        problems = compare_pipeline(name, old_jobs, new_jobs)
        if problems:
            differing += 1
            report.append("%s: DIFFERENT" % name)
            report.extend(problems)
        else:
            commands = sum(1 for _, _, words in old_jobs if words)
            report.append("%s: identical (%d jobs, %d commands, %d words)"
                          % (name, len(old_jobs), commands, sum(len(w) for _, _, w in old_jobs)))
    return report, differing


def main(argv):
    if len(argv) != 3:
        print("Usage: compare_dryrun.py BEFORE_DIR AFTER_DIR", file=sys.stderr)
        return 2
    try:
        report, differing = compare_folders(argv[1], argv[2])
    except Unreadable as err:
        print("compare_dryrun: %s" % err, file=sys.stderr)
        return 2
    print("\n".join(report))
    if differing:
        print("\n%d pipeline(s) differ." % differing)
        return 1
    print("\nNo differences.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
