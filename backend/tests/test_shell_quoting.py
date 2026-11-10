"""Tests that every value in the pipeline's shell commands is quoted, and that the commands
that go over ssh quote the paths they carry.

Snakemake pastes a rule's values into its shell string before bash sees it, so a value with a quote,
a space or a $( ) in it changes the command unless it is quoted. `{params.x:q}` quotes it; a quote
written by hand around `{params.x}` does not, since the value can contain that quote. These tests
scan every `shell:` string in backend/workflow/rules (and FINAL_SHELL, which final.smk shares between
two rules) and refuse a placeholder that

  * has no :q, unless it is one of the few listed in ALLOWED_UNQUOTED with the reason, or
  * sits inside quotes the string writes by hand ('{params.x}', "{params.x}", "a/{params.x}/b"),
    or inside a comment, where a newline in the value would end the comment and start a command.

Quoting is not the same as keeping the commands as they were. That is checked by
dryrun_capture.sh + compare_dryrun.py, which run every pipeline's dry run and compare the words of
each command before and after; the tests at the bottom cover those two tools and compare_archives.py.

Run from the repository root, in an image with bash >= 4:
    python3 -m unittest backend/tests/test_shell_quoting.py -v
"""
import ast
import contextlib
import io
import json
import os
import re
import shlex
import shutil
import string
import subprocess
import sys
import tarfile
import tempfile
import tokenize
import unittest
import zipfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TESTS = os.path.join(REPO, "backend", "tests")
RULES = os.path.join(REPO, "backend", "workflow", "rules")
sys.path.insert(0, TESTS)

import compare_archives  # noqa: E402
import compare_dryrun  # noqa: E402

# Placeholders that may stand without :q. Everything else that names a value of the run needs it.
ALLOWED_UNQUOTED = {
    "threads": "an integer Snakemake works out from the rule's threads: line",
    "params.rev_primers": "several words on purpose (--inner-reverse-primer SEQ, repeated); every "
                          "primer is checked against [ACGTacgt]{15,40} when the Snakefile loads",
    "_NP_PATH": "a constant: the PATH prefix built from pipelines.yaml, not from the submitted run",
    "_SCRIPTS": "a constant: the scripts folder, 'workflow/scripts'",
}

# The first part of a placeholder's name that marks a value the rule was given.
VALUE_NAMES = ("params", "input", "output", "log", "wildcards", "resources")

# The fewest `shell:` directives there are today. Fewer found means the scan above missed some.
SHELL_DIRECTIVES = 23


# ------------------------------------------------------------------ finding the shell strings
def _literal_after(source, pos):
    """What follows `shell:` or `NAME =` at source[pos:]: ('string', value, line) for a string
    literal (line is its first line, 1-based, in source), ('name', NAME, None) for a name standing
    in for one. Anything else (a call, a concatenation) cannot be scanned, so it raises."""
    base = source.count("\n", 0, pos)
    tokens = tokenize.generate_tokens(io.StringIO(source[pos:]).readline)
    skip = (tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.COMMENT)
    for token in tokens:
        if token.type in skip:
            continue
        following = next(tokens)
        if following.type != tokenize.NEWLINE:
            raise AssertionError("the shell string at line %d is not a plain string: it is followed "
                                 "by %r" % (base + token.start[0], following.string))
        if token.type == tokenize.STRING:
            return "string", ast.literal_eval(token.string), base + token.start[0]
        if token.type == tokenize.NAME:
            return "name", token.string, None
        break
    raise AssertionError("no string after the shell directive at line %d" % (base + 1))


def shell_blocks():
    """[(label, text, first line)] for every shell string in the rule files, and the number of
    `shell:` directives met (final.smk's two refer to FINAL_SHELL, which is listed once)."""
    blocks = []
    directives = 0
    for name in sorted(os.listdir(RULES)):
        if not name.endswith(".smk"):
            continue
        with open(os.path.join(RULES, name), encoding="utf-8") as handle:
            source = handle.read()
        for match in re.finditer(r"^[ \t]*shell:[ \t]*", source, re.M):
            directives += 1
            kind, value, line = _literal_after(source, match.end())
            rule = re.findall(r"^rule (\w+):", source[:match.start()], re.M)[-1]
            if kind == "string":
                blocks.append(("%s rule %s" % (name, rule), value, line))
        for match in re.finditer(r"^(\w*SHELL\w*)[ \t]*=[ \t]*", source, re.M):
            kind, value, line = _literal_after(source, match.end())
            if kind == "string":
                blocks.append(("%s %s" % (name, match.group(1)), value, line))
    return blocks, directives


# ------------------------------------------------------------------ where a placeholder sits
class ShellPosition:
    """Follows a shell string character by character and says where the next placeholder lands:
    'plain' (a word of a command), 'single' / 'double' (inside quotes the string wrote), or
    'comment'. $( ) and ( ) start a fresh plain context, even inside double quotes, as in bash.
    It reads the way these rules are written, not all of shell: no here-documents, no backticks."""

    def __init__(self):
        self.stack = ["plain"]
        self.escaped = False
        self.comment = False
        self.previous = "\n"

    def feed(self, text):
        for char in text:
            if self.comment:
                self.comment = char != "\n"
            elif self.stack[-1] == "single":
                if char == "'":
                    self.stack.pop()
            elif self.escaped:
                self.escaped = False
            elif char == "\\":
                self.escaped = True
            elif self.stack[-1] == "double":
                if char == '"':
                    self.stack.pop()
                elif char == "(" and self.previous == "$":
                    self.stack.append("plain")
            elif char == "'":
                self.stack.append("single")
            elif char == '"':
                self.stack.append("double")
            elif char == "(":
                self.stack.append("plain")
            elif char == ")":
                if len(self.stack) > 1:
                    self.stack.pop()
            elif char == "#" and (self.previous.isspace() or self.previous in ";&|("):
                self.comment = True
            self.previous = char

    def where(self):
        return "comment" if self.comment else self.stack[-1]

    def balanced(self):
        return self.stack == ["plain"] and not self.escaped


def check_shell_string(text):
    """The problems with one shell string, as (line in the string, what is wrong) pairs."""
    problems = []
    position = ShellPosition()
    line = 1
    for literal, field, spec, _ in string.Formatter().parse(text):
        position.feed(literal)
        line += literal.count("\n")
        if field is None:
            continue
        if field in ALLOWED_UNQUOTED:
            continue
        if re.split(r"[.\[]", field, 1)[0] not in VALUE_NAMES:
            problems.append((line, "{%s} is not a value of the rule and is not in ALLOWED_UNQUOTED" % field))
            continue
        if spec != "q":
            problems.append((line, "{%s} has no :q" % field))
        if position.where() != "plain":
            problems.append((line, "{%s%s} sits inside %s quotes or a comment written by hand; put "
                             "it outside them" % (field, ":" + spec if spec else "", position.where())))
    if not position.balanced():
        problems.append((line, "the quotes or brackets of this string do not balance, so the scan "
                         "cannot be trusted"))
    return problems


class TestScanner(unittest.TestCase):
    """The scan has to be right before its verdict on the rules means anything."""

    def problems(self, text):
        return [what for _, what in check_shell_string(text)]

    def test_a_bare_placeholder_is_refused(self):
        self.assertEqual(self.problems("echo {params.x}"), ["{params.x} has no :q"])
        self.assertEqual(self.problems("cat {input[0]} > {output}"),
                         ["{input[0]} has no :q", "{output} has no :q"])
        self.assertEqual(self.problems("echo {log.std}"), ["{log.std} has no :q"])

    def test_quotes_written_by_hand_are_refused(self):
        self.assertEqual(len(self.problems("echo '{params.x:q}'")), 1)
        self.assertEqual(len(self.problems('echo "{params.x:q}"')), 1)
        self.assertEqual(len(self.problems("echo '{params.x}'")), 2)
        self.assertEqual(len(self.problems('W="{params.a:q}/work/{params.b:q}"')), 2)
        self.assertEqual(len(self.problems('echo "$W/m_{params.b:q}_D"')), 1)

    def test_a_placeholder_in_a_comment_is_refused(self):
        self.assertEqual(len(self.problems("# the value is {params.x:q}\necho hi")), 1)

    def test_quoted_placeholders_outside_any_quotes_pass(self):
        for text in (
            "echo {params.x:q}",
            "W={params.a:q}/work/mouse/{params.b:q}",
            'f "$W/mouse_"{params.b:q}"_D.nhr"',
            'mkdir -p "$(dirname {output.report:q})"',
            "cp {params.prefix:q}.tsv {output.table:q}",
            "n=$(grep -c '^>' {output.final:q} || true)",
            "echo {threads} {_SCRIPTS}/x.py {params.rev_primers}",
            "x=$(( $y + 1 )); echo {params.x:q}",
            "ls {params.dir:q}*.fastq | xargs -I{{}} echo {{}}",
        ):
            self.assertEqual(self.problems(text), [], text)

    def test_an_apostrophe_in_a_comment_does_not_confuse_the_scan(self):
        self.assertEqual(self.problems("# the 5' end\necho {params.x:q}\n"), [])
        self.assertEqual(self.problems("echo hi # don't\necho {params.x:q}\n"), [])

    def test_a_word_that_is_not_a_value_of_the_rule_is_refused(self):
        self.assertEqual(len(self.problems("echo {something_else}")), 1)

    def test_doubled_braces_are_shell_not_placeholders(self):
        self.assertEqual(self.problems('f() {{ echo "${{1:-x}}"; }}\nf {params.x:q}'), [])

    def test_a_string_that_does_not_balance_is_reported(self):
        self.assertEqual(len(self.problems("echo 'unfinished {params.x:q}")), 2)


# ------------------------------------------------------------------ the rules themselves
class TestRuleQuoting(unittest.TestCase):
    def test_the_scan_sees_every_shell_string(self):
        blocks, directives = shell_blocks()
        self.assertGreaterEqual(directives, SHELL_DIRECTIVES)
        self.assertTrue(any("FINAL_SHELL" in label for label, _, _ in blocks))
        self.assertGreaterEqual(len(blocks), SHELL_DIRECTIVES - 2)

    def test_no_unquoted_value_in_shell_strings(self):
        blocks, _ = shell_blocks()
        found = []
        for label, text, first in blocks:
            for line, what in check_shell_string(text):
                shown = text.splitlines()[line - 1].strip() if line <= len(text.splitlines()) else ""
                found.append("%s, line %d: %s\n      %s" % (label, first + line - 1, what, shown[:130]))
        self.assertEqual(found, [], "\n%d placeholder(s) not safely quoted:\n%s"
                         % (len(found), "\n".join(found)))

    def test_every_allowed_placeholder_still_has_a_reason_and_a_use(self):
        """An allow-list entry nothing uses any more is a hole left open for the next edit."""
        blocks, _ = shell_blocks()
        every_field = set()
        for _, text, _ in blocks:
            every_field.update(f for _, f, _, _ in string.Formatter().parse(text) if f)
        for name, reason in ALLOWED_UNQUOTED.items():
            self.assertTrue(reason.strip(), name)
            self.assertIn(name, every_field, "%s is allowed but no shell string uses it" % name)


# ------------------------------------------------------------------ commands sent over ssh
# Paths that break a command written with a hand-made quote, in two ways: one leaves the quotes
# unbalanced (the command is a syntax error), the other balances them and runs `touch`.
HOSTILE_PATHS = (
    "/srv/it's here/$(touch PWN)/a;touch PWN2 `touch PWN3` \"q\"",
    "/srv/x'; touch PWN; echo '/y",
)

# What a remote shell finds: every program is a stub that notes its arguments, and anything else
# (touch, here) is noted as UNEXPECTED and not run, because PATH is empty.
STUBS = r"""
rec() { { printf '%s\0' "$@"; printf '\036'; } >> "$RECORD"; }
mkdir() { rec mkdir "$@"; }
ls() { rec ls "$@"; }
find() { rec find "$@"; }
awk() { rec awk "$@"; }
command_not_found_handle() { rec UNEXPECTED "$@"; return 127; }
"""

HARNESS = r"""
ssh() {
    local cmd="${!#}"
    PATH=/nonexistent "$BASH" --noprofile --norc -c "source \"$STUBS\"; $cmd" || true
}
"""


def run_remote(script_body, hostile, **env):
    """Run bash code that calls ssh. Returns the programs the 'remote' shell ran, each as the list
    [program, arguments...], and the finished process."""
    with tempfile.TemporaryDirectory() as tmp:
        stubs = os.path.join(tmp, "stubs.sh")
        record = os.path.join(tmp, "record")
        with open(stubs, "w") as handle:
            handle.write(STUBS)
        open(record, "w").close()
        env = dict(os.environ, STUBS=stubs, RECORD=record, HOSTILE=hostile, **env)
        result = subprocess.run(["bash", "-c", HARNESS + script_body], env=env, cwd=tmp,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(record, "rb") as handle:
            raw = handle.read()
    calls = [part.split(b"\0")[:-1] for part in raw.split(b"\x1e") if part]
    return [[word.decode("utf-8") for word in call] for call in calls], result


def matching_lines(path, pattern):
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if re.search(pattern, line)]


class TestRemoteStrings(unittest.TestCase):
    """A path inside a command string that a remote shell reads has to be quoted for that shell."""

    def message(self, result):
        return result.stderr.decode("utf-8", "replace")

    def test_final_script_quotes_the_remote_path(self):
        script = os.path.join(REPO, "backend", "workflow", "scripts", "final.script.sh")
        lines = matching_lines(script, r'^\s*ssh "\$\{SSH_ARGS\[@\]\}" "\$\{REMOTE_HOST\}" "(mkdir|ls) ')
        self.assertEqual(len(lines), 2, "expected the mkdir and ls commands sent to the results server")
        body = ('SSH_ARGS=(); REMOTE_HOST=host\nREMOTE_PATH="$HOSTILE/"\n'
                'REMOTE_FILE_PATH="${REMOTE_PATH%/}/archive.zip"\n' + "\n".join(lines) + "\n")
        for hostile in HOSTILE_PATHS:
            with self.subTest(path=hostile):
                calls, result = run_remote(body, hostile)
                self.assertEqual(calls, [["mkdir", "-p", hostile], ["ls", "-lh", hostile + "/archive.zip"]],
                                 self.message(result))

    def test_ensure_remote_directory_quotes_the_remote_path(self):
        common = os.path.join(REPO, "frontend", "scripts", "common.sh")
        body = 'source "$COMMON"\nset +eu\nensure_remote_directory "user@host:$HOSTILE/"\n'
        for hostile in HOSTILE_PATHS:
            with self.subTest(path=hostile):
                calls, result = run_remote(body, hostile, COMMON=common)
                self.assertEqual(calls, [["mkdir", "-p", hostile]], self.message(result))

    def test_goserver_quotes_the_path_in_the_size_probe(self):
        common = os.path.join(REPO, "frontend", "scripts", "common.sh")
        script = os.path.join(REPO, "frontend", "scripts", "goServer.bash")
        lines = matching_lines(script, r'^INPUT_BYTES="\$\(ssh "\$\{DTN_HOST\}" "find ')
        self.assertEqual(len(lines), 1, "expected the one size probe that runs find on the DTN")
        body = ('source "$COMMON"\nset +eu\n'
                'DTN_HOST=dtn; PIPELINE_ROOT="$HOSTILE"; USER_ID=u; RUN_NAME=r\n' + lines[0] + "\n")
        for hostile in HOSTILE_PATHS:
            with self.subTest(path=hostile):
                calls, result = run_remote(body, hostile, COMMON=common)
                self.assertEqual(calls, [["find", hostile + "/data/runs/u/r", "-maxdepth", "2", "-type",
                                          "f", "!", "-path", "*/@eaDir/*", "-printf", "%f\\t%s\\n"],
                                         ["awk", "-F", "\\t",
                                          "$1 ~ /\\.fastq(\\.gz)?$/ {sum += $2} END {print sum + 0}"]],
                                 self.message(result))


class TestSbatchAccount(unittest.TestCase):
    """The SLURM account comes from config/cluster.json (__default__.account), not from the script.

    snakemake.sh passes `--account <value>` to sbatch only when the value is not empty: the lab that
    wrote the workflow charges one account, and a cluster whose partitions need none leaves it empty.
    Snakemake runs the --cluster string through a shell, so the value is checked and quoted on the way.
    These tests run the real snakemake.sh with a stand-in `snakemake` that records its arguments."""

    # the --cluster string as it was when the account was written into the script
    CLUSTER = ("sbatch {account}--partition {{cluster.partition}} --job-name {{cluster.name}} "
               "-o {{cluster.output}} -e {{cluster.error}} --ntasks-per-node 1 --time {{cluster.time}} "
               "--mem {{cluster.mem}} --cpus-per-task {{cluster.cpus_per_task}}")

    def run_script(self, cluster_json=None, executor=None):
        """(process, the arguments snakemake was given or None, a canary file existed) for one run
        of snakemake.sh in a scratch backend whose config/cluster.json is `cluster_json` (a dict, a
        string written as it is, or None for no file)."""
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(os.path.join(REPO, "backend", "snakemake.sh"), os.path.join(tmp, "snakemake.sh"))
            os.makedirs(os.path.join(tmp, "config"))
            if cluster_json is not None:
                with open(os.path.join(tmp, "config", "cluster.json"), "w", encoding="utf-8") as handle:
                    handle.write(cluster_json if isinstance(cluster_json, str) else json.dumps(cluster_json))
            bin_dir = os.path.join(tmp, "bin")
            os.makedirs(bin_dir)
            recorded = os.path.join(tmp, "snakemake.args")
            fake = os.path.join(bin_dir, "snakemake")
            with open(fake, "w") as handle:
                handle.write('#!/bin/sh\nprintf \'%s\\0\' "$@" > "$RECORDED_ARGS"\n')
            os.chmod(fake, 0o755)
            env = {k: v for k, v in os.environ.items() if not k.startswith("TBRC_")}
            env.update(PATH=bin_dir + os.pathsep + env.get("PATH", ""), RECORDED_ARGS=recorded)
            if executor:
                env["TBRC_EXECUTOR"] = executor
            result = subprocess.run(["bash", os.path.join(tmp, "snakemake.sh")], env=env, cwd=tmp,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            args = None
            if os.path.exists(recorded):
                with open(recorded, "rb") as handle:
                    args = [a.decode("utf-8") for a in handle.read().split(b"\0")[:-1]]
            return result, args, os.path.exists(os.path.join(tmp, "canary"))

    def cluster_string(self, args):
        self.assertIsNotNone(args, "snakemake was not run")
        self.assertEqual(args.count("--cluster"), 1, args)
        return args[args.index("--cluster") + 1]

    def test_a_named_account_is_charged_on_the_line_as_it_always_was(self):
        result, args, _ = self.run_script({"__default__": {"account": "lab_account_1", "partition": "p"}})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.cluster_string(args), self.CLUSTER.format(account="--account lab_account_1 "))

    def test_an_empty_or_missing_account_adds_no_option(self):
        for config in ({"__default__": {"account": "", "partition": "p"}}, {"__default__": {"partition": "p"}},
                       {"__default__": {"account": None, "partition": "p"}}):
            with self.subTest(config=config):
                result, args, _ = self.run_script(config)
                self.assertEqual(result.returncode, 0, result.stderr)
                line = self.cluster_string(args)
                self.assertEqual(line, self.CLUSTER.format(account=""))
                self.assertNotIn("--account", line)

    def test_a_value_that_is_not_a_plain_name_is_refused_before_snakemake_starts(self):
        hostile = ["x; touch canary", "$(touch canary)", "`touch canary`", "a b", "-x", "a\nb", "x'y", 'x"y',
                   "a*", "a|b", "a&&b", "x\\y", "caf\u00e9", 7, ["a"]]
        for value in hostile:
            with self.subTest(account=value):
                result, args, canary = self.run_script({"__default__": {"account": value, "partition": "p"}})
                self.assertNotEqual(result.returncode, 0)
                self.assertIsNone(args, "snakemake must not start with this account")
                self.assertFalse(canary, "the value reached a shell")
                message = result.stderr.decode("utf-8", "replace")
                self.assertIn("account", message)
                self.assertNotIn("canary", message, "the message must not repeat the value")

    def test_an_unreadable_cluster_json_stops_the_run_rather_than_dropping_the_account(self):
        for config in (None, "{not json", "[]", {"partition": "p"}):
            with self.subTest(config=config):
                result, args, _ = self.run_script(config)
                self.assertNotEqual(result.returncode, 0)
                self.assertIsNone(args)
                self.assertIn("cluster.json", result.stderr.decode("utf-8", "replace"))

    def test_a_run_on_this_host_does_not_read_the_account(self):
        # TBRC_EXECUTOR=local has no sbatch, so a broken or missing cluster.json must not matter
        for config in (None, {"__default__": {"account": "x; touch canary"}}):
            with self.subTest(config=config):
                result, args, canary = self.run_script(config, executor="local")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("--cluster", args)
                self.assertFalse(canary)

    def test_the_shipped_cluster_json_gives_the_line_for_its_own_account(self):
        path = os.path.join(REPO, "backend", "config", "cluster.json")
        with open(path, encoding="utf-8") as handle:
            shipped = json.load(handle)
        account = shipped["__default__"].get("account")
        self.assertIsInstance(account, str, "cluster.json __default__ needs an account key (empty is allowed)")
        result, args, _ = self.run_script(shipped)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.cluster_string(args),
                         self.CLUSTER.format(account="--account %s " % account if account else ""))


class TestPandaseqIsRepeatable(unittest.TestCase):
    """classic_luka's output has to be the same every run for the same reads. pandaseq's default
    threading writes the same reads in a different order each run, and fastx_collapser breaks ties
    by input order, so `head` kept different equally common reads. One thread makes the order fixed."""

    def pandaseq_command(self):
        blocks, _ = shell_blocks()
        found = [text for label, text, _ in blocks if label == "pandaseq.smk rule pandaseq"]
        self.assertEqual(len(found), 1, "expected one shell string in rule pandaseq")
        return found[0]

    def test_pandaseq_runs_on_one_thread(self):
        words = shlex.split(self.pandaseq_command())
        self.assertEqual(words[0], "pandaseq")
        self.assertIn("-T", words)
        self.assertEqual(words[words.index("-T") + 1], "1", "-T must be followed by 1: %r" % (words,))

    def test_the_rule_says_why(self):
        with open(os.path.join(RULES, "pandaseq.smk"), encoding="utf-8") as handle:
            comments = [line for line in handle if line.lstrip().startswith("#")]
        self.assertTrue(any("fastx_collapser" in line and "-T 1" in line for line in comments),
                        "a comment in pandaseq.smk should name -T 1 and fastx_collapser: %r" % (comments,))


# ------------------------------------------------------------------ the comparison tools
def dryrun(rule, command, header="input: a, b"):
    return ("Building DAG of jobs...\n\n[Mon Nov  9 04:53:39 2026]\nrule %s:\n    %s\n    output: c\n"
            "    jobid: 1\n    reason: Missing output files: c\n    resources: tmpdir=/tmp\n\n%s\n\n"
            "Job stats:\njob count\n" % (rule, header, command))


class TestCompareDryrun(unittest.TestCase):
    def compare(self, before, after):
        old = compare_dryrun.parse_dryrun(dryrun("r", before), "before")
        new = compare_dryrun.parse_dryrun(dryrun("r", after), "after")
        return compare_dryrun.compare_pipeline("p", old, new)

    def test_quoting_style_alone_is_not_a_difference(self):
        self.assertEqual(self.compare("run -s 'a' -d \"b\"", "run -s a -d b"), [])
        self.assertEqual(self.compare('cp "$W/m_"igh"_D.nhr" x', 'cp "$W/m_igh_D.nhr" x'), [])

    def test_an_argument_lost_is_a_difference(self):
        problems = self.compare("run -t '' -y x", "run -t -y x")
        self.assertTrue(problems and "command differs" in problems[0], problems)
        self.assertTrue(any("''" in line for line in problems[1:]), problems)

    def test_a_different_value_is_a_difference(self):
        self.assertTrue(self.compare("run -s a", "run -s b"))

    def test_comments_with_an_apostrophe_are_skipped(self):
        self.assertEqual(self.compare("# the 5' end\nrun a", "run a"), [])

    def test_temp_file_notes_are_not_the_command_and_their_order_does_not_matter(self):
        def run(notes):
            return compare_dryrun.parse_dryrun(dryrun("r", "run a\n\n" + notes), "x")
        first = run("Would remove temporary output f1\nWould remove temporary output f2")
        again = run("Would remove temporary output f2\nWould remove temporary output f1")
        other = run("Would remove temporary output f1\nWould remove temporary output f3")
        self.assertEqual(first[0][2], ["run", "a"])
        self.assertEqual(compare_dryrun.compare_pipeline("p", first, again), [])
        self.assertTrue(any("temporary differs" in line
                            for line in compare_dryrun.compare_pipeline("p", first, other)))

    def test_a_changed_input_list_is_a_difference(self):
        old = compare_dryrun.parse_dryrun(dryrun("r", "run"), "b")
        new = compare_dryrun.parse_dryrun(dryrun("r", "run", header="input: a, c"), "a")
        self.assertTrue(any("input differs" in line for line in compare_dryrun.compare_pipeline("p", old, new)))

    def test_folders_exit_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            for side, command in (("before", "run -t ''"), ("after", "run -t")):
                os.makedirs(os.path.join(tmp, side))
                with open(os.path.join(tmp, side, "p.txt"), "w") as handle:
                    handle.write(dryrun("r", command))
            same, differing = compare_dryrun.compare_folders(os.path.join(tmp, "before"),
                                                             os.path.join(tmp, "before"))
            self.assertEqual(differing, 0, same)
            report, differing = compare_dryrun.compare_folders(os.path.join(tmp, "before"),
                                                               os.path.join(tmp, "after"))
            self.assertEqual(differing, 1, report)
            with contextlib.redirect_stdout(io.StringIO()):
                status = compare_dryrun.main(["x", os.path.join(tmp, "before"), os.path.join(tmp, "after")])
            self.assertEqual(status, 1)


class TestCompareArchives(unittest.TestCase):
    def test_lists_missing_and_changed_files_in_zip_and_tar(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = os.path.join(tmp, "a.zip")
            second = os.path.join(tmp, "b.tar.gz")
            with zipfile.ZipFile(first, "w") as archive:
                archive.writestr("run/same.txt", "same")
                archive.writestr("run/changed.txt", "one")
                archive.writestr("run/only_a.txt", "a")
            os.makedirs(os.path.join(tmp, "src", "run"))
            for name, text in (("same.txt", "same"), ("changed.txt", "two"), ("only_b.txt", "b")):
                with open(os.path.join(tmp, "src", "run", name), "w") as handle:
                    handle.write(text)
            with tarfile.open(second, "w:gz") as archive:
                archive.add(os.path.join(tmp, "src", "run"), arcname="./run")
            only_a, only_b, different = compare_archives.compare(compare_archives.digests(first),
                                                                 compare_archives.digests(second))
            self.assertEqual((only_a, only_b, different),
                             (["run/only_a.txt"], ["run/only_b.txt"], ["run/changed.txt"]))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(compare_archives.main(["x", first, second]), 1)
                self.assertEqual(compare_archives.main(["x", first, first]), 0)


if __name__ == "__main__":
    unittest.main()
