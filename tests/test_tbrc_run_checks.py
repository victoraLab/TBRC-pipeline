"""`./tbrc run` applies the workflow's own sample.json checks (backend/workflow/sample_rules.py) to the
row it has built, before it stages a FASTQ or reaches a remote machine, and stops with the message that
names the field. Until then the CLI accepted a lower-case trim sequence, 600 contigs or an address
starting with "-" and the workflow refused each one only after the files were copied.

The tests run `tbrc` in-process against a scratch tree (config, run folders), with the commands that
would touch a machine replaced by ones that fail the test if they are reached.

Run from the repository root:  python3 -m unittest tests/test_tbrc_run_checks.py -v
"""
import contextlib
import importlib.machinery
import importlib.util
import io
import json
import os
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TBRC = os.path.join(REPO, "tbrc")
CONFIG_FILE = os.path.join(REPO, "frontend", "server_config.json")
EXAMPLE_CONFIG = os.path.join(REPO, "frontend", "server_config.example.json")


def load_tbrc():
    loader = importlib.machinery.SourceFileLoader("tbrc_cli_for_run_checks", TBRC)
    spec = importlib.util.spec_from_loader("tbrc_cli_for_run_checks", loader)
    module = importlib.util.module_from_spec(spec)
    keep = sys.dont_write_bytecode
    sys.dont_write_bytecode = True  # no tbrc .pyc dropped into the repository root
    try:
        loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = keep
    return module


class RunChecksTest(unittest.TestCase):
    def setUp(self):
        config_path = CONFIG_FILE if os.path.exists(CONFIG_FILE) else EXAMPLE_CONFIG
        if not os.path.exists(config_path):
            self.skipTest("this tree has neither frontend/server_config.json nor server_config.example.json")
        with open(config_path, encoding="utf-8") as handle:
            self.config = json.load(handle)
        self.server = sorted(self.config["storage_server"])[0]
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = self._tmp.name
        self.cli = load_tbrc()
        # a scratch backend that shares the real workflow code, with no trusted list and no run folders
        self.backend = os.path.join(self.tmp, "backend")
        os.makedirs(os.path.join(self.backend, "config"))
        os.symlink(os.path.join(REPO, "backend", "workflow"), os.path.join(self.backend, "workflow"))
        self.runs = os.path.join(self.tmp, "runs")
        self.cli.BACKEND = self.backend
        self.cli.RUNS_DIR = self.runs
        self.write_config(self.config)
        self.cli.local_runner = lambda a: ("native", None)
        for name in ("run", "ssh", "stage_fastqs"):
            setattr(self.cli, name, self.unreachable(name))

    def unreachable(self, name):
        def call(*args, **kwargs):
            raise AssertionError("tbrc reached %s() after a sample.json the workflow would refuse" % name)
        return call

    def write_config(self, config):
        self.cli.CONFIG_FILE = os.path.join(self.tmp, "server_config.json")
        with open(self.cli.CONFIG_FILE, "w", encoding="utf-8") as handle:
            json.dump(config, handle)

    def run_cli(self, *args, **kwargs):
        """(exit status, printed text). The input is a local folder unless `inputs` gives another."""
        inputs = kwargs.pop("inputs", ["--fastq-dir", self.tmp])
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = self.cli.main(["run", "--user", "alice", "--run-name", "plate7", "--primer-set", "B2"]
                                 + inputs + list(args))
        return code, out.getvalue()

    def assertRefused(self, args, field, option=None, secret=None, **kwargs):
        code, text = self.run_cli(*args, **kwargs)
        self.assertEqual(code, 2, text)
        self.assertIn("error: sample.json refused", text)
        self.assertIn("\n  %s: " % field, text)
        self.assertIn("hint:", text)
        if option:
            self.assertIn(option, text.split("hint:", 1)[1])
        if secret:
            self.assertNotIn(secret, text, "the message must not repeat the value")
        self.assertFalse(os.path.exists(self.runs), "a run folder was made")
        self.assertFalse(os.path.exists(os.path.join(self.backend, "data")), "a run folder was made")
        return text

    # ---------------------------------------------------------- what the workflow refuses, refused here
    def test_a_lower_case_or_short_or_long_trim_sequence_is_refused(self):
        for value in ("acgtacgtacgtacgt", "GGGAATTCG", "A" * 65, "GGGAATTCGAGG TGC"):
            with self.subTest(trim_seq=value):
                text = self.assertRefused(["--local", "--trim-seq", value], "trim_seq", "--trim-seq", value)
                self.assertIn("10 to 64 letters", text)

    def test_contigs_outside_1_to_500_are_refused(self):
        for value in ("0", "501", "600", "-3"):
            with self.subTest(contigs=value):
                text = self.assertRefused(["--local", "--contigs=" + value], "sequence_number", "--contigs")
                self.assertIn("1 to 500", text.split("hint:", 1)[1])

    def test_an_email_the_workflow_would_refuse_is_refused(self):
        for value in ("-x@lab.org", "a b@lab.org", "no-at-sign", "x@lab", "café@lab.org"):
            with self.subTest(email=value):
                self.assertRefused(["--local", "--email=" + value], "user_email", "--email", value)

    def test_names_are_one_path_segment(self):
        # a set with no entry in barcode_defaults needs --trim-seq to get as far as the workflow's rule
        self.assertRefused(["--local", "--primer-set", "B2;x", "--trim-seq", "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"],
                           "primer_set", "--primer-set", "B2;x")
        self.assertRefused(["--local", "--sequencing-kit", "a b"], "sequencing_kit", "--sequencing-kit")
        self.assertRefused(["--local", "--run-group", "../x"], "run_group", "--run-group")

    def test_a_storage_server_whose_paths_are_not_plain_is_refused(self):
        bad = json.loads(json.dumps(self.config))
        bad["storage_server"][self.server]["results_path"] = "me@host:/my results/"
        self.write_config(bad)
        self.assertRefused([], "server_results_path", "server_config.json",
                           inputs=["--server", self.server, "--folder", "plate7"])

    def test_a_server_that_is_not_in_the_trusted_list_is_refused(self):
        # a list that names a different server: the list is on, so the name must be in it
        trusted = {"servers": {"elsewhere": {"results_path": "me@host:/res/", "ssh_port": 22, "transfer_mode": "push"}},
                   "cluster": {}, "igblast": {k: [""] for k in (
                       "igblast_bin", "igblast_data", "igblast_auxiliary_data", "igblast_db_v", "igblast_db_d",
                       "igblast_db_j")}}
        trusted["igblast"]["igblast_bin"] = [self.config["cluster"].get("igblast_bin", "igblastn")]
        with open(os.path.join(self.backend, "config", "trusted_servers.json"), "w") as handle:
            json.dump(trusted, handle)
        self.assertRefused(["--local"], "server_name", "--server", inputs=["--server", self.server, "--folder", "plate7"])

    def test_several_problems_are_listed_together(self):
        code, text = self.run_cli("--local", "--trim-seq", "acgt", "--contigs=600", "--email=-x@lab.org")
        self.assertEqual(code, 2)
        self.assertIn("3 problems", text)
        for field in ("trim_seq", "sequence_number", "user_email"):
            self.assertIn("\n  %s: " % field, text)

    # ---------------------------------------------------------- and what passes still passes
    def test_ordinary_runs_pass_the_check(self):
        for args in (["--local", "--dry-run"],
                     ["--dry-run"],
                     ["--local", "--dry-run", "--trim-seq", "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"],
                     ["--local", "--dry-run", "--contigs=500", "--email", "a.b+c@lab.example.org"],
                     ["--local", "--dry-run", "--no-trim-primers", "--pipeline", "classic_luka", "--species", "human"]):
            with self.subTest(args=args):
                code, text = self.run_cli(*args)
                self.assertEqual(code, 0, text)

    def test_a_run_that_names_a_storage_server_passes_the_check(self):
        for extra in ([], ["--local"]):
            with self.subTest(extra=extra):
                code, text = self.run_cli("--dry-run", *extra, inputs=["--server", self.server, "--folder", "plate7"])
                self.assertEqual(code, 0, text)

    def test_a_copy_of_tbrc_without_the_workflow_code_still_runs_and_says_it_did_not_check(self):
        # the launcher tests run tbrc from a folder holding only a config and a stub launcher
        os.unlink(os.path.join(self.backend, "workflow"))
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code, text = self.run_cli("--dry-run", "--trim-seq", "acgt")
        self.assertEqual(code, 0, text)
        self.assertIn("was not checked", err.getvalue())

    def test_the_check_runs_before_the_run_folder_is_looked_at(self):
        # an existing run folder is reported after the sample.json check, not instead of it
        os.makedirs(os.path.join(self.backend, "data", "runs", "alice", "plate7"))
        with open(os.path.join(self.backend, "data", "runs", "alice", "plate7", "x"), "w") as handle:
            handle.write("x")
        code, text = self.run_cli("--local", "--trim-seq", "acgt")
        self.assertEqual(code, 2)
        self.assertIn("sample.json refused", text)


if __name__ == "__main__":
    unittest.main()
