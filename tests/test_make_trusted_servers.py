"""Tests for scripts/make_trusted_servers.py, which turns frontend/server_config.json into the
trusted list the cluster checks a sample.json against.

Run from the repository root:  python3 -m unittest tests/test_make_trusted_servers.py -v
"""
import copy
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPT = os.path.join(REPO, "scripts", "make_trusted_servers.py")
sys.path.insert(0, os.path.join(REPO, "scripts"))

import make_trusted_servers  # noqa: E402

REAL_CONFIG = os.path.join(REPO, "frontend", "server_config.json")
REAL_LIST = os.path.join(REPO, "backend", "config", "trusted_servers.json")
# The public export has no server_config.json and no list (the list is optional); it ships this example.
EXAMPLE_CONFIG = os.path.join(REPO, "frontend", "server_config.example.json")

# a config of the same shape as the real one, with values no real host uses
FAKE_CONFIG = {
    "storage_server": {
        "boxA": {
            "rawdata_path": "me@boxA:/vol/raw/",
            "results_path": "me@boxA:/vol/res/",
            "ssh_port": 2222,
            "web_port": 4000,
            "transfer_mode": "push",
        },
        "boxB": {
            "rawdata_path": "you@boxB:/data/raw/",
            "results_path": "you@boxB:/data/res/",
            "ssh_port": 22,
            "web_port": 5000,
            "transfer_mode": "pull",
        },
    },
    "cluster": {
        "dtn_host": "u@dtn.example",
        "path_prefix": "/opt/conda/condabin:/usr/bin",
        "snakemake_env": "smk",
        "igblast_bin": "igblastn",
        "igblast_data": "/pipelines/igblast/internal_data",
        "igblast_organism": "human",
        "igblast_db_v": "",
        "igblast_db_d": "",
        "igblast_db_j": "",
        "igblast_auxiliary_data": "/pipelines/igblast/optional_file",
    },
    "metrics": {"manual_minutes_per_100mb": 29},
}

LEGACY_IGBLAST_DATA = [
    "/path/to/home/Downloads/ncbi-igblast-1.22.0/internal_data",
    "/path/to/scratch/pipelines/TBRC/igblast/internal_data",
]
LEGACY_IGBLAST_AUX = [
    "/path/to/home/Downloads/ncbi-igblast-1.22.0/optional_file",
    "/path/to/scratch/pipelines/TBRC/igblast/optional_file",
]
# A database prefix two finished runs (a mouse IgH D set) still carry; the current config has none.
LEGACY_IGBLAST_DB_D = ["/path/to/scratch/pipelines/TBRC/igblast/work/mouse/igh/mouse_igh_D"]


def run_script(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


class TrustedServersTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = self._tmp.name

    def path(self, name):
        return os.path.join(self.tmp, name)

    def write_config(self, config, name="server_config.json"):
        p = self.path(name)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(config, f)
        return p

    def generate(self, config_path, out_path):
        result = run_script("--config", config_path, "--out", out_path)
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(out_path, encoding="utf-8") as f:
            return json.load(f)

    def test_generates_exact_fields_for_each_server(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        trusted = self.generate(cfg_path, self.path("trusted.json"))
        self.assertEqual(sorted(trusted["servers"]), ["boxA", "boxB"])
        for name, entry in FAKE_CONFIG["storage_server"].items():
            got = trusted["servers"][name]
            self.assertEqual(
                sorted(got), ["rawdata_path", "results_path", "ssh_port", "transfer_mode"],
                "web_port and anything else the cluster never sees must not be listed")
            for field in ("rawdata_path", "results_path", "ssh_port", "transfer_mode"):
                self.assertEqual(got[field], entry[field], "%s.%s" % (name, field))
        self.assertIsInstance(trusted["servers"]["boxA"]["ssh_port"], int)
        self.assertEqual(trusted["cluster"], {
            "path_prefix": FAKE_CONFIG["cluster"]["path_prefix"],
            "snakemake_env": "smk",
        })

    def test_build_in_process_matches_the_command_line(self):
        trusted = self.generate(self.write_config(FAKE_CONFIG), self.path("trusted.json"))
        self.assertEqual(make_trusted_servers.build(FAKE_CONFIG), trusted)

    def real_pair(self):
        """(config, list) as parsed JSON: this tree's own server_config.json and trusted_servers.json, or,
        in a tree that has no list (the public export), the example config and the list a new installation
        generates from it with this script. Skips when the tree has neither a config nor the example."""
        config_path = REAL_CONFIG if os.path.exists(REAL_CONFIG) else EXAMPLE_CONFIG
        if not os.path.exists(config_path):
            self.skipTest("this tree has neither frontend/server_config.json nor server_config.example.json")
        list_path = REAL_LIST
        if not os.path.exists(list_path):
            list_path = self.path("generated_trusted_servers.json")
            result = run_script("--config", config_path, "--out", list_path)
            self.assertEqual(result.returncode, 0, result.stderr)
        with open(config_path, encoding="utf-8") as f:
            config = json.load(f)
        with open(list_path, encoding="utf-8") as f:
            return config, json.load(f)

    def test_real_config_matches_the_committed_file(self):
        """The committed list is what the real config generates, field for field."""
        real, committed = self.real_pair()
        self.assertEqual(sorted(committed["servers"]), sorted(real["storage_server"]))
        for name, entry in real["storage_server"].items():
            for field in ("rawdata_path", "results_path", "ssh_port", "transfer_mode"):
                self.assertEqual(committed["servers"][name][field], entry[field], "%s.%s" % (name, field))
        self.assertEqual(committed["cluster"]["path_prefix"], real["cluster"]["path_prefix"])
        self.assertEqual(committed["cluster"]["snakemake_env"], real["cluster"]["snakemake_env"])

    def test_igblast_lists_hold_current_and_legacy_paths(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        trusted = self.generate(cfg_path, self.path("trusted.json"))
        ig = trusted["igblast"]
        self.assertEqual(sorted(ig), ["igblast_auxiliary_data", "igblast_bin", "igblast_data",
                                      "igblast_db_d", "igblast_db_j", "igblast_db_v"])
        for key, values in ig.items():
            self.assertIsInstance(values, list, key)
        # today's value comes first, then the two older install locations, then the empty value that
        # means "the copy bundled with the backend" (CLI --local and benchmark runs write null there)
        self.assertEqual(ig["igblast_data"], ["/pipelines/igblast/internal_data"] + LEGACY_IGBLAST_DATA + [""])
        self.assertEqual(ig["igblast_auxiliary_data"],
                         ["/pipelines/igblast/optional_file"] + LEGACY_IGBLAST_AUX + [""])
        self.assertEqual(ig["igblast_bin"], ["igblastn"])
        for key in ("igblast_db_v", "igblast_db_j"):
            self.assertEqual(ig[key], [""], key)
        self.assertEqual(ig["igblast_db_d"], [""] + LEGACY_IGBLAST_DB_D)

    def test_igblast_lists_for_the_real_config(self):
        """The real file lists today's paths, the legacy locations and the empty value, each exactly once."""
        config, committed = self.real_pair()
        real, ig = config["cluster"], committed["igblast"]
        self.assertEqual(ig["igblast_data"], [real["igblast_data"]] + LEGACY_IGBLAST_DATA + [""])
        self.assertEqual(ig["igblast_auxiliary_data"],
                         [real["igblast_auxiliary_data"]] + LEGACY_IGBLAST_AUX + [""])
        self.assertEqual(ig["igblast_bin"], [real["igblast_bin"]])
        self.assertEqual(ig["igblast_db_d"], [real["igblast_db_d"]] + LEGACY_IGBLAST_DB_D)
        self.assertEqual(ig["igblast_db_v"], [real["igblast_db_v"]])
        self.assertEqual(ig["igblast_db_j"], [real["igblast_db_j"]])
        for values in ig.values():
            self.assertEqual(len(values), len(set(values)))

    def test_a_current_path_equal_to_a_legacy_one_is_listed_once(self):
        cfg = copy.deepcopy(FAKE_CONFIG)
        cfg["cluster"]["igblast_data"] = LEGACY_IGBLAST_DATA[1]
        trusted = self.generate(self.write_config(cfg), self.path("trusted.json"))
        self.assertEqual(trusted["igblast"]["igblast_data"],
                         [LEGACY_IGBLAST_DATA[1], LEGACY_IGBLAST_DATA[0], ""])

    def test_output_is_sorted_and_stable(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        out = self.path("trusted.json")
        self.generate(cfg_path, out)
        with open(out, encoding="utf-8") as f:
            first = f.read()
        self.generate(cfg_path, out)
        with open(out, encoding="utf-8") as f:
            self.assertEqual(f.read(), first)
        self.assertTrue(first.endswith("\n"))
        self.assertEqual(first, json.dumps(json.loads(first), indent=2, sort_keys=True) + "\n")

    def test_check_mode_passes_when_in_step(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        out = self.path("trusted.json")
        self.generate(cfg_path, out)
        result = run_script("--check", "--config", cfg_path, "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_check_mode_fails_when_file_is_stale(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        out = self.path("trusted.json")
        self.generate(cfg_path, out)
        stale = copy.deepcopy(FAKE_CONFIG)
        stale["storage_server"]["boxA"]["results_path"] = "me@boxA:/vol/elsewhere/"
        result = run_script("--check", "--config", self.write_config(stale), "--out", out)
        self.assertEqual(result.returncode, 1)
        self.assertIn("make_trusted_servers.py", result.stderr)
        # check mode never rewrites the file
        with open(out, encoding="utf-8") as f:
            self.assertIn("/vol/res/", f.read())

    def test_check_mode_fails_when_file_is_missing(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        result = run_script("--check", "--config", cfg_path, "--out", self.path("absent.json"))
        self.assertEqual(result.returncode, 1)
        self.assertFalse(os.path.exists(self.path("absent.json")))

    def test_check_mode_fails_when_the_file_was_edited_by_hand(self):
        cfg_path = self.write_config(FAKE_CONFIG)
        out = self.path("trusted.json")
        self.generate(cfg_path, out)
        with open(out, encoding="utf-8") as f:
            edited = f.read().replace('"boxB"', '"boxC"')
        with open(out, "w", encoding="utf-8") as f:
            f.write(edited)
        result = run_script("--check", "--config", cfg_path, "--out", out)
        self.assertEqual(result.returncode, 1)

    def test_old_string_form_fails_clearly(self):
        cfg = copy.deepcopy(FAKE_CONFIG)
        cfg["storage_server"]["boxA"] = "me@boxA:/vol/raw/"
        out = self.path("trusted.json")
        result = run_script("--config", self.write_config(cfg), "--out", out)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("boxA", result.stderr)
        self.assertFalse(os.path.exists(out), "nothing is written from a config it cannot read")

    def test_a_missing_field_fails_clearly(self):
        for field in ("rawdata_path", "results_path", "ssh_port", "transfer_mode"):
            cfg = copy.deepcopy(FAKE_CONFIG)
            del cfg["storage_server"]["boxB"][field]
            result = run_script("--config", self.write_config(cfg), "--out", self.path("t.json"))
            self.assertNotEqual(result.returncode, 0, field)
            self.assertIn("boxB", result.stderr, field)
            self.assertIn(field, result.stderr, field)

    def test_a_wrongly_typed_field_fails_clearly(self):
        cfg = copy.deepcopy(FAKE_CONFIG)
        cfg["storage_server"]["boxA"]["ssh_port"] = "2222"
        result = run_script("--config", self.write_config(cfg), "--out", self.path("t.json"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ssh_port", result.stderr)

    def test_a_missing_cluster_field_fails_clearly(self):
        for field in ("path_prefix", "snakemake_env", "igblast_bin", "igblast_data", "igblast_auxiliary_data",
                      "igblast_db_v", "igblast_db_d", "igblast_db_j"):
            cfg = copy.deepcopy(FAKE_CONFIG)
            del cfg["cluster"][field]
            result = run_script("--config", self.write_config(cfg), "--out", self.path("t.json"))
            self.assertNotEqual(result.returncode, 0, field)
            self.assertIn(field, result.stderr, field)

    def test_no_servers_fails_clearly(self):
        cfg = copy.deepcopy(FAKE_CONFIG)
        cfg["storage_server"] = {}
        result = run_script("--config", self.write_config(cfg), "--out", self.path("t.json"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("storage_server", result.stderr)

    def test_unreadable_config_fails_clearly(self):
        result = run_script("--config", self.path("nope.json"), "--out", self.path("t.json"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nope.json", result.stderr)

    def test_defaults_work_from_any_directory(self):
        """With no options --check reads the repository's own config and list, wherever it is run from.

        A tree without a list (the public export) is given one: a scratch repository holding this script,
        the example config as its server_config.json and the list generated from it."""
        script = SCRIPT
        if not (os.path.exists(REAL_CONFIG) and os.path.exists(REAL_LIST)):
            if not os.path.exists(EXAMPLE_CONFIG):
                self.skipTest("this tree has no trusted list and no server_config.example.json to build one from")
            scratch = self.path("scratch_repo")
            for folder in ("scripts", "frontend", os.path.join("backend", "config")):
                os.makedirs(os.path.join(scratch, folder))
            script = os.path.join(scratch, "scripts", "make_trusted_servers.py")
            shutil.copy(SCRIPT, script)
            shutil.copy(EXAMPLE_CONFIG, os.path.join(scratch, "frontend", "server_config.json"))
            made = subprocess.run([sys.executable, script], capture_output=True, text=True, cwd=self.tmp)
            self.assertEqual(made.returncode, 0, made.stderr)
        result = subprocess.run([sys.executable, script, "--check"], capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)


class ValidationCasesTest(unittest.TestCase):
    """The shared case table has the shape Tasks 3, 4 and 6 read."""

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(REPO, "tests", "validation_cases.json"), encoding="utf-8") as f:
            cls.cases = json.load(f)

    def test_top_level_keys(self):
        self.assertEqual(sorted(self.cases), ["email", "hostile", "method", "path_segment", "trim_sequence"])

    def test_path_segment_groups(self):
        seg = self.cases["path_segment"]
        self.assertEqual(sorted(seg), ["bad", "bad_untrimmed", "good"])
        self.assertIn("231003_M03851_0285_000000000-DLPWD", seg["good"])
        self.assertEqual(seg["bad_untrimmed"], ["a\n", " a"])
        for needed in ("", ".", "..", "-x", "a\nb", "café", "１２"):
            self.assertIn(needed, seg["bad"])

    def test_trim_sequence_groups(self):
        trim = self.cases["trim_sequence"]
        self.assertEqual(sorted(trim), ["bad", "good"])
        lengths = sorted(len(s) for s in trim["bad"])
        self.assertIn(9, lengths)
        self.assertIn(65, lengths)

    def test_method_groups(self):
        self.assertEqual(self.cases["method"]["good"],
                         ["classic_luka", "tiago_vsearch", "tiago_nanopore", "tiago_sulkseq"])
        self.assertTrue(self.cases["method"]["bad"])

    def test_email_groups(self):
        email = self.cases["email"]
        self.assertEqual(sorted(email), ["bad", "good"])
        self.assertTrue(email["good"])
        for needed in ("-a@b.co", "--delete@b.co"):     # a local part may not start with "-"
            self.assertIn(needed, email["bad"])
        self.assertIn("a-b@c.co", email["good"])        # a "-" after the first character is fine

    def test_hostile_values(self):
        self.assertEqual(len(self.cases["hostile"]), 8)
        self.assertIn("\n touch canary", self.cases["hostile"])
        self.assertIn("--delete", self.cases["hostile"])

    def test_the_documented_rules_agree_with_the_table(self):
        """Path segment `^[A-Za-z0-9._-]+$` (not . or .., no leading -); trim `^[ACGTRYSWKMBDHVN]{10,64}$`."""
        seg_ok = lambda s: (re.fullmatch(r"[A-Za-z0-9._-]+", s, re.ASCII) is not None
                            and s not in (".", "..") and not s.startswith("-"))
        trim_ok = lambda s: re.fullmatch(r"[ACGTRYSWKMBDHVN]{10,64}", s, re.ASCII) is not None
        seg = self.cases["path_segment"]
        self.assertTrue(all(seg_ok(s) for s in seg["good"]))
        self.assertFalse(any(seg_ok(s) for s in seg["bad"] + seg["bad_untrimmed"]))
        trim = self.cases["trim_sequence"]
        self.assertTrue(all(trim_ok(s) for s in trim["good"]))
        self.assertFalse(any(trim_ok(s) for s in trim["bad"]))
        self.assertFalse(any(seg_ok(s) for s in self.cases["hostile"]))
        self.assertFalse(any(trim_ok(s) for s in self.cases["hostile"]))
        self.assertFalse(any(s in self.cases["method"]["good"] for s in self.cases["hostile"]))


if __name__ == "__main__":
    unittest.main()
