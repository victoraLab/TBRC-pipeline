"""Tests for backend/workflow/sample_rules.py, which the Snakefile runs on every row of
sample.json before any job starts.

It assumes the web container may have been hijacked and hands the cluster a doctored sample.json,
so a value that says where data goes (server_results_path, port, mode) or which program runs
(igblast_*) has to equal what the cluster's own trusted list says, not just look well formed.

Run from the repository root:  python3 -m unittest backend/tests/test_sample_rules.py -v
"""
import copy
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WORKFLOW = os.path.join(REPO, "backend", "workflow")
sys.path.insert(0, WORKFLOW)

import sample_rules  # noqa: E402
from sample_rules import (FIELD_RULES, SampleRejected, load_trusted, validate_run_folder,  # noqa: E402
                          validate_sample, validate_sample_name)

with open(os.path.join(REPO, "tests", "validation_cases.json")) as _handle:
    CASES = json.load(_handle)

REAL_TRUSTED_FILE = os.path.join(REPO, "backend", "config", "trusted_servers.json")
# The public export has neither frontend/server_config.json nor the trusted list (the list is optional);
# it ships this example instead, and make_trusted_servers.py turns a config into a list.
EXAMPLE_CONFIG = os.path.join(REPO, "frontend", "server_config.example.json")
sys.path.insert(0, os.path.join(REPO, "scripts"))


def trusted_list_to_test(directory):
    """The path of a trusted list to test against: this tree's own backend/config/trusted_servers.json,
    or, in a tree without one, a list generated from frontend/server_config.example.json into `directory`
    the way a new installation would (python3 scripts/make_trusted_servers.py). None if neither source
    is there."""
    if os.path.exists(REAL_TRUSTED_FILE):
        return REAL_TRUSTED_FILE
    if not os.path.exists(EXAMPLE_CONFIG):
        return None
    import make_trusted_servers
    with open(EXAMPLE_CONFIG, encoding="utf-8") as handle:
        text = make_trusted_servers.render(make_trusted_servers.build(json.load(handle)))
    path = os.path.join(directory, "trusted_servers.json")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path

# Servers and program paths no real host uses, in the shape of backend/config/trusted_servers.json.
TRUSTED = {
    "cluster": {"path_prefix": "/x", "snakemake_env": "snakemake"},
    "igblast": {
        "igblast_auxiliary_data": ["/opt/tbrc/igblast/optional_file", "/old/install/optional_file"],
        "igblast_bin": ["igblastn"],
        "igblast_data": ["/opt/tbrc/igblast/internal_data", "/old/install/internal_data"],
        "igblast_db_d": [""],
        "igblast_db_j": [""],
        "igblast_db_v": [""],
    },
    "servers": {
        "boxA": {"rawdata_path": "me@boxA:/vol/raw/", "results_path": "me@boxA:/vol/res/",
                 "ssh_port": 2222, "transfer_mode": "push"},
        "boxB": {"rawdata_path": "you@boxB:/data/raw/", "results_path": "you@boxB:/data/res/",
                 "ssh_port": 22, "transfer_mode": "push"},
    },
}


def good_row(**changes):
    """A sample.json row of the shape the form writes, for the trusted server boxA."""
    row = {
        "user_id": "u",
        "user_email": "u@example.org",
        "sample_name": "run1",
        "input_method": "Server",
        "server_name": "boxA",
        "server_path": "me@boxA:/vol/raw/",
        "server_results_path": "me@boxA:/vol/res/",
        "server_ssh_port": 2222,
        "server_transfer_mode": "push",
        "primer_set": "B2",
        "trim_primers": "TRUE",
        "sequence_number": 12,
        "pipeline": "classic_luka",
        "luka_light_protocol": "standard",
        "run_igblast": True,
        "run_clonality": True,
        "archive_format": "zip",
        "igblast_species": "mouse",
        "igblast_panel": "ig_all",
        "mibs_junction_max_distance": 4,
        "mibs_gene_grouping": "strict",
        "mibs_clone_scope": "sample",
        "mibs_min_reads_per_well": 6,
        "mibs_min_wells_per_clone": 2,
        "mibs_productive_only": False,
        "mibs_min_reads_per_cluster": 2,
        "platemap_attached": False,
        "sequencing_kit": "miseq_nano_2x250",
        "receptor_check": "auto",
        "igblast_bin": "igblastn",
        "igblast_organism": "mouse",
        "igblast_db_v": None,
        "igblast_db_d": None,
        "igblast_db_j": None,
        "igblast_data": "/opt/tbrc/igblast/internal_data",
        "igblast_auxiliary_data": "/opt/tbrc/igblast/optional_file",
        "trim_seq": "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG",
        "keep_raw": "FALSE",
        "archive_scope": "light",
        "run_group": "",
        "source": "u/run1",
        "dest": "results/u/run1",
    }
    row.update(changes)
    return row


def with_user(user, run="run1", **changes):
    return good_row(user_id=user, source="{}/{}".format(user, run),
                    dest="results/{}/{}".format(user, run), **changes)


def local_row(**changes):
    """A CLI run: nothing is pushed, there is no trusted entry."""
    base = dict(server_name="local", server_path="", server_results_path="", server_ssh_port=22,
                server_transfer_mode="none")
    base.update(changes)
    return good_row(**base)


NAN = float("nan")

# Letters and digits that a \w or \d in a pattern would let through, in every kind of field.
NON_ASCII = ["café", "１２"]
SNAKEFILE = os.path.join(WORKFLOW, "Snakefile")


class SampleRulesCase(unittest.TestCase):
    def problems(self, row, trusted=TRUSTED):
        """The fields refused, or [] when the row passes."""
        try:
            validate_sample(row, trusted)
        except SampleRejected as err:
            return [p.split(":", 1)[0] for p in err.problems]
        return []

    def assertAccepted(self, row, trusted=TRUSTED, msg=None):
        self.assertEqual(self.problems(row, trusted), [], msg)

    def assertRefused(self, row, field, trusted=TRUSTED, msg=None):
        found = self.problems(row, trusted)
        self.assertIn(field, found, msg or "expected {} to be refused, got {}".format(field, found))


class TestBaseline(SampleRulesCase):
    def test_a_good_row_passes_with_and_without_the_trusted_list(self):
        self.assertAccepted(good_row())
        self.assertAccepted(good_row(), trusted=None)

    def test_rejected_is_a_value_error_with_a_problem_per_field(self):
        with self.assertRaises(SampleRejected) as caught:
            validate_sample(good_row(user_id="a b", pipeline="nope"), TRUSTED)
        self.assertIsInstance(caught.exception, ValueError)
        self.assertIsInstance(caught.exception.problems, list)
        text = str(caught.exception)
        for problem in caught.exception.problems:
            self.assertIn(problem, text)

    def test_messages_never_repeat_the_value(self):
        marker = "ZZ hostile;marker"
        for field in ("user_id", "pipeline", "primer_set", "igblast_bin", "server_results_path", "dest"):
            with self.assertRaises(SampleRejected) as caught:
                validate_sample(good_row(**{field: marker}), TRUSTED)
            self.assertNotIn(marker, str(caught.exception), field)

    def test_a_row_that_is_not_a_dict_is_refused(self):
        for bad in (None, [], "x", 3):
            with self.assertRaises(SampleRejected):
                validate_sample(bad, TRUSTED)


class TestSharedCases(SampleRulesCase):
    """tests/validation_cases.json holds the cases the form's tests use too."""

    SEGMENT_FIELDS = {"user_id": False, "primer_set": False,   # field: blank allowed
                      "sample_name": True, "run_group": True, "sequencing_kit": True}

    def row_with(self, field, value):
        if field == "user_id":
            return with_user(value)
        return good_row(**{field: value})

    def test_shared_cases(self):
        for trusted in (TRUSTED, None):
            for field, blank_ok in sorted(self.SEGMENT_FIELDS.items()):
                for good in CASES["path_segment"]["good"]:
                    self.assertAccepted(self.row_with(field, good), trusted, "{} <- {!r}".format(field, good))
                for bad in CASES["path_segment"]["bad"] + CASES["path_segment"]["bad_untrimmed"]:
                    if bad == "" and blank_ok:
                        continue
                    self.assertRefused(self.row_with(field, bad), field, trusted, "{} <- {!r}".format(field, bad))

            for good in CASES["trim_sequence"]["good"]:
                self.assertAccepted(good_row(trim_seq=good), trusted, repr(good))
            for bad in CASES["trim_sequence"]["bad"]:
                self.assertRefused(good_row(trim_seq=bad), "trim_seq", trusted, repr(bad))

            for good in CASES["method"]["good"]:
                self.assertAccepted(good_row(pipeline=good), trusted, repr(good))
            for bad in CASES["method"]["bad"]:
                self.assertRefused(good_row(pipeline=bad), "pipeline", trusted, repr(bad))

            for good in CASES["email"]["good"]:
                self.assertAccepted(good_row(user_email=good), trusted, repr(good))
            for bad in CASES["email"]["bad"]:
                self.assertRefused(good_row(user_email=bad), "user_email", trusted, repr(bad))

    def test_a_hostile_value_is_refused_in_every_field(self):
        values = CASES["hostile"] + CASES["path_segment"]["bad_untrimmed"] + NON_ASCII
        for trusted in (TRUSTED, None):
            for field in FIELD_RULES:
                for value in values:
                    row = with_user(value) if field == "user_id" else good_row(**{field: value})
                    self.assertRefused(row, field, trusted, "{} <- {!r} (trusted={})".format(
                        field, value, trusted is not None))

    def test_a_hostile_value_is_refused_in_every_field_of_a_local_run(self):
        for field in FIELD_RULES:
            for value in CASES["hostile"] + NON_ASCII:
                row = with_user(value, server_name="local", server_path="", server_results_path="",
                                server_ssh_port=22, server_transfer_mode="none") if field == "user_id" \
                    else local_row(**{field: value})
                self.assertRefused(row, field, TRUSTED, "{} <- {!r}".format(field, value))

    def test_trailing_newline_refused(self):
        self.assertRefused(with_user("tiago\n"), "user_id")
        self.assertRefused(good_row(primer_set="B2\n"), "primer_set")
        self.assertRefused(good_row(pipeline="classic_luka\n"), "pipeline")
        self.assertRefused(good_row(trim_seq="GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG\n"), "trim_seq")
        self.assertRefused(good_row(user_email="u@example.org\n"), "user_email")
        self.assertRefused(good_row(server_results_path="me@boxA:/vol/res/\n"), "server_results_path")
        self.assertRefused(good_row(igblast_bin="igblastn\n"), "igblast_bin")
        self.assertRefused(good_row(archive_format="zip\n"), "archive_format")

    def test_no_trimming_or_case_folding(self):
        # the form cleans these before it writes; a value is right as written or refused
        self.assertRefused(good_row(trim_seq="gggaattcgaggtgcagctgcaggagtctgg"), "trim_seq")
        self.assertRefused(good_row(trim_seq=" GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"), "trim_seq")
        self.assertRefused(with_user(" u"), "user_id")
        self.assertRefused(good_row(pipeline="Classic_Luka"), "pipeline")
        self.assertRefused(good_row(archive_format="ZIP"), "archive_format")
        self.assertRefused(good_row(igblast_species="Mouse"), "igblast_species")


def snakefile_columns():
    """Every sample.json column the workflow reads, from the source, guarded reads included."""
    patterns = [
        r"samples_df\['([A-Za-z_]+)'\]",
        r"\bdf\['([A-Za-z_]+)'\]",
        r"'([A-Za-z_]+)' (?:not )?in (?:samples_df|df)\.columns",
        r"_sulk_\w+\((?:samples_df|df), '([A-Za-z_]+)'",
    ]
    columns = set()
    paths = [os.path.join(WORKFLOW, "Snakefile")] + sorted(glob.glob(os.path.join(WORKFLOW, "rules", "*.smk")))
    for path in paths:
        with open(path) as handle:
            for line in handle:
                code = line.split("#", 1)[0]
                for pattern in patterns:
                    columns.update(re.findall(pattern, code))
    return columns


class TestEveryColumnHasARule(SampleRulesCase):
    def test_every_snakefile_column_has_a_rule(self):
        columns = snakefile_columns()
        # the scan itself must keep working, or the check below proves nothing
        for known in ("user_id", "dest", "server_results_path", "run_group", "mibs_clone_scope",
                      "platemap_attached", "igblast_data", "sample_name", "trim_seq"):
            self.assertIn(known, columns)
        self.assertGreaterEqual(len(columns), 35)
        self.assertEqual(sorted(columns - set(FIELD_RULES)), [])

    def test_every_form_and_cli_field_has_a_rule(self):
        for name in ("user_email", "input_method", "server_path", "sequencing_kit", "igblast_organism",
                     "mibs_productive_only", "mibs_min_reads_per_cluster", "source"):
            self.assertIn(name, FIELD_RULES)

    def test_a_missing_column_the_snakefile_indexes_is_reported(self):
        row = good_row()
        del row["trim_seq"]
        del row["igblast_bin"]
        found = self.problems(row)
        self.assertIn("trim_seq", found)
        self.assertIn("igblast_bin", found)

    def test_columns_the_snakefile_guards_may_be_missing(self):
        row = good_row()
        for name in ("pipeline", "sample_name", "archive_format", "archive_scope", "run_group",
                     "receptor_check", "igblast_data", "platemap_attached", "mibs_clone_scope",
                     "mibs_min_reads_per_well", "sequencing_kit", "input_method"):
            del row[name]
        self.assertAccepted(row)


class TestDestination(SampleRulesCase):
    def test_results_path_must_match_trusted(self):
        for bad in ("user@storage1:/path/to/storage", "evil@outside.example:/x/",
                    "me@boxA:/vol/res", "me@boxA:/vol/res/../x/", "me@boxB:/data/res/"):
            self.assertRefused(good_row(server_results_path=bad), "server_results_path", msg=bad)
        self.assertAccepted(good_row(server_results_path="me@boxA:/vol/res/"))
        self.assertAccepted(good_row(server_name="boxB", server_results_path="you@boxB:/data/res/",
                                     server_ssh_port=22))

    def test_results_path_is_a_trusted_entry_of_its_own_server(self):
        # boxB's path with boxA's name is a well-formed path to the wrong folder
        self.assertRefused(good_row(server_results_path="you@boxB:/data/res/"), "server_results_path")

    def test_port_must_match_trusted_for_push(self):
        for bad in (22, 2223, 0, "2222; touch x", True, NAN, None, ""):
            self.assertRefused(good_row(server_ssh_port=bad), "server_ssh_port", msg=repr(bad))
        self.assertAccepted(good_row(server_ssh_port=2222))
        # compared as integers: a float a table may hand over, or a plain digit string
        self.assertAccepted(good_row(server_ssh_port=2222.0))
        self.assertAccepted(good_row(server_ssh_port="2222"))
        self.assertRefused(good_row(server_ssh_port=2222.5), "server_ssh_port")

    def test_mode_must_match_trusted_for_push(self):
        self.assertAccepted(good_row(server_transfer_mode="push"))
        for bad in ("rm -rf", "Push", "", None, NAN, True):
            self.assertRefused(good_row(server_transfer_mode=bad), "server_transfer_mode", msg=repr(bad))

    def test_mode_push_is_refused_for_a_server_the_list_says_only_pulls(self):
        trusted = copy.deepcopy(TRUSTED)
        trusted["servers"]["boxP"] = {"rawdata_path": "p@boxP:/raw/", "results_path": "p@boxP:/res/",
                                      "ssh_port": 22, "transfer_mode": "pull"}
        row = good_row(server_name="boxP", server_results_path="p@boxP:/res/", server_ssh_port=22)
        self.assertRefused(dict(row, server_transfer_mode="push"), "server_transfer_mode", trusted)
        self.assertAccepted(dict(row, server_transfer_mode="pull"), trusted)

    def test_empty_destination_only_with_pull_or_none(self):
        for empty in ("", None, NAN, "NA"):
            self.assertRefused(good_row(server_results_path=empty), "server_results_path", msg=repr(empty))
            for mode in ("pull", "none"):
                self.assertAccepted(good_row(server_results_path=empty, server_transfer_mode=mode),
                                    msg="{!r} with {}".format(empty, mode))

    def test_nonempty_destination_with_none_must_be_trusted(self):
        # final.script.sh pushes unless the mode is "pull", so "none" is no protection
        for mode in ("pull", "none"):
            self.assertAccepted(good_row(server_transfer_mode=mode))
            self.assertRefused(good_row(server_transfer_mode=mode,
                                        server_results_path="evil@outside.example:/x/"),
                               "server_results_path", msg=mode)
            self.assertRefused(good_row(server_transfer_mode=mode,
                                        server_results_path="user@storage1:/path/to/storage"),
                               "server_results_path", msg=mode)

    def test_port_must_match_trusted_whenever_the_mode_is_not_pull_and_the_destination_is_set(self):
        # final.script.sh pushes for every mode but pull once the destination is not empty, so with
        # mode "none" the trusted path is pushed to on whatever port the file names
        for mode in ("push", "none"):
            self.assertAccepted(good_row(server_transfer_mode=mode, server_ssh_port=2222), msg=mode)
            for bad in (1, 22, 2223, 65535, "", None, NAN, "NA", "2222; x", True, 2222.5):
                self.assertRefused(good_row(server_transfer_mode=mode, server_ssh_port=bad),
                                   "server_ssh_port", msg="{} {!r}".format(mode, bad))
        # the same holds for a mode that is none of the three, which pushes the same way
        self.assertRefused(good_row(server_transfer_mode="rm", server_ssh_port=1), "server_ssh_port")
        # the other trusted server has its own port
        boxb = dict(server_name="boxB", server_results_path="you@boxB:/data/res/", server_transfer_mode="none")
        self.assertAccepted(good_row(server_ssh_port=22, **boxb))
        for bad in (2222, 1, "", "NA"):
            self.assertRefused(good_row(server_ssh_port=bad, **boxb), "server_ssh_port", msg=repr(bad))

    def test_port_is_free_when_nothing_can_be_pushed(self):
        # pull leaves the archive staged; mode none with an empty destination has nowhere to push
        for bad_port in (1, 22, 65535, "", None, NAN, "NA"):
            self.assertAccepted(good_row(server_transfer_mode="pull", server_ssh_port=bad_port),
                                msg="pull {!r}".format(bad_port))
            for empty in ("", None, NAN, "NA"):
                self.assertAccepted(good_row(server_transfer_mode="none", server_results_path=empty,
                                             server_ssh_port=bad_port),
                                    msg="none {!r} {!r}".format(empty, bad_port))
        # but it is still a port number when it is there
        for bad in (0, 65536, "22; x", "２２２２", True):
            self.assertRefused(good_row(server_transfer_mode="pull", server_ssh_port=bad), "server_ssh_port",
                               msg=repr(bad))
            self.assertRefused(good_row(server_transfer_mode="none", server_results_path="",
                                        server_ssh_port=bad), "server_ssh_port", msg=repr(bad))

    def test_port_digits_are_ascii(self):
        # int("２２２２") is 2222, so a pattern that lets full-width digits through would pass the
        # trusted port
        for trusted in (TRUSTED, None):
            for bad in ("２２２２", "２２", "２２２２ ", "²²²²"):
                self.assertRefused(good_row(server_ssh_port=bad), "server_ssh_port", trusted, repr(bad))
                self.assertRefused(good_row(server_ssh_port=bad, server_transfer_mode="pull"),
                                   "server_ssh_port", trusted, repr(bad))
        self.assertRefused(local_row(server_ssh_port="２２"), "server_ssh_port")

    def test_pull_still_takes_the_port_as_a_number_or_nothing(self):
        self.assertAccepted(good_row(server_transfer_mode="pull", server_ssh_port=22))
        self.assertAccepted(good_row(server_transfer_mode="pull", server_ssh_port=""))
        self.assertRefused(good_row(server_transfer_mode="pull", server_ssh_port="22; x"), "server_ssh_port")

    def test_unknown_server_is_refused(self):
        for bad in ("storage1", "boxC", None, NAN, "BoxA", "boxA ", "local ", "Local", "upload\n"):
            self.assertRefused(good_row(server_name=bad), "server_name", msg=repr(bad))

    def test_empty_server_name_is_treated_like_local(self):
        # `tbrc run --fastq-dir` writes "" when there is no storage server
        self.assertAccepted(local_row(server_name=""))
        for empty in ("", None, NAN, "NA"):
            self.assertAccepted(local_row(server_name="", server_results_path=empty), msg=repr(empty))
        for mode in ("push", "pull", "none"):
            self.assertAccepted(local_row(server_name="", server_transfer_mode=mode), msg=mode)
        for port in (22, "", None, NAN):
            self.assertAccepted(local_row(server_name="", server_ssh_port=port), msg=repr(port))
        # a destination is refused whatever the mode, including the trusted server's own path
        for path in ("me@boxA:/vol/res/", "evil@outside.example:/x/", "/tmp/out/", "x"):
            for mode in ("push", "pull", "none"):
                self.assertRefused(local_row(server_name="", server_results_path=path, server_transfer_mode=mode),
                                   "server_results_path", msg="{} {}".format(path, mode))
        self.assertRefused(local_row(server_name="", server_transfer_mode="rm"), "server_transfer_mode")
        self.assertRefused(local_row(server_name="", server_ssh_port="22; x"), "server_ssh_port")
        # only the empty string, not a missing value or one with spaces around it
        for bad in (None, NAN, " ", "\n"):
            self.assertRefused(local_row(server_name=bad), "server_name", msg=repr(bad))

    def test_local_server_needs_empty_destination(self):
        self.assertAccepted(local_row())
        for empty in ("", None, NAN, "NA"):
            self.assertAccepted(local_row(server_results_path=empty), msg=repr(empty))
        for port in (22, "", None, NAN, 2112):
            self.assertAccepted(local_row(server_ssh_port=port), msg=repr(port))
        for mode in ("push", "pull", "none"):
            self.assertAccepted(local_row(server_transfer_mode=mode), msg=mode)
        for path in ("me@boxA:/vol/res/", "evil@outside.example:/x/", "/tmp/out/", "x"):
            for mode in ("push", "pull", "none"):
                self.assertRefused(local_row(server_results_path=path, server_transfer_mode=mode),
                                   "server_results_path", msg="{} {}".format(path, mode))
        self.assertRefused(local_row(server_transfer_mode="rm"), "server_transfer_mode")
        self.assertRefused(local_row(server_ssh_port="22; x"), "server_ssh_port")

    def test_upload_server_is_treated_like_local(self):
        row = good_row(server_name="upload", input_method="Upload", server_path=NAN,
                       server_results_path=NAN, server_ssh_port=NAN, server_transfer_mode="push")
        self.assertAccepted(row)
        self.assertRefused(good_row(server_name="upload", server_results_path="me@boxA:/vol/res/"),
                           "server_results_path")

    def test_dest_must_match_user(self):
        self.assertAccepted(good_row(dest="results/u/run1"))
        for bad in ("results/other/run", "results/u/../x", "results/u/a/b", "results/u/", "results/u",
                    "results//run", "/results/u/run", "results/u/run/", "other/u/run", "results/u/.",
                    "results/u/..", "results/u/-x", "results/u/a b", "results/u/run\n", "", None, NAN):
            self.assertRefused(good_row(dest=bad), "dest", msg=repr(bad))

    def test_source_is_user_and_run_or_one_segment(self):
        for good in ("u/run1", "run1", "RUN15", "231003_M03851_0285_000000000-DLPWD"):
            self.assertAccepted(good_row(source=good), msg=good)
        for bad in ("other/run1", "u/../x", "u/a/b", "u/", "/u/run1", "", "..", "-x", "u/run1\n", None, NAN):
            self.assertRefused(good_row(source=bad), "source", msg=repr(bad))

    def test_server_path_is_format_checked(self):
        self.assertAccepted(good_row(server_path="me@boxA:/vol/raw/"))
        for empty in ("", None, NAN, "NA"):
            self.assertAccepted(good_row(server_path=empty), msg=repr(empty))
        for bad in ("a b", "me@boxA:/vol/raw/;x", "$(x)", "../x", "-x"):
            self.assertRefused(good_row(server_path=bad), "server_path", msg=bad)


class TestIgblast(SampleRulesCase):
    def test_igblast_bin_must_match(self):
        self.assertAccepted(good_row(igblast_bin="igblastn"))
        for bad in ("/tmp/evil", "igblastn ", "igblastx", "", None, NAN, "igblastn;x"):
            self.assertRefused(good_row(igblast_bin=bad), "igblast_bin", msg=repr(bad))

    def test_legacy_igblast_dirs_accepted(self):
        for data in TRUSTED["igblast"]["igblast_data"]:
            self.assertAccepted(good_row(igblast_data=data), msg=data)
        for aux in TRUSTED["igblast"]["igblast_auxiliary_data"]:
            self.assertAccepted(good_row(igblast_auxiliary_data=aux), msg=aux)

    def test_igblast_paths_must_be_on_the_list(self):
        for bad in ("/tmp/evil/internal_data", "/opt/tbrc/igblast/internal_data/", "/opt/tbrc/igblast/optional_file"):
            self.assertRefused(good_row(igblast_data=bad), "igblast_data", msg=bad)
        self.assertRefused(good_row(igblast_auxiliary_data="/opt/tbrc/igblast/internal_data"),
                           "igblast_auxiliary_data")
        # a list holds the strings that are accepted, none of them a prefix match
        self.assertRefused(good_row(igblast_data="/opt/tbrc/igblast"), "igblast_data")

    def test_igblast_data_and_aux_empty_only_if_the_list_holds_empty(self):
        for empty in (None, NAN, ""):
            self.assertRefused(good_row(igblast_data=empty), "igblast_data", msg=repr(empty))
            self.assertRefused(good_row(igblast_auxiliary_data=empty), "igblast_auxiliary_data", msg=repr(empty))
        allowing = copy.deepcopy(TRUSTED)
        allowing["igblast"]["igblast_data"].append("")
        for empty in (None, NAN, ""):
            self.assertAccepted(good_row(igblast_data=empty), allowing, msg=repr(empty))

    def test_igblast_db_prefix(self):
        for empty in (None, NAN, ""):
            for field in ("igblast_db_v", "igblast_db_d", "igblast_db_j"):
                self.assertAccepted(good_row(**{field: empty}), msg="{} {!r}".format(field, empty))
        for field in ("igblast_db_v", "igblast_db_d", "igblast_db_j"):
            self.assertRefused(good_row(**{field: "/tmp/evil/db"}), field)
        strict = copy.deepcopy(TRUSTED)
        strict["igblast"]["igblast_db_v"] = ["/opt/db/v"]
        self.assertAccepted(good_row(igblast_db_v="/opt/db/v"), strict)
        for empty in (None, NAN, ""):
            self.assertRefused(good_row(igblast_db_v=empty), "igblast_db_v", strict, repr(empty))

    def test_igblast_panel_takes_every_choice_and_the_legacy_one(self):
        for panel in ("igh", "ig_all", "tra", "trb", "tcr_all", "all_receptors", "ig_b6", "ig_light"):
            self.assertAccepted(good_row(igblast_panel=panel), msg=panel)
        for bad in ("ig_light ", "IG_ALL", "all", "", None, NAN, "ig_all;x"):
            self.assertRefused(good_row(igblast_panel=bad), "igblast_panel", msg=repr(bad))

    def test_species_and_organism(self):
        for field in ("igblast_species", "igblast_organism"):
            for good in ("mouse", "human"):
                self.assertAccepted(good_row(**{field: good}))
            for bad in ("rat", "", None, NAN):
                self.assertRefused(good_row(**{field: bad}), field, msg=repr(bad))


class TestValueShapes(SampleRulesCase):
    def test_nan_and_null_for_optional_pass(self):
        for blank in (NAN, None):
            for field in ("run_group", "igblast_db_v", "sample_name", "archive_scope", "receptor_check",
                          "archive_format", "sequencing_kit", "platemap_attached", "mibs_gene_grouping",
                          "mibs_clone_scope", "mibs_junction_max_distance", "mibs_productive_only", "pipeline"):
                self.assertAccepted(good_row(**{field: blank}), msg="{} {!r}".format(field, blank))
        for field in ("run_group", "sample_name", "archive_scope", "receptor_check", "archive_format",
                      "sequencing_kit", "mibs_clone_scope", "platemap_attached"):
            self.assertAccepted(good_row(**{field: ""}), msg="{} empty".format(field))
        self.assertAccepted(good_row(user_email=""))
        self.assertAccepted(good_row(user_email=NAN))
        self.assertAccepted(good_row(user_email=None))

    def test_the_numeric_na_string_jsonlite_writes_is_empty_for_the_optional_numbers(self):
        for field in ("mibs_junction_max_distance", "mibs_min_reads_per_well", "mibs_min_wells_per_clone",
                      "mibs_min_reads_per_cluster"):
            self.assertAccepted(good_row(**{field: "NA"}), msg=field)
        self.assertRefused(good_row(sequence_number="NA"), "sequence_number")

    def test_nan_for_required_is_reported(self):
        self.assertRefused(good_row(trim_seq=NAN, trim_primers="TRUE"), "trim_seq")
        for field in ("user_id", "primer_set", "trim_primers", "sequence_number", "luka_light_protocol",
                      "run_igblast", "run_clonality", "keep_raw", "igblast_species", "igblast_panel",
                      "igblast_organism", "igblast_bin", "server_name", "dest", "source"):
            # server_name "" is a local run (CLI, no storage server); a missing value is still refused
            for blank in ((NAN, None) if field == "server_name" else (NAN, None, "")):
                row = with_user(blank) if field == "user_id" else good_row(**{field: blank})
                self.assertRefused(row, field, msg="{} {!r}".format(field, blank))

    def test_trim_seq_is_empty_only_when_trimming_is_off(self):
        for blank in (NAN, None, ""):
            self.assertAccepted(good_row(trim_primers="FALSE", trim_seq=blank), msg=repr(blank))
            self.assertAccepted(good_row(trim_primers=False, trim_seq=blank), msg=repr(blank))
            self.assertRefused(good_row(trim_primers="TRUE", trim_seq=blank), "trim_seq", msg=repr(blank))
            self.assertRefused(good_row(trim_primers=True, trim_seq=blank), "trim_seq", msg=repr(blank))
        # a sequence that is there is checked whether or not the trim is on
        self.assertRefused(good_row(trim_primers="FALSE", trim_seq="ACGT;touch x"), "trim_seq")
        self.assertAccepted(good_row(trim_primers="FALSE"))

    def test_legacy_boolean_strings_accepted(self):
        for field in ("run_igblast", "run_clonality", "keep_raw", "trim_primers", "platemap_attached",
                      "mibs_productive_only"):
            for value in (True, False, "TRUE", "FALSE", "true", "false"):
                self.assertAccepted(good_row(**{field: value}), msg="{} {!r}".format(field, value))
            for bad in ("True", "yes", "1", 1, 0, "T", "TRUE ", "TRUE\n", "on"):
                self.assertRefused(good_row(**{field: bad}), field, msg="{} {!r}".format(field, bad))

    def test_sequence_number(self):
        for good in (1, 6, 12, 200, 1000, 12.0):
            self.assertAccepted(good_row(sequence_number=good), msg=repr(good))
        for bad in (0, -2, 1001, 12.5, "12", "12; x", True, [12], {"a": 1}, NAN, None):
            self.assertRefused(good_row(sequence_number=bad), "sequence_number", msg=repr(bad))

    def test_mibs_ranges(self):
        for field, (low, high) in (("mibs_min_reads_per_well", (1, 10000)),
                                   ("mibs_min_wells_per_clone", (1, 96)),
                                   ("mibs_junction_max_distance", (0, 30)),
                                   ("mibs_min_reads_per_cluster", (1, 1000))):
            self.assertAccepted(good_row(**{field: low}))
            self.assertAccepted(good_row(**{field: high}))
            self.assertRefused(good_row(**{field: low - 1}), field)
            self.assertRefused(good_row(**{field: high + 1}), field)
            self.assertRefused(good_row(**{field: "4"}), field)
            self.assertRefused(good_row(**{field: True}), field)

    def test_choices(self):
        for field, good, bad in (
                ("input_method", ("Server", "Upload"), ("server", "upload", "Server ")),
                ("luka_light_protocol", ("standard", "luka_light"), ("Standard", "light")),
                ("archive_format", ("zip", "tar.gz"), ("tgz", "targz", "TAR.GZ")),
                ("archive_scope", ("light", "full"), ("auto", "Light")),
                ("receptor_check", ("auto", "warn", "off"), ("on", "AUTO")),
                ("mibs_gene_grouping", ("strict", "transitive"), ("loose",)),
                ("mibs_clone_scope", ("sample", "global"), ("run",))):
            for value in good:
                self.assertAccepted(good_row(**{field: value}), msg="{} {}".format(field, value))
            for value in bad:
                self.assertRefused(good_row(**{field: value}), field, msg="{} {}".format(field, value))

    def test_wrong_types_are_refused_in_every_field(self):
        for trusted in (TRUSTED, None):
            for field in FIELD_RULES:
                for value in ([], {}, [1], {"a": "b"}, object(), b"x"):
                    self.assertRefused(good_row(**{field: value}), field, trusted,
                                       "{} <- {} (trusted={})".format(field, type(value).__name__,
                                                                      trusted is not None))
                    self.assertRefused(local_row(**{field: value}), field, trusted,
                                       "{} <- {} (local, trusted={})".format(field, type(value).__name__,
                                                                             trusted is not None))

    def test_odd_scalars_are_refused_or_accepted_but_never_raise(self):
        for field in FIELD_RULES:
            for value in (1.5, -1, True, 0, float("inf"), -0.0):
                try:
                    validate_sample(good_row(**{field: value}), TRUSTED)
                except SampleRejected:
                    pass

    def test_email(self):
        for good in ("a@b.co", "first.last+tag@lab.example.org", "x_y%z@a-b.example.com"):
            self.assertAccepted(good_row(user_email=good), msg=good)
        for bad in ("a b@c.co", "a@b", "a@b.c", "'a@b.co", "a@b.co;x", "$(x)@b.co", "a@b.co\nBcc: x@y.co",
                    "a" * 250 + "@b.co", 5):
            self.assertRefused(good_row(user_email=bad), "user_email", msg=repr(bad))
        # non-ASCII letters and digits in the local part, the host and the last label
        for bad in ("café@example.org", "１２@example.org", "a@exämple.org", "a@ex１２.org", "a@example.cö",
                    "a@example.１２", "a@example.ｏｒｇ"):
            self.assertRefused(good_row(user_email=bad), "user_email", msg=repr(bad))

    def test_all_problems_are_reported_together(self):
        row = good_row(user_id="bad id", pipeline="nope", trim_seq="ACGT;x", server_results_path="evil@x:/y/",
                       igblast_bin="/tmp/evil", run_igblast="maybe")
        with self.assertRaises(SampleRejected) as caught:
            validate_sample(row, TRUSTED)
        fields = [p.split(":", 1)[0] for p in caught.exception.problems]
        for expected in ("user_id", "pipeline", "trim_seq", "server_results_path", "igblast_bin", "run_igblast"):
            self.assertIn(expected, fields)
        self.assertGreaterEqual(len(caught.exception.problems), 6)


class TestNoTrustedFile(SampleRulesCase):
    def test_no_trusted_file_means_format_checks_only(self):
        # a well-formed destination to any host passes, which is the point of "until its turn comes"
        self.assertAccepted(good_row(server_name="storage1", server_results_path="user@storage1:/a/b/",
                                     server_ssh_port=2112), trusted=None)
        self.assertAccepted(good_row(igblast_bin="/usr/local/bin/igblastn", igblast_data="/x/internal_data",
                                     igblast_auxiliary_data="/x/optional_file", igblast_db_v="/x/db/v"),
                            trusted=None)
        self.assertAccepted(good_row(igblast_data=None, igblast_auxiliary_data=NAN), trusted=None)
        self.assertAccepted(local_row(), trusted=None)

    def test_format_still_applies_without_it(self):
        for bad in ("evil host:/x/", "a@b:/x/;touch x", "a@b:/../x/", "-oProxyCommand=x", "a@b:/x/ -e x",
                    "$(x)", "`x`"):
            self.assertRefused(good_row(server_results_path=bad), "server_results_path", None, bad)
        for bad in ("/tmp/evil bin", "-x", "../x", "x;y"):
            self.assertRefused(good_row(igblast_bin=bad), "igblast_bin", None, bad)
        for bad in (0, 65536, "22; x", True):
            self.assertRefused(good_row(server_ssh_port=bad), "server_ssh_port", None, repr(bad))
        self.assertRefused(good_row(server_transfer_mode="rm"), "server_transfer_mode", None)
        self.assertRefused(good_row(server_name="a b"), "server_name", None)
        self.assertRefused(good_row(server_name=None), "server_name", None)
        # "" is a local run here too, so it takes no destination
        self.assertAccepted(local_row(server_name=""), trusted=None)
        self.assertRefused(good_row(server_name=""), "server_results_path", None)

    def test_local_with_a_destination_is_refused_even_without_it(self):
        self.assertRefused(local_row(server_results_path="a@b:/x/"), "server_results_path", None)


class TestSampleNameDerivedFromReads(SampleRulesCase):
    """With no sample_name in sample.json the Snakefile takes the name from the FASTQ file names, which
    come from disk, not from the file the row check read; that name is checked on its own."""

    def test_a_plain_name_passes(self):
        for good in CASES["path_segment"]["good"] + ["RG-78_S1", "RUN15"]:
            self.assertIsNone(validate_sample_name(good), good)

    def test_a_name_that_is_not_one_path_segment_is_refused(self):
        bad_names = (CASES["path_segment"]["bad"] + CASES["path_segment"]["bad_untrimmed"] + CASES["hostile"]
                     + NON_ASCII + ["a;b", "tiago\n", "a\n", " a", "a b", "", "-x", ".", ".."])
        for bad in bad_names:
            with self.assertRaises(SampleRejected, msg=repr(bad)):
                validate_sample_name(bad)

    def test_something_that_is_not_text_is_refused(self):
        for bad in (None, NAN, 5, b"x", [], {}, object()):
            with self.assertRaises(SampleRejected, msg=repr(bad)):
                validate_sample_name(bad)

    def test_the_message_names_the_field_and_not_the_value(self):
        for bad in ("a;b$(touch canary)", "tiago\n", "café"):
            with self.assertRaises(SampleRejected) as caught:
                validate_sample_name(bad)
            self.assertIn("sample_name", str(caught.exception))
            self.assertEqual(len(caught.exception.problems), 1)
            for piece in ("touch", "canary", "tiago", "café", "a;b"):
                self.assertNotIn(piece, str(caught.exception))


class TestRowsBelongToTheirRunFolder(SampleRulesCase):
    """A sample.json that sits in data/runs/<user>/<run>/ may only describe that run: dest drives the
    cleanup in final.script.sh, and the web app's key can put any file in any run folder."""

    TABLE = "data/runs/u/run1/sample.json"

    def refused(self, rows, table=TABLE):
        """The fields refused, or [] when the rows pass."""
        try:
            validate_run_folder(rows, table)
        except SampleRejected as err:
            return [p.split(":", 1)[0] for p in err.problems]
        return []

    def test_a_row_that_matches_its_folder_passes(self):
        self.assertIsNone(validate_run_folder([good_row()], self.TABLE))

    def test_every_path_shape_the_pipeline_receives_is_matched(self):
        # relative inside .parallel/<key> (data/ is a symlink to the root's), from the root,
        # absolute on the cluster or in a container, with a legacy table name, and with odd slashes
        for table in ("data/runs/u/run1/sample.json", "./data/runs/u/run1/sample.json",
                      "data//runs/u/run1/sample.json", "data/runs/u/run1//sample.json",
                      "/srv/pipeline/TBRC/data/runs/u/run1/sample.json", "/work/data/runs/u/run1/sample.json",
                      ".parallel/u-run1/data/runs/u/run1/sample.json", "data/runs/u/run1/sample.tsv"):
            self.assertIsNone(validate_run_folder([good_row()], table), table)
            self.assertEqual(self.refused([good_row(user_id="v")], table), ["user_id"], table)

    def test_a_row_for_another_user_is_refused(self):
        self.assertEqual(self.refused([good_row(user_id="v")]), ["user_id"])
        # all three describe v, none of them the folder's u
        self.assertEqual(self.refused([with_user("v")]), ["user_id", "dest", "source"])

    def test_a_dest_for_another_run_is_refused(self):
        # dest is exact: it drives the cleanup, and unlike source it has no older one-segment form
        for dest in ("run1", "results/u/run2", "results/u/run1/", "results/u/run1/x", "results/u", "results/u/",
                     "results/v/run1", "results//u/run1", "results/u/../run1", "u/run1", "/results/u/run1",
                     "results/u/run1 ", "results/u/RUN1", "results/u/run1\n", ""):
            self.assertEqual(self.refused([good_row(dest=dest)]), ["dest"], repr(dest))

    def test_a_source_may_be_user_slash_run_or_the_older_run_alone(self):
        # Older sample.json files, written before the form added the user, hold just the run; no rule
        # reads source (SOUR in the Snakefile is not used), so re-running such a run must keep working.
        for table in ("data/runs/u/run1/sample.json", "/srv/pipeline/TBRC/data/runs/u/run1/sample.json",
                      ".parallel/u-run1/data/runs/u/run1/sample.json"):
            for source in ("u/run1", "run1"):
                self.assertEqual(self.refused([good_row(source=source)], table), [], (table, source))
        # the same two shapes for another folder
        self.assertEqual(self.refused([good_row(user_id="alice", source="plate7", dest="results/alice/plate7")],
                                      "data/runs/alice/plate7/sample.json"), [])
        self.assertEqual(self.refused([with_user("alice", "plate7")], "data/runs/alice/plate7/sample.json"), [])

    def test_an_older_style_row_passes_the_whole_check(self):
        row = good_row(source="run1")
        self.assertAccepted(row)
        self.assertIsNone(validate_run_folder([row], self.TABLE))

    def test_a_source_for_another_run_is_refused(self):
        # another run's name, another user's folder, the user alone, and anything not exactly one of the
        # two accepted shapes
        for source in ("u/run2", "run2", "v/run1", "v/run2", "v", "u", "other/u/run1", "u/run1/", "run1/",
                       "/run1", "u/run1/x", "u/", "/u/run1", "u//run1", "u/run1 ", " run1", "run1 ", "run1\n",
                       "U/run1", "RUN1", "../run1", "u/../run1", "results/u/run1", ""):
            self.assertEqual(self.refused([good_row(source=source)]), ["source"], repr(source))

    def test_values_that_are_not_text_are_refused_not_raised(self):
        for field in ("user_id", "dest", "source"):
            for bad in (None, NAN, 5, True, b"u", ["u"], {}, object()):
                self.assertEqual(self.refused([good_row(**{field: bad})]), [field], "{} {!r}".format(field, bad))

    def test_a_missing_field_is_refused(self):
        for field in ("user_id", "dest", "source"):
            row = good_row()
            del row[field]
            self.assertEqual(self.refused([row]), [field], field)

    def test_every_row_is_checked(self):
        self.assertEqual(self.refused([good_row(), good_row(dest="results/u/other")]), ["dest"])
        self.assertEqual(self.refused([good_row(dest="results/u/other"), good_row(source="u/other")]),
                         ["dest", "source"])
        self.assertEqual(self.refused([good_row(), good_row()]), [])

    def test_a_row_that_is_not_a_row_is_refused(self):
        for bad in (None, "u", ["u"], 5):
            self.assertTrue(self.refused([bad]), repr(bad))

    def test_a_table_anywhere_else_is_left_alone(self):
        wrong = [good_row(user_id="v", dest="results/x/y", source="z")]
        for table in ("/fixtures/classic_luka.sample.json", "sample.json", "data/runs/sample.json",
                      "data/runs/u/sample.json", "runs/u/run1/sample.json", "data/run/u/run1/sample.json",
                      "other/runs/u/run1/sample.json", "data/runs/u/run1/deeper/sample.json",
                      "data/runs/u/run1", "/tmp/data/results/u/run1/sample.json", "", None):
            self.assertIsNone(validate_run_folder(wrong, table), repr(table))

    def test_the_message_names_the_fields_and_not_the_values(self):
        row = good_row(user_id="mallory", dest="results/mallory/loot", source="mallory/loot")
        with self.assertRaises(SampleRejected) as caught:
            validate_run_folder([row], "data/runs/victim/prize/sample.json")
        text = str(caught.exception) + " ".join(caught.exception.problems)
        for field in ("user_id", "dest", "source"):
            self.assertIn(field, text)
        for value in ("mallory", "loot", "victim", "prize"):
            self.assertNotIn(value, text, value)
        self.assertEqual(len(caught.exception.problems), 3)

    def test_the_message_does_not_echo_awkward_values_either(self):
        for value in ("a;b$(touch canary)", "tiago\n", "caf\u00e9"):
            with self.assertRaises(SampleRejected) as caught:
                validate_run_folder([good_row(user_id=value, dest=value, source=value)], self.TABLE)
            for piece in ("touch", "canary", "tiago", "caf"):
                self.assertNotIn(piece, str(caught.exception))

    def test_the_benchmark_rows_pass(self):
        # benchmark/run_benchmark.sh puts its table at data/runs/bench/<name>/ with user bench, so the
        # rule applies to it; the row and the table path are read from the script so they stay in step
        script = os.path.join(REPO, "benchmark", "run_benchmark.sh")
        if not os.path.exists(script):
            self.skipTest("no benchmark/ in this tree")
        with open(script) as handle:
            text = handle.read()
        name = "D1"
        user = re.search(r'"user_id": "([^"]+)"', text).group(1)
        source = re.search(r'"source": f"([^"]+)"', text).group(1).replace("{name}", name)
        dest = re.search(r'"dest": f"([^"]+)"', text).group(1).replace("{name}", name)
        table = re.search(r'DATA_TABLE="([^"]+)"', text).group(1).replace("${name}", name)
        self.assertEqual(table, "data/runs/bench/D1/sample.json")
        self.assertIsNone(validate_run_folder([good_row(user_id=user, source=source, dest=dest)], table))
        self.assertIsNone(validate_run_folder([good_row(user_id="bench", source="bench/D1", dest="results/bench/D1")],
                                              "data/runs/bench/D1/sample.json"))

    def test_the_tbrc_self_test_rows_match_their_folders(self):
        # `tbrc test` writes each dry-run table to data/runs/tbrc_test/<name>_<pipeline>/ and used to
        # build the row for a run called "selftest", which this rule refuses. No exemption: the row is
        # built for the folder's own name at both call sites (the local and the cluster dry run).
        tbrc_path = os.path.join(REPO, "tbrc")
        if not os.path.exists(tbrc_path):
            self.skipTest("no tbrc in this tree")
        with open(tbrc_path) as handle:
            text = handle.read()
        writes = re.findall(r"write_sample\([^\n]*\)\)\n", text)
        dry_runs = [w for w in writes if "build_sample(cfg, fake," in w]
        self.assertEqual(len(dry_runs), 2, dry_runs)
        for call in dry_runs:
            self.assertIn('build_sample(cfg, fake, name + "_" + pipeline)', call)
        # and the row that call builds agrees with the folder it is written to
        import importlib.machinery
        import importlib.util
        import argparse
        loader = importlib.machinery.SourceFileLoader("tbrc_cli_for_test", tbrc_path)
        spec = importlib.util.spec_from_loader("tbrc_cli_for_test", loader)
        cli = importlib.util.module_from_spec(spec)
        keep = sys.dont_write_bytecode
        sys.dont_write_bytecode = True  # no tbrc .pyc dropped into the repository root
        try:
            loader.exec_module(cli)
        finally:
            sys.dont_write_bytecode = keep
        # the lab's own config when the tree has one; a tree without it (the public export) has the example
        config_path = cli.CONFIG_FILE if os.path.exists(cli.CONFIG_FILE) else EXAMPLE_CONFIG
        if not os.path.exists(config_path):
            self.skipTest("this tree has neither frontend/server_config.json nor server_config.example.json")
        with open(config_path, encoding="utf-8") as handle:
            cfg = json.load(handle)
        name = "_tbrc_selftest_20261109_000000"
        for pipeline in cli.PIPELINES:
            fake = argparse.Namespace(**cli.DEFAULTS, user="tbrc_test", email="", pipeline=pipeline,
                                      primer_set="B2", trim_seq="", platemap="", push=False)
            row = cli.build_sample(cfg, fake, name + "_" + pipeline)
            table = "data/runs/tbrc_test/{}_{}/sample.json".format(name, pipeline)
            self.assertIsNone(validate_run_folder([row], table), pipeline)
            # while the old "selftest" row is refused there (the rule has no special case for it)
            old = cli.build_sample(cfg, fake, "selftest")
            self.assertEqual(self.refused([old], table), ["dest", "source"], pipeline)

    def test_a_row_the_form_writes_passes_both_checks(self):
        row = with_user("alice", "plate7")
        self.assertAccepted(row)
        self.assertIsNone(validate_run_folder([row], "data/runs/alice/plate7/sample.json"))


def order_problems(path):
    """What is wrong with where the Snakefile at `path` runs its checks; [] when it is right.

    The row check has to come before anything reads samples_df, the rows have to be tied to their run
    folder right after it (after the validate_sample loop, before the sample name is resolved or a rule
    file is included), the sample name derived from the FASTQ files has to be checked before a rule file
    is included, and each of those has to be there."""
    with open(path) as handle:
        code = "\n".join(line for line in handle.read().splitlines() if not line.lstrip().startswith("#"))
    problems = []

    def where(needle, start=0):
        return code.find(needle, start)

    first_include = where("include:")
    load = where('load_trusted("config/trusted_servers.json")')
    rows = where("validate_sample(")
    folder = where("validate_run_folder(")
    assign = where("samples_df = load_samples_table(")
    derived = where("SAMPLE = resolve_sample_name(")
    first_resolve = where("resolve_sample_name(samples_df")
    name_check = where("validate_sample_name(SAMPLE)")
    check_name = where("check_sample_name(SAMPLE")

    for what, spot in (('load_trusted("config/trusted_servers.json")', load), ("validate_sample(", rows),
                       ("validate_run_folder(samples_df.to_dict(\"records\"), DATA_TABLE)", where(
                           'validate_run_folder(samples_df.to_dict("records"), DATA_TABLE)')),
                       ("validate_sample_name(SAMPLE)", name_check), ("an include:", first_include),
                       ("SAMPLE = resolve_sample_name(", derived), ("samples_df = load_samples_table(", assign)):
        if spot < 0:
            problems.append("missing: " + what)
    if problems:
        return problems

    if not load < rows:
        problems.append("validate_sample( comes before the trusted list is loaded")
    for what, spot in (("load_trusted", load), ("validate_sample(", rows), ("validate_run_folder(", folder)):
        if spot > first_resolve:
            problems.append(what + " comes after the first resolve_sample_name(samples_df")
        if spot > first_include:
            problems.append(what + " comes after the first include:")
    # the first thing that reads samples_df after it is loaded must be the row check
    first_read = where("samples_df", assign + len("samples_df"))
    if not code.startswith('samples_df.to_dict("records")', first_read):
        problems.append("samples_df is read before the row check")
    elif first_read > rows:
        problems.append("samples_df.to_dict(\"records\") is not what validate_sample( is given")
    if not rows < folder:
        problems.append("validate_run_folder( does not come after the validate_sample( loop")
    if not derived < name_check < check_name:
        problems.append("validate_sample_name(SAMPLE) is not right after SAMPLE = resolve_sample_name(")
    if name_check > first_include:
        problems.append("validate_sample_name(SAMPLE) comes after the first include:")
    return problems


class TestFinalScriptChecksTheRecipient(unittest.TestCase):
    """final.script.sh re-checks the address before it hands it to a mail program, because the
    script should not rely on its caller. This runs the script's own send_notification_email()
    with a stub `mail` on the path and asks which of the shared email cases reach it."""

    SCRIPT = os.path.join(REPO, "backend", "workflow", "scripts", "final.script.sh")

    def function_text(self):
        lines, inside = [], False
        with open(self.SCRIPT, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("send_notification_email() {"):
                    inside = True
                if inside:
                    lines.append(line)
                    if line.rstrip("\n") == "}":
                        break
        self.assertTrue(lines, "send_notification_email() not found in final.script.sh")
        return "".join(lines)

    def run_with_stub_mail(self, recipient):
        """Returns (the arguments the stub mail was called with, or None; stderr text)."""
        with tempfile.TemporaryDirectory() as tmp:
            stub = os.path.join(tmp, "mail")
            record = os.path.join(tmp, "mail.args")
            with open(stub, "w") as handle:
                handle.write('#!/bin/sh\nfor a in "$@"; do printf \'%s\\n\' "$a"; done > "$MAIL_RECORD"\ncat > /dev/null\n')
            os.chmod(stub, 0o755)
            env = {"PATH": tmp + os.pathsep + "/usr/bin" + os.pathsep + "/bin", "MAIL_RECORD": record,
                   "RECIPIENT": recipient}
            body = self.function_text() + 'send_notification_email "$RECIPIENT" "a subject" "a body"\n'
            result = subprocess.run(["bash", "-c", body], env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            called = None
            if os.path.exists(record):
                with open(record, encoding="utf-8") as handle:
                    called = handle.read().split("\n")[:-1]
            return called, result.stderr.decode("utf-8", "replace")

    def test_the_shared_email_cases(self):
        for good in CASES["email"]["good"]:
            with self.subTest(good=good):
                called, err = self.run_with_stub_mail(good)
                self.assertEqual(called, ["-s", "a subject", good], err)
        for bad in CASES["email"]["bad"]:
            with self.subTest(bad=bad):
                called, err = self.run_with_stub_mail(bad)
                self.assertIsNone(called, "the mail program was run with %r" % (bad,))
                self.assertIn("not a valid address", err)


class TestSnakefileRunsTheChecksFirst(SampleRulesCase):
    """The checks only protect anything if they run before the data is used. There is no way to run a
    Snakefile here, so the order is read from its text."""

    ROW_BLOCK = re.compile(r"_trusted = load_trusted\([^\n]*\)\nfor _row in [^\n]*\n    validate_sample\([^\n]*\)\n")
    NAME_LINE = "validate_sample_name(SAMPLE)\n"

    def snakefile_text(self):
        with open(SNAKEFILE) as handle:
            return handle.read()

    def problems_of(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "Snakefile")
            with open(path, "w") as handle:
                handle.write(text)
            return order_problems(path)

    def test_the_real_snakefile_checks_before_it_uses_anything(self):
        self.assertEqual(order_problems(SNAKEFILE), [])

    def test_the_scan_finds_what_it_looks_for(self):
        text = self.snakefile_text()
        self.assertEqual(len(self.ROW_BLOCK.findall(text)), 1)
        self.assertEqual(text.count(self.NAME_LINE), 1)

    def test_rows_checked_after_the_sample_name_is_resolved_is_found(self):
        text = self.snakefile_text()
        block = self.ROW_BLOCK.search(text).group(0)
        anchor = "SAMPLE = resolve_sample_name(samples_df, READ1_GZ, READ2_GZ)\n"
        moved = text.replace(block, "", 1).replace(anchor, anchor + block, 1)
        self.assertNotEqual(moved, text)
        self.assertIn("load_trusted comes after the first resolve_sample_name(samples_df", self.problems_of(moved))
        self.assertIn("validate_sample( comes after the first resolve_sample_name(samples_df", self.problems_of(moved))

    def test_rows_checked_after_the_rule_files_are_included_is_found(self):
        text = self.snakefile_text()
        block = self.ROW_BLOCK.search(text).group(0)
        moved = text.replace(block, "", 1).replace('include: "rules/final.smk"',
                                                   'include: "rules/final.smk"\n' + block, 1)
        found = self.problems_of(moved)
        self.assertIn("validate_sample( comes after the first include:", found)
        self.assertIn("load_trusted comes after the first include:", found)

    def test_a_read_of_samples_df_before_the_row_check_is_found(self):
        text = self.snakefile_text()
        block = self.ROW_BLOCK.search(text).group(0)
        early = text.replace(block, "_early = list(samples_df['dest'])\n" + block, 1)
        self.assertIn("samples_df is read before the row check", self.problems_of(early))

    def test_a_missing_check_is_found(self):
        text = self.snakefile_text()
        no_load = text.replace('_trusted = load_trusted("config/trusted_servers.json")\n', "_trusted = None\n", 1)
        self.assertIn('missing: load_trusted("config/trusted_servers.json")', self.problems_of(no_load))
        no_rows = text.replace("    validate_sample(_row, _trusted)\n", "    pass\n", 1)
        self.assertIn("missing: validate_sample(", self.problems_of(no_rows))
        no_name = text.replace(self.NAME_LINE, "", 1)
        self.assertIn("missing: validate_sample_name(SAMPLE)", self.problems_of(no_name))

    def test_the_derived_name_checked_after_the_includes_is_found(self):
        text = self.snakefile_text()
        moved = text.replace(self.NAME_LINE, "", 1).replace('include: "rules/final.smk"',
                                                            'include: "rules/final.smk"\n' + self.NAME_LINE, 1)
        found = self.problems_of(moved)
        self.assertIn("validate_sample_name(SAMPLE) is not right after SAMPLE = resolve_sample_name(", found)
        self.assertIn("validate_sample_name(SAMPLE) comes after the first include:", found)

    def test_the_derived_name_checked_before_it_is_derived_is_found(self):
        text = self.snakefile_text()
        anchor = "SAMPLE = resolve_sample_name(samples_df, READ1_GZ, READ2_GZ)\n"
        moved = text.replace(self.NAME_LINE, "", 1).replace(anchor, self.NAME_LINE + anchor, 1)
        self.assertIn("validate_sample_name(SAMPLE) is not right after SAMPLE = resolve_sample_name(",
                      self.problems_of(moved))

    FOLDER_CALL = 'validate_run_folder(samples_df.to_dict("records"), DATA_TABLE)\n'

    def test_the_run_folder_check_is_in_the_real_snakefile_once(self):
        self.assertEqual(self.snakefile_text().count(self.FOLDER_CALL), 1)

    def test_a_missing_run_folder_check_is_found(self):
        gone = self.snakefile_text().replace(self.FOLDER_CALL, "", 1)
        self.assertIn('missing: validate_run_folder(samples_df.to_dict("records"), DATA_TABLE)',
                      self.problems_of(gone))

    def test_a_run_folder_check_that_is_not_given_the_table_path_is_found(self):
        other = self.snakefile_text().replace(self.FOLDER_CALL, 'validate_run_folder(samples_df.to_dict("records"), None)\n', 1)
        self.assertIn('missing: validate_run_folder(samples_df.to_dict("records"), DATA_TABLE)',
                      self.problems_of(other))

    def test_the_run_folder_check_after_the_sample_name_is_resolved_is_found(self):
        text = self.snakefile_text()
        anchor = "SAMPLE = resolve_sample_name(samples_df, READ1_GZ, READ2_GZ)\n"
        moved = text.replace(self.FOLDER_CALL, "", 1).replace(anchor, anchor + self.FOLDER_CALL, 1)
        self.assertNotEqual(moved, text)
        self.assertIn("validate_run_folder( comes after the first resolve_sample_name(samples_df",
                      self.problems_of(moved))

    def test_the_run_folder_check_after_the_rule_files_are_included_is_found(self):
        text = self.snakefile_text()
        moved = text.replace(self.FOLDER_CALL, "", 1).replace('include: "rules/final.smk"',
                                                              'include: "rules/final.smk"\n' + self.FOLDER_CALL, 1)
        self.assertIn("validate_run_folder( comes after the first include:", self.problems_of(moved))

    def test_the_run_folder_check_before_the_row_check_is_found(self):
        text = self.snakefile_text()
        block = self.ROW_BLOCK.search(text).group(0)
        moved = text.replace(self.FOLDER_CALL, "", 1).replace(block, self.FOLDER_CALL + block, 1)
        self.assertIn("validate_run_folder( does not come after the validate_sample( loop",
                      self.problems_of(moved))

    def test_the_checker_ignores_comments(self):
        # a comment that names the calls must not stand in for them
        text = self.snakefile_text()
        commented = text.replace(self.NAME_LINE, "# " + self.NAME_LINE, 1)
        self.assertIn("missing: validate_sample_name(SAMPLE)", self.problems_of(commented))
        folder_commented = text.replace(self.FOLDER_CALL, "# " + self.FOLDER_CALL, 1)
        self.assertIn('missing: validate_run_folder(samples_df.to_dict("records"), DATA_TABLE)',
                      self.problems_of(folder_commented))


class TestLoadTrusted(SampleRulesCase):
    def write(self, directory, payload, name="trusted_servers.json"):
        path = os.path.join(directory, name)
        with open(path, "w") as handle:
            handle.write(payload if isinstance(payload, str) else json.dumps(payload))
        return path

    def test_absent_file_gives_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(load_trusted(os.path.join(tmp, "nope.json")))

    def test_a_valid_file_is_returned(self):
        with tempfile.TemporaryDirectory() as tmp:
            loaded = load_trusted(self.write(tmp, TRUSTED))
        self.assertEqual(loaded["servers"]["boxA"]["ssh_port"], 2222)
        self.assertEqual(loaded["igblast"]["igblast_bin"], ["igblastn"])

    def test_a_broken_file_stops_the_run_rather_than_switching_the_check_off(self):
        broken = [
            "{not json",
            "[]",
            {},
            {"servers": {}, "igblast": {}},
            dict(TRUSTED, servers={"boxA": {"results_path": "x", "ssh_port": "2222", "transfer_mode": "push"}}),
            dict(TRUSTED, servers={"boxA": {"results_path": "x", "ssh_port": 2222, "transfer_mode": "rm"}}),
            dict(TRUSTED, servers={"boxA": {"results_path": "", "ssh_port": 2222, "transfer_mode": "push"}}),
            dict(TRUSTED, igblast=dict(TRUSTED["igblast"], igblast_bin="igblastn")),
            dict(TRUSTED, igblast=dict(TRUSTED["igblast"], igblast_bin=[1])),
            dict(TRUSTED, igblast={k: v for k, v in TRUSTED["igblast"].items() if k != "igblast_data"}),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for payload in broken:
                with self.assertRaises(SampleRejected, msg=str(payload)[:60]):
                    load_trusted(self.write(tmp, payload))

    def shipped(self):
        """The trusted list of this tree, or one generated from the example config when it has none."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = trusted_list_to_test(tmp.name)
        if path is None:
            self.skipTest("this tree has neither backend/config/trusted_servers.json nor "
                          "frontend/server_config.example.json to build one from")
        return load_trusted(path)

    def test_the_shipped_file_loads_and_accepts_a_row_built_from_it(self):
        trusted = self.shipped()
        self.assertIsNotNone(trusted)
        server = sorted(trusted["servers"])[0]
        entry = trusted["servers"][server]
        igb = trusted["igblast"]
        row = good_row(server_name=server, server_results_path=entry["results_path"],
                       server_ssh_port=entry["ssh_port"], server_transfer_mode=entry["transfer_mode"],
                       igblast_bin=igb["igblast_bin"][0], igblast_data=igb["igblast_data"][0],
                       igblast_auxiliary_data=igb["igblast_auxiliary_data"][0])
        self.assertAccepted(row, trusted)
        self.assertRefused(dict(row, server_results_path="user@storage1:/path/to/storage"),
                           "server_results_path", trusted)
        # the shipped list has the old install locations that finished runs still carry
        self.assertGreaterEqual(len(igb["igblast_data"]), 2)

    def test_the_shipped_file_takes_the_rows_the_other_writers_make(self):
        trusted = self.shipped()
        igb = trusted["igblast"]
        # the lab's storage server when the list has it, else the first one (an example config has one)
        server = "storage1" if "storage1" in trusted["servers"] else sorted(trusted["servers"])[0]
        # a finished run whose igblast_db_d is the older IgH D database prefix
        legacy_d = [v for v in igb["igblast_db_d"] if v]
        self.assertEqual(len(legacy_d), 1)
        shipped = good_row(server_name=server, igblast_data=igb["igblast_data"][0],
                           igblast_auxiliary_data=igb["igblast_auxiliary_data"][0],
                           server_results_path=trusted["servers"][server]["results_path"],
                           server_ssh_port=trusted["servers"][server]["ssh_port"])
        self.assertAccepted(dict(shipped, igblast_db_d=legacy_d[0]), trusted)
        self.assertAccepted(shipped, trusted)
        self.assertRefused(dict(shipped, igblast_db_d=legacy_d[0] + "x"), "igblast_db_d", trusted)
        # benchmark/run_benchmark.sh: a trusted server name, nothing to push, no IgBLAST paths
        benchmark = good_row(server_name=server, server_path="", server_results_path="", server_ssh_port="",
                             server_transfer_mode="pull", igblast_data=None, igblast_auxiliary_data=None,
                             user_email="", sequence_number=6, run_igblast="FALSE", run_clonality="FALSE")
        self.assertAccepted(benchmark, trusted)
        # `tbrc run --local` and `tbrc run --fastq-dir`: no server name, null IgBLAST paths
        cli = dict(benchmark, server_name="", server_transfer_mode="pull")
        self.assertAccepted(cli, trusted)
        self.assertAccepted(dict(cli, server_transfer_mode="push"), trusted)
        # what is blank is only what the list says: the program and the destination still are not
        self.assertRefused(dict(benchmark, igblast_bin=None), "igblast_bin", trusted)
        self.assertRefused(dict(benchmark, server_results_path="evil@outside.example:/x/"),
                           "server_results_path", trusted)


def r_values(text, name):
    block = re.search(r"^{} <- c\((.*?)\n\)".format(re.escape(name)), text, re.S | re.M)
    assert block, name
    return re.findall(r'=\s*"([^"]*)"', block.group(1))


class TestChoicesMatchTheForm(SampleRulesCase):
    """The spec asks for two copies of the choice sets, R's and Python's. This keeps them equal."""

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(REPO, "frontend", "validate_input.R")) as handle:
            cls.validate_r = handle.read()
        with open(os.path.join(REPO, "frontend", "global.R")) as handle:
            cls.global_r = handle.read()

    def python_choices(self, field):
        return set(FIELD_RULES[field]["choices"])

    def test_choice_sets(self):
        for r_name, field in (("METHOD_CHOICES", "pipeline"),
                              ("ARCHIVE_SCOPE_CHOICES", "archive_scope"),
                              ("ARCHIVE_FORMAT_CHOICES", "archive_format"),
                              ("IGBLAST_SPECIES_CHOICES", "igblast_species"),
                              ("RECEPTOR_CHECK_CHOICES", "receptor_check"),
                              ("LUKA_LIGHT_CHOICES", "luka_light_protocol"),
                              ("MIBS_GENE_GROUPING_CHOICES", "mibs_gene_grouping"),
                              ("MIBS_CLONE_SCOPE_CHOICES", "mibs_clone_scope")):
            self.assertEqual(self.python_choices(field), set(r_values(self.validate_r, r_name)), r_name)
        self.assertEqual(self.python_choices("igblast_organism"), set(r_values(self.validate_r, "IGBLAST_SPECIES_CHOICES")))

    def test_panel_choices_are_the_form_plus_the_legacy_one(self):
        form = set(r_values(self.global_r, "igblast_panel_choices"))
        self.assertEqual(self.python_choices("igblast_panel"), form | {"ig_light"})

    def test_numeric_ranges(self):
        for name in ("mibs_min_reads_per_well", "mibs_min_wells_per_clone", "mibs_junction_max_distance",
                     "mibs_min_reads_per_cluster"):
            match = re.search(r"{}\s*=\s*c\((\d+),\s*(\d+)\)".format(name), self.validate_r)
            self.assertTrue(match, name)
            self.assertEqual(tuple(FIELD_RULES[name]["range"]), (int(match.group(1)), int(match.group(2))), name)
        n_read = re.search(r"^N_READ_RANGE <- c\((\d+),\s*(\d+)\)", self.validate_r, re.M)
        # sequence_number is twice n_read; the pipeline also takes the odd and old values
        self.assertEqual(FIELD_RULES["sequence_number"]["range"][1], 2 * int(n_read.group(2)))
        self.assertLessEqual(FIELD_RULES["sequence_number"]["range"][0], 2 * int(n_read.group(1)))

    def test_trim_sequence_pattern(self):
        match = re.search(r"TRIM_SEQUENCE_PATTERN <- \"\^\[([A-Z]+)\]\{(\d+),(\d+)\}", self.validate_r)
        self.assertTrue(match)
        self.assertEqual(sample_rules.TRIM_ALPHABET, match.group(1))
        self.assertEqual((sample_rules.TRIM_MIN, sample_rules.TRIM_MAX), (int(match.group(2)), int(match.group(3))))


try:
    import pandas  # noqa: E402
except ImportError:  # pragma: no cover
    pandas = None


@unittest.skipIf(pandas is None, "pandas is not installed here (it is in the backend image)")
class TestRowsAsTheSnakefileGetsThem(SampleRulesCase):
    """The Snakefile hands over samples_df.to_dict("records"), so NaN and None arrive as pandas makes them."""

    def records(self, *rows):
        return pandas.DataFrame(list(rows)).to_dict("records")

    def test_a_single_row_frame(self):
        for record in self.records(good_row()):
            self.assertAccepted(record)

    def test_nulls_and_missing_values_across_rows(self):
        a = good_row()
        b = good_row(run_group=None, archive_scope=None, igblast_db_v=None)
        del b["receptor_check"]          # becomes NaN in the table
        for record in self.records(a, b):
            self.assertAccepted(record)

    def test_numbers_that_became_floats(self):
        a = good_row()
        b = good_row(mibs_junction_max_distance=None)
        for record in self.records(a, b):
            self.assertAccepted(record)

    def test_a_bad_row_in_a_table_is_found(self):
        records = self.records(good_row(), good_row(user_id="x y"))
        self.assertAccepted(records[0])
        self.assertRefused(records[1], "user_id")


if __name__ == "__main__":
    unittest.main()
