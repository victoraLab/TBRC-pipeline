"""Tests for scripts/deploy_all.sh and its check of the trusted-server list.

The check (make_trusted_servers.py --check) is for an installation that has switched the list on, by
generating backend/config/trusted_servers.json. One that has not, which is every fresh clone of the
public repository, has no such file and must be able to deploy: the script says the list is off and
carries on. With the file present, drift from frontend/server_config.json still stops the deploy
before anything is copied.

The script runs in a scratch tree against a stand-in rsync that only records its arguments.

Run from the repository root:  python3 -m unittest tests/test_deploy_all.py -v
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))

import make_trusted_servers  # noqa: E402

CONFIG = {
    "storage_server": {"boxA": {"rawdata_path": "me@boxA:/vol/raw/", "results_path": "me@boxA:/vol/res/",
                                "ssh_port": 2222, "web_port": 0, "transfer_mode": "push"}},
    "cluster": {"path_prefix": "/opt/conda/condabin:/usr/bin", "snakemake_env": "smk", "igblast_bin": "igblastn",
                "igblast_data": "/pipelines/igblast/internal_data", "igblast_db_v": "", "igblast_db_d": "",
                "igblast_db_j": "", "igblast_auxiliary_data": "/pipelines/igblast/optional_file"},
}


class DeployAllTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = os.path.join(self._tmp.name, "repo")
        for folder in ("scripts", "backend/config", "frontend", "bin"):
            os.makedirs(os.path.join(self.root, folder))
        for name in ("deploy_all.sh", "make_trusted_servers.py"):
            shutil.copy(os.path.join(REPO, "scripts", name), os.path.join(self.root, "scripts", name))
        with open(os.path.join(self.root, "backend", "workflow.txt"), "w") as handle:
            handle.write("x\n")
        with open(os.path.join(self.root, "scripts", "deploy_config.sh"), "w") as handle:
            handle.write('HPC_TARGET="user@cluster:/path/to/TBRC/backend/"\n')
        self.record = os.path.join(self.root, "rsync.args")
        fake = os.path.join(self.root, "bin", "rsync")
        with open(fake, "w") as handle:
            handle.write('#!/bin/sh\nprintf "%s\\n" "$@" > "$RSYNC_RECORD"\n')
        os.chmod(fake, 0o755)

    def write_config(self, config=CONFIG):
        with open(os.path.join(self.root, "frontend", "server_config.json"), "w") as handle:
            json.dump(config, handle)

    def write_list(self, config=CONFIG):
        path = os.path.join(self.root, "backend", "config", "trusted_servers.json")
        with open(path, "w") as handle:
            handle.write(make_trusted_servers.render(make_trusted_servers.build(config)))
        return path

    def deploy(self, *args):
        env = dict(os.environ, PATH=os.path.join(self.root, "bin") + os.pathsep + os.environ["PATH"],
                   RSYNC_RECORD=self.record)
        result = subprocess.run(["bash", os.path.join(self.root, "scripts", "deploy_all.sh"), *args],
                                env=env, cwd=self.root, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        return result, os.path.exists(self.record)

    def test_without_a_list_it_says_the_list_is_off_and_deploys(self):
        # a fresh clone: no server_config.json to compare with either
        result, copied = self.deploy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(copied, "rsync was not run")
        lines = [l for l in (result.stdout + result.stderr).splitlines() if "trusted" in l.lower()]
        self.assertEqual(len(lines), 1, lines)
        self.assertIn("off", lines[0])
        self.assertIn("backend/config/trusted_servers.json", lines[0])
        self.assertNotIn("Run: python3 scripts/make_trusted_servers.py", result.stderr)

    def test_without_a_list_the_deploy_does_not_need_a_server_config(self):
        self.assertFalse(os.path.exists(os.path.join(self.root, "frontend", "server_config.json")))
        result, copied = self.deploy("--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(copied)
        with open(self.record) as handle:
            self.assertIn("--dry-run", handle.read().split("\n"))

    def test_with_a_list_in_step_with_the_config_it_deploys_quietly(self):
        self.write_config()
        self.write_list()
        result, copied = self.deploy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(copied)
        self.assertNotIn("trusted", (result.stdout + result.stderr).lower())

    def test_with_a_list_that_has_drifted_it_stops_before_copying_anything(self):
        self.write_config()
        self.write_list()
        changed = json.loads(json.dumps(CONFIG))
        changed["storage_server"]["boxA"]["results_path"] = "me@boxA:/vol/elsewhere/"
        self.write_config(changed)
        result, copied = self.deploy()
        self.assertEqual(result.returncode, 1)
        self.assertFalse(copied, "rsync ran although the list is out of step")
        self.assertIn("make_trusted_servers.py", result.stderr)

    def test_a_list_with_no_config_to_compare_with_still_stops_the_deploy(self):
        # the file is there, so the installation uses the list; an unreadable config is not "off"
        self.write_list()
        result, copied = self.deploy()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(copied)


if __name__ == "__main__":
    unittest.main()
