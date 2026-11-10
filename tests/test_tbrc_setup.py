"""Tests for `./tbrc setup`, which asks a few plain questions and writes frontend/server_config.json (and,
for the web app, deploy/docker/.env), builds the Docker images when they are missing, and starts the web app.

Everything runs in-process against a scratch checkout: the trees it writes into are temporary, every
command that would reach Docker, SSH or the pipeline is replaced by a recorder, and answers come from an
--answers file or from a scripted terminal. Nothing here needs Docker, a network or a real config.

Run from the repository root:  python3 -m unittest tests/test_tbrc_setup.py -v
"""
import contextlib
import importlib.machinery
import importlib.util
import io
import json
import os
import re
import shutil
import socket
import stat
import sys
import tempfile
import unittest
from unittest import mock

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TBRC = os.path.join(REPO, "tbrc")
REAL_EXAMPLE = os.path.join(REPO, "frontend", "server_config.example.json")

PASSWORD = "correct horse battery staple"
# a bcrypt hash of the shape caddy hash-password prints: $2a$, two digits, 53 characters
FAKE_HASH = "$2a$14$" + "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ012345678"[:53]
STAMP = "20261109-120000"
# words the shipped files must not contain, written so that this file does not contain them either
OLD_KIT_FOLDER = "mac" + "mini"
HOME_FOLDERS = "/us" + "ers/"
BACKEND_IMAGE = "tbrc-backend:arm64"
UNIFIED_IMAGE = "tbrc-unified:arm64"

# deploy/docker/env.example, as it ships
ENV_EXAMPLE = """# Copy to .env in this folder and edit. Used as: docker compose --env-file .env up -d

# Image built by ./tbrc setup (docker build --target unified).
TBRC_IMAGE=tbrc-unified:arm64

# Login for the web page. Hash from:
#   docker run --rm caddy:2-alpine caddy hash-password --plaintext 'the-password'
# Double every $ in the hash ($$2a$$14$$...), or compose will eat it.
TBRC_USER=tbrc
TBRC_PASSWORD_HASH=

# 127.0.0.1 = only this Mac can open the app (http://localhost:8080).
# 0.0.0.0   = anyone who can reach this Mac on the network (keep the password set).
TBRC_BIND=127.0.0.1
TBRC_PORT=8080

# Local time zone of the lab, for run timestamps (e.g. America/New_York, Europe/Lisbon).
TBRC_TZ=UTC

# Upload limit in MB for the app's upload mode.
TBRC_MAX_UPLOAD_MB=8192
"""

# the commands of HANDOFF.md step 6, which build the IgBLAST databases in the image
STEP6 = ("cp -r /opt/conda/share/igblast/internal_data . "
         "&& cp --update=none /opt/conda/share/igblast/optional_file/*_gl.aux optional_file/ "
         "&& ./setup_igblast_databases.sh all all "
         '&& for p in "mouse ig_all" "mouse tcr_all" "human ig_all" "human tcr_all"; '
         "do bash check_igblast_setup.sh $p | tail -1; done")

# a config of the shape of the public example, for trees that do not ship the example
FIXTURE_TEMPLATE = {
    "storage_server": {"storage1": {
        "rawdata_path": "user@storage1.example.org:/path/to/pipelines/tbrc/rawdata/",
        "results_path": "user@storage1.example.org:/path/to/pipelines/tbrc/results/",
        "ssh_port": 22, "web_port": 0, "transfer_mode": "push"}},
    "cluster": {"dtn_host": "user@dtn.hpc.example.org", "ssh_host": "login.hpc.example.org",
                "pipeline_root": "/path/to/TBRC/backend", "snakemake_env": "snakemake",
                "path_prefix": "/path/to/miniforge3/condabin:/usr/local/bin:/usr/bin",
                "igblast_bin": "igblastn", "igblast_organism": "human",
                "igblast_data": "/path/to/TBRC/backend/igblast/internal_data",
                "igblast_db_v": "", "igblast_db_d": "", "igblast_db_j": "",
                "igblast_auxiliary_data": "/path/to/TBRC/backend/igblast/optional_file",
                "notification_admin_email": "admin@example.org"},
    "metrics": {"manual_minutes_per_100mb": 29, "run_index_offset": 0},
    "barcode_defaults": {"B2": {"chain": "B", "trim_seq": "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"},
                         "T": {"chain": "T", "trim_seq": "CCAGGGTTTTCCCAGTCACGAC"}},
    "costs": {"currency": "USD"},
}

LOCAL_UPLOAD_TERMINAL = {"runs_on": "local", "data": "upload", "interface": "terminal"}
SERVER = {"data": "server", "server_name": "labnas", "server_raw": "me@nas.example.org:/vol/raw/",
          "server_results": "me@nas.example.org:/vol/results/", "server_port": 2112}
CLUSTER = {"runs_on": "cluster", "login_host": "me@login.cluster.example.org",
           "dtn_host": "me@dtn.cluster.example.org", "backend_path": "/home/me/TBRC/backend",
           "conda_env": "snakemake", "conda_path": "/home/me/miniforge3/condabin",
           "slurm_account": "lab_account", "slurm_partition": "compute"}


def load_tbrc():
    loader = importlib.machinery.SourceFileLoader("tbrc_cli_for_setup", TBRC)
    spec = importlib.util.spec_from_loader("tbrc_cli_for_setup", loader)
    module = importlib.util.module_from_spec(spec)
    keep = sys.dont_write_bytecode
    sys.dont_write_bytecode = True  # no tbrc .pyc dropped into the repository root
    try:
        loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = keep
    return module


def snapshot(root):
    """{relative path: (size, mtime_ns, bytes)} of every file under root."""
    seen = {}
    for folder, _, files in os.walk(root):
        for name in files:
            path = os.path.join(folder, name)
            if os.path.islink(path):
                continue
            info = os.stat(path)
            with open(path, "rb") as handle:
                seen[os.path.relpath(path, root)] = (info.st_size, info.st_mtime_ns, handle.read())
    return seen


class SetupCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = self._tmp.name
        self.root = os.path.join(self.tmp, "TBRC")
        for folder in ("frontend", "backend/config", "backend/igblast", "deploy/docker"):
            os.makedirs(os.path.join(self.root, folder))
        os.symlink(os.path.join(REPO, "backend", "workflow"), os.path.join(self.root, "backend", "workflow"))
        template = FIXTURE_TEMPLATE
        if os.path.exists(REAL_EXAMPLE):
            with open(REAL_EXAMPLE, encoding="utf-8") as handle:
                template = json.load(handle)
        self.write_json(os.path.join("frontend", "server_config.example.json"), template)
        self.template = template
        shutil.copy(os.path.join(REPO, "backend", "config", "cluster.json"),
                    os.path.join(self.root, "backend", "config", "cluster.json"))
        self.write_text(os.path.join("deploy", "docker", "env.example"), ENV_EXAMPLE)
        self.write_text(os.path.join("deploy", "docker", "docker-compose.yml"), "services: {}\n")
        self.write_text(os.path.join("deploy", "docker", "Dockerfile"), "FROM scratch AS backend\nFROM backend AS unified\n")
        self.write_text(os.path.join("deploy", "docker", "seed-barcodes.sh"), "#!/bin/bash\n")

        self.cli = load_tbrc()
        cli = self.cli
        cli.ROOT = self.root
        cli.FRONTEND = os.path.join(self.root, "frontend")
        cli.BACKEND = os.path.join(self.root, "backend")
        cli.CONFIG_FILE = os.path.join(cli.FRONTEND, "server_config.json")
        cli.RUNS_DIR = os.path.join(cli.FRONTEND, "runs")
        self.calls, self.checks, self.test_runs = [], [], []
        self.check_status, self.test_run_status = 0, 0
        cli.run = self.fake_run
        # Docker has both images and nothing runs natively, unless a test says otherwise
        self.images, self.native, self.docker_missing = {BACKEND_IMAGE, UNIFIED_IMAGE}, set(), False
        cli.docker_installed_images = lambda: None if self.docker_missing else set(self.images)
        cli.which = lambda tool: tool in self.native
        self.real_check = cli.cmd_check
        cli.cmd_check = self.fake_check
        cli.cmd_run = self.fake_cmd_run
        cli.timestamp = lambda: STAMP

    # ------------------------------------------------------------ doubles
    def fake_run(self, cmd, cwd=None, env=None, stream=False, timeout=None, log=None, stdin=None):
        self.calls.append({"cmd": list(cmd), "stdin": stdin, "stream": stream, "cwd": cwd})
        if cmd[:2] == ["docker", "build"]:
            self.images.add(cmd[cmd.index("-t") + 1])
        if "hash-password" in cmd:
            return 0, FAKE_HASH + "\n", ""
        return 0, "", ""

    def fake_check(self, a):
        self.checks.append(dict(vars(a)))
        return self.check_status

    def fake_cmd_run(self, a):
        self.test_runs.append(dict(vars(a)))
        return self.test_run_status

    # ------------------------------------------------------------ helpers
    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def write_text(self, relative, text, mode=None):
        path = self.path(relative)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        if mode:
            os.chmod(path, mode)

    def write_json(self, relative, data):
        self.write_text(relative, json.dumps(data, indent=2) + "\n")

    def read_text(self, *parts):
        with open(self.path(*parts), encoding="utf-8") as handle:
            return handle.read()

    def config(self):
        return json.loads(self.read_text("frontend", "server_config.json"))

    def env_lines(self):
        return dict(line.split("=", 1) for line in self.read_text("deploy", "docker", ".env").splitlines()
                    if line and not line.startswith("#") and "=" in line)

    def answers_file(self, answers):
        path = os.path.join(self.tmp, "answers.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(answers, handle)
        return path

    def setup(self, answers=None, *flags, **kwargs):
        """Run `tbrc setup`; returns (exit status, stdout, stderr). With `answers` they come from a file."""
        argv = ["setup"] + list(flags)
        if answers is not None:
            argv += ["--answers", self.answers_file(answers)]
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = self.cli.main(argv)
        return code, out.getvalue(), err.getvalue()

    def terminal(self, typed, hidden=()):
        """Patch the terminal: `typed` answers the visible questions in order, `hidden` the password ones.
        Returns the list of questions asked (the prompts)."""
        asked = []
        visible, secret = list(typed), list(hidden)

        def fake_input(prompt=""):
            asked.append(prompt)
            if not visible:
                raise EOFError
            return visible.pop(0)

        def fake_getpass(prompt="Password: ", stream=None):
            asked.append(prompt)
            if not secret:
                raise EOFError
            return secret.pop(0)

        patches = (mock.patch("builtins.input", fake_input), mock.patch.object(self.cli.getpass, "getpass", fake_getpass))
        for p in patches:
            p.start()
            self.addCleanup(p.stop)
        return asked

    def web(self, audience="me", **more):
        return dict(interface="web", audience=audience, password=PASSWORD, **more)

    def docker_builds(self):
        return [c for c in self.calls if c["cmd"][:2] == ["docker", "run"] and "bash" in c["cmd"]]

    def image_builds(self):
        return [c for c in self.calls if c["cmd"][:2] == ["docker", "build"]]

    def compose_calls(self):
        return [c for c in self.calls if c["cmd"][:2] == ["docker", "compose"]]

    def seed_calls(self):
        return [c for c in self.calls if c["cmd"][:1] == ["bash"] and "seed-barcodes.sh" in c["cmd"]]

    def hash_calls(self):
        return [c for c in self.calls if "hash-password" in c["cmd"]]


# ====================================================================== the four combinations
class TestCombinations(SetupCase):
    def test_this_computer_upload_terminal(self):
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 0, out + err)
        config = self.config()
        self.assertEqual(list(config["storage_server"]), list(self.template["storage_server"]))
        for entry in config["storage_server"].values():
            self.assertEqual(entry["transfer_mode"], "pull")
        self.assertFalse(os.path.exists(self.path("deploy", "docker", ".env")))
        self.assertEqual(self.hash_calls(), [], "no password is needed without the web app")
        # the databases are built in the image, then the local check runs
        self.assertEqual(len(self.docker_builds()), 1)
        self.assertEqual([c["local"] for c in self.checks], [True])
        self.assertEqual(self.test_runs, [])

    def test_this_computer_file_server_web_for_me(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **SERVER, **self.web("me")))
        self.assertEqual(code, 0, out + err)
        config = self.config()
        self.assertEqual(config["storage_server"], {"labnas": {
            "rawdata_path": "me@nas.example.org:/vol/raw/", "results_path": "me@nas.example.org:/vol/results/",
            "ssh_port": 2112, "web_port": 0, "transfer_mode": "pull"}})
        env = self.env_lines()
        self.assertEqual(env["TBRC_BIND"], "127.0.0.1")
        self.assertEqual(env["TBRC_USER"], "tbrc")
        self.assertEqual(env["TBRC_PORT"], "8080")      # the rest of env.example is kept
        self.assertEqual(len(self.docker_builds()), 1)

    def test_this_computer_upload_web_for_the_whole_lab(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web("lab")))
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.env_lines()["TBRC_BIND"], "0.0.0.0")
        self.assertIn("0.0.0.0", out, "the choice of the whole lab is said out loud")

    def test_cluster_file_server_terminal(self):
        code, out, err = self.setup(dict(CLUSTER, **SERVER, interface="terminal"), "--force")
        self.assertEqual(code, 0, out + err)
        config = self.config()
        cluster = config["cluster"]
        self.assertEqual(cluster["ssh_host"], "me@login.cluster.example.org")
        self.assertEqual(cluster["dtn_host"], "me@dtn.cluster.example.org")
        self.assertEqual(cluster["pipeline_root"], "/home/me/TBRC/backend")
        self.assertEqual(cluster["snakemake_env"], "snakemake")
        self.assertEqual(cluster["igblast_data"], "/home/me/TBRC/backend/igblast/internal_data")
        self.assertEqual(cluster["igblast_auxiliary_data"], "/home/me/TBRC/backend/igblast/optional_file")
        self.assertTrue(cluster["path_prefix"].startswith("/home/me/miniforge3/condabin:"))
        # a real file server on a cluster keeps the lab's usual flow: the archive is pushed to it
        self.assertEqual(config["storage_server"]["labnas"]["transfer_mode"], "push")
        slurm = json.loads(self.read_text("backend", "config", "cluster.json"))["__default__"]
        self.assertEqual((slurm["account"], slurm["partition"]), ("lab_account", "compute"))
        self.assertEqual(self.docker_builds(), [], "nothing is built on this computer for a cluster")
        self.assertEqual([c["local"] for c in self.checks], [False])
        self.assertFalse(os.path.exists(self.path("deploy", "docker", ".env")))
        # what remains to be done on the cluster is printed, with the values the person gave
        self.assertIn("/home/me/TBRC/backend", out)
        self.assertIn("rsync", out)
        self.assertIn("environment.yml", out)
        self.assertIn("setup_igblast_databases.sh", out)

    def test_cluster_upload_web_for_the_whole_lab(self):
        code, out, err = self.setup(dict(CLUSTER, data="upload", **self.web("lab")), "--force")
        self.assertEqual(code, 0, out + err)
        config = self.config()
        # no real server: the archive stays on the cluster, and the placeholder is not checked for reachability
        for entry in config["storage_server"].values():
            self.assertEqual(entry["transfer_mode"], "pull")
        self.assertTrue(self.checks[0].get("no_storage"))
        self.assertEqual(self.env_lines()["TBRC_BIND"], "0.0.0.0")
        self.assertIn("TBRC_SSH_DIR", out, "the app needs a key folder to reach the cluster")

    def test_a_cluster_answer_without_an_account_leaves_the_account_empty(self):
        code, out, err = self.setup(dict(CLUSTER, **SERVER, interface="terminal", slurm_account=""), "--force")
        self.assertEqual(code, 0, out + err)
        self.assertEqual(json.loads(self.read_text("backend", "config", "cluster.json"))["__default__"]["account"], "")

    def test_the_other_settings_of_the_example_are_kept(self):
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **SERVER))
        config = self.config()
        for key in self.template:
            if key not in ("storage_server",):
                self.assertEqual(config[key], self.template[key], key)


class TestCheckWithoutStorage(SetupCase):
    def test_uploads_only_skips_the_storage_server_checks(self):
        self.write_json(os.path.join("frontend", "server_config.json"), self.template)
        reached = []

        def ssh(host, command, port=None, timeout=60, stream=False, log=None):
            reached.append(host)
            return 0, "ROOT=ok\nDRIVER=ok\nENV=ok\nSNAKEMAKE=7.32.4\nPANDAS=ok\nACTIVE=0\nFREE_KB=99999999\n", ""
        self.cli.ssh = ssh
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.real_check(self.cli.argparse.Namespace(local=False, all=False, server="", no_storage=True))
        self.assertIn("[SKIP] storage servers", out.getvalue())
        self.assertTrue(reached, "the cluster is still checked")
        self.assertFalse([h for h in reached if "storage1" in h], reached)


class TestNeverTheGate(SetupCase):
    def test_use_gate_is_never_written(self):
        for answers in (LOCAL_UPLOAD_TERMINAL, dict(LOCAL_UPLOAD_TERMINAL, **SERVER, **self.web()),
                        dict(CLUSTER, **SERVER, interface="terminal"), dict(CLUSTER, data="upload", **self.web("lab"))):
            with self.subTest(answers=sorted(answers)):
                self.setup(answers, "--force")
                self.assertNotIn("use_gate", self.read_text("frontend", "server_config.json"))

    def test_use_gate_in_the_template_is_dropped(self):
        template = json.loads(json.dumps(self.template))
        template.setdefault("cluster", {})["use_gate"] = True
        self.write_json(os.path.join("frontend", "server_config.example.json"), template)
        self.setup(dict(CLUSTER, **SERVER, interface="terminal"), "--force")
        self.assertNotIn("use_gate", self.read_text("frontend", "server_config.json"))

    def test_a_local_run_is_never_pushed(self):
        # with "push" a local run tries to copy its archive to the file server and is recorded as failed
        for extra in ({}, SERVER):
            with self.subTest(server=bool(extra)):
                self.setup(dict(LOCAL_UPLOAD_TERMINAL, **extra), "--force")
                self.assertEqual({e["transfer_mode"] for e in self.config()["storage_server"].values()}, {"pull"})


class TestTheEnvFile(SetupCase):
    def test_mode_600_and_every_dollar_doubled(self):
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        path = self.path("deploy", "docker", ".env")
        self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), 0o600)
        value = self.env_lines()["TBRC_PASSWORD_HASH"]
        self.assertEqual(value, FAKE_HASH.replace("$", "$$"))
        self.assertIsNone(re.search(r"(?<!\$)\$(?!\$)", value), "a single $ would be eaten by compose")
        self.assertEqual(value.count("$$"), FAKE_HASH.count("$"))

    def test_the_password_is_hashed_on_stdin_and_never_on_a_command_line(self):
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        hashes = self.hash_calls()
        self.assertEqual(len(hashes), 1)
        self.assertEqual(hashes[0]["cmd"], ["docker", "run", "--rm", "-i", "caddy:2-alpine", "caddy", "hash-password"])
        self.assertEqual(hashes[0]["stdin"], PASSWORD + "\n")
        for call in self.calls:
            self.assertFalse(any(PASSWORD in word for word in call["cmd"]), call["cmd"])
            self.assertFalse(any("--plaintext" == word for word in call["cmd"]))

    def test_the_password_and_its_hash_are_never_printed_or_stored_elsewhere(self):
        for flags in ((), ("--dry-run",)):
            with self.subTest(flags=flags):
                code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **SERVER, **self.web("lab")), "--force", *flags)
                self.assertEqual(code, 0, out + err)
                for secret in (PASSWORD, FAKE_HASH, FAKE_HASH.replace("$", "$$"), FAKE_HASH[7:]):
                    self.assertNotIn(secret, out + err)
        for relative in ("frontend/server_config.json", "backend/config/cluster.json"):
            self.assertNotIn(PASSWORD, self.read_text(*relative.split("/")))

    def test_a_failed_hash_is_reported_without_the_password_and_writes_nothing(self):
        def failing(cmd, **kwargs):
            if "hash-password" in cmd:
                return 1, "", "Error: could not hash %s" % kwargs.get("stdin", "").strip()
            return 0, "", ""
        self.cli.run = failing
        before = snapshot(self.root)
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 2)
        self.assertNotIn(PASSWORD, out + err)
        self.assertEqual(snapshot(self.root), before)

    def test_a_world_readable_answers_file_with_a_password_is_warned_about(self):
        answers = dict(LOCAL_UPLOAD_TERMINAL, **self.web())
        path = self.answers_file(answers)
        for mode, warned in ((0o644, True), (0o600, False)):
            os.chmod(path, mode)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.cli.main(["setup", "--answers", path, "--force"])
            self.assertEqual("holds a password in the clear" in out.getvalue(), warned, oct(mode))
            self.assertNotIn(PASSWORD, out.getvalue())

    def test_a_password_hash_can_be_given_instead_of_a_password(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me", password_hash=FAKE_HASH))
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.hash_calls(), [])
        self.assertEqual(self.env_lines()["TBRC_PASSWORD_HASH"], FAKE_HASH.replace("$", "$$"))

    def test_a_password_is_always_required(self):
        before = snapshot(self.root)
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me"))
        self.assertEqual(code, 2)
        self.assertIn("password", out + err)
        self.assertEqual(snapshot(self.root), before)
        for bad in ("", "short", "bad\nnewline", "nul\x00byte"):
            with self.subTest(password=bad):
                code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me", password=bad))
                self.assertEqual(code, 2)
                self.assertEqual(snapshot(self.root), before)

    def test_a_hash_that_is_not_bcrypt_is_refused(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me", password_hash="$1$x$y"))
        self.assertEqual(code, 2)
        self.assertIn("password_hash", out + err)

    def test_without_deploy_docker_it_says_the_folder_is_part_of_the_repository_and_goes_on(self):
        shutil.rmtree(self.path("deploy"))
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me"))
        self.assertEqual(code, 0, out + err)
        self.assertIn("deploy/docker", out)
        self.assertIn("git pull", out)
        self.assertTrue(os.path.exists(self.path("frontend", "server_config.json")))
        self.assertFalse(os.path.exists(self.path("deploy", "docker", ".env")))
        self.assertEqual(self.hash_calls(), [])
        self.assertEqual(self.image_builds(), [], "no Dockerfile to build the web app image from")

    def test_with_deploy_docker_the_env_file_is_written_there_with_mode_600(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 0, out + err)
        path = self.path("deploy", "docker", ".env")
        self.assertTrue(os.path.isfile(path))
        self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), 0o600)
        self.assertEqual(self.env_lines()["TBRC_IMAGE"], UNIFIED_IMAGE)

    def test_the_next_step_for_the_web_app_is_printed(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertIn("docker compose --env-file .env up -d", out)
        self.assertIn("seed-barcodes.sh", out)


# ====================================================================== refusing bad answers
class TestRefusals(SetupCase):
    def refused(self, answers, key):
        before = snapshot(self.root)
        code, out, err = self.setup(answers)
        self.assertEqual(code, 2, out + err)
        self.assertIn(key, out + err)
        self.assertEqual(snapshot(self.root), before, "nothing may be written when an answer is refused")
        self.assertEqual(self.calls, [], "nothing may be run when an answer is refused")
        return out + err

    def test_a_bad_host_is_refused(self):
        for bad in ("-oProxyCommand=touch@x", "-o", "host with space", "a;b", "$(x)", "me@", "@host", "me@@host",
                    "host..name", "", "host\nname", "caf\u00e9", "x" * 400):
            with self.subTest(login_host=bad):
                self.refused(dict(CLUSTER, **SERVER, interface="terminal", login_host=bad), "login_host")
        for bad in ("-x", "a b", "me@-x"):
            with self.subTest(dtn_host=bad):
                self.refused(dict(CLUSTER, **SERVER, interface="terminal", dtn_host=bad), "dtn_host")

    def test_a_bad_file_server_path_is_refused(self):
        for key in ("server_raw", "server_results"):
            for bad in ("me@nas:/my folder/", "me@nas:/a/../b/", "nas:/a/", "me@nas:relative/", "me@nas:/~/x",
                        "me@nas:/-rf/", "me@nas:/a;b/", "me@nas", "/just/a/path/", "", "me@nas:/a$(x)/", "me@nas:/a\nb/"):
                with self.subTest(**{key: bad}):
                    self.refused(dict(LOCAL_UPLOAD_TERMINAL, **dict(SERVER, **{key: bad})), key)

    def test_a_bad_server_name_is_refused(self):
        for bad in ("", "local", "upload", "LOCAL", "a b", "a/b", "..", ".", "-x", "a;b", "x" * 80):
            with self.subTest(server_name=bad):
                self.refused(dict(LOCAL_UPLOAD_TERMINAL, **dict(SERVER, server_name=bad)), "server_name")

    def test_a_bad_port_is_refused(self):
        for bad in ("0", "65536", "-1", "22; ls", "http", "2.5", 22.5, True, [22]):
            with self.subTest(server_port=bad):
                self.refused(dict(LOCAL_UPLOAD_TERMINAL, **dict(SERVER, server_port=bad)), "server_port")

    def test_bad_cluster_answers_are_refused(self):
        cases = (("backend_path", ("", "relative/path", "/a b/c", "/a/../b", "/-x", "~/x", "/a;b", "/")),
                 ("conda_env", ("a b", "-x", "a;b", "..")),
                 ("conda_path", ("relative", "/a b", "/a/../b", "/a:/b", "/a;b")),
                 ("slurm_account", ("a b", "-x", "a;b", "$(x)", "caf\u00e9")),
                 ("slurm_partition", ("", "a b", "-x", "a;b", "x|y")))
        for key, values in cases:
            for bad in values:
                with self.subTest(**{key: bad}):
                    self.refused(dict(CLUSTER, **SERVER, interface="terminal", **{key: bad}), key)

    def test_a_choice_that_is_not_one_of_the_choices_is_refused(self):
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, runs_on="the moon"), "runs_on")
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, data="carrier pigeon"), "data")
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, interface="telepathy"), "interface")
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, **self.web("the moon")), "audience")

    def test_a_missing_or_unknown_answer_is_named(self):
        self.refused({"data": "upload", "interface": "terminal"}, "runs_on")
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, **dict(SERVER, server_name=None)), "server_name")
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, flavour="vanilla"), "flavour")

    def test_a_test_folder_that_is_not_a_folder_is_refused(self):
        self.refused(dict(LOCAL_UPLOAD_TERMINAL, test_data=os.path.join(self.tmp, "nope")), "test_data")

    def test_a_broken_answers_file_is_refused(self):
        path = os.path.join(self.tmp, "broken.json")
        for text in ("{not json", "[]", '"x"'):
            with open(path, "w") as handle:
                handle.write(text)
            code, out, err = self.setup(None, "--answers", path)
            self.assertEqual(code, 2, text)
            self.assertIn("answers", out + err)
        code, out, err = self.setup(None, "--answers", os.path.join(self.tmp, "absent.json"))
        self.assertEqual(code, 2)

    def test_the_validators_agree_with_the_workflow_rules(self):
        """Names and paths are judged by backend/workflow/sample_rules.py, so the shared table of good and
        bad path segments (tests/validation_cases.json, also read by the form's tests) holds here too."""
        with open(os.path.join(REPO, "tests", "validation_cases.json"), encoding="utf-8") as handle:
            cases = json.load(handle)["path_segment"]
        for good in cases["good"]:
            self.assertEqual(self.cli.check_segment(good), good)
        for bad in cases["bad"] + cases["bad_untrimmed"]:
            with self.assertRaises(ValueError, msg=repr(bad)):
                self.cli.check_segment(bad)


# ====================================================================== existing files
class TestExistingFiles(SetupCase):
    def test_an_overwritten_file_is_kept_as_a_dated_backup(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--force")
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.read_text("frontend", "server_config.json.bak-" + STAMP), '{"old": true}\n')
        self.assertIn("storage_server", self.config())
        self.assertIn("server_config.json.bak-" + STAMP, out)

    def test_a_second_backup_in_the_same_second_does_not_replace_the_first(self):
        self.write_text("frontend/server_config.json", '{"old": 1}\n')
        self.setup(LOCAL_UPLOAD_TERMINAL, "--force")
        self.setup(LOCAL_UPLOAD_TERMINAL, "--force")
        self.assertEqual(self.read_text("frontend", "server_config.json.bak-" + STAMP), '{"old": 1}\n')
        backups = sorted(n for n in os.listdir(self.path("frontend")) if ".bak-" in n)
        self.assertEqual(len(backups), 2, backups)

    def test_a_backup_of_the_env_file_is_private_too(self):
        self.write_text("deploy/docker/.env", "TBRC_USER=old\n", mode=0o600)
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--force")
        self.assertEqual(code, 0, out + err)
        backup = self.path("deploy", "docker", ".env.bak-" + STAMP)
        self.assertEqual(stat.S_IMODE(os.stat(backup).st_mode), 0o600)
        self.assertEqual(stat.S_IMODE(os.stat(self.path("deploy", "docker", ".env")).st_mode), 0o600)
        self.assertNotEqual(self.env_lines()["TBRC_USER"], "old")

    def test_a_cluster_setup_needs_force_without_a_terminal_because_it_changes_cluster_json(self):
        before = snapshot(self.root)
        code, out, err = self.setup(dict(CLUSTER, **SERVER, interface="terminal"))
        self.assertEqual(code, 2)
        self.assertIn("cluster.json", out + err)
        self.assertEqual(snapshot(self.root), before)

    def test_the_cluster_settings_file_is_backed_up_too(self):
        original = self.read_text("backend", "config", "cluster.json")
        self.setup(dict(CLUSTER, **SERVER, interface="terminal"), "--force")
        self.assertEqual(self.read_text("backend", "config", "cluster.json.bak-" + STAMP), original)

    def test_without_force_a_non_interactive_run_stops_before_writing_anything(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        before = snapshot(self.root)
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 2)
        self.assertIn("--force", out + err)
        self.assertIn("server_config.json", out + err)
        self.assertEqual(snapshot(self.root), before, "not even the .env may be written when the config is refused")
        self.assertEqual(self.hash_calls(), [])

    def test_in_a_terminal_it_asks_and_a_no_keeps_the_file(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        asked = self.terminal(["1", "1", "1", "", "n"])     # runs here, upload, terminal; no test data; keep the file
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.read_text("frontend", "server_config.json"), '{"old": true}\n')
        self.assertEqual([n for n in os.listdir(self.path("frontend")) if ".bak-" in n], [])
        self.assertTrue(any("server_config.json" in q for q in asked))

    def test_in_a_terminal_it_asks_and_a_yes_replaces_it_with_a_backup(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        self.terminal(["1", "1", "1", "", "y"])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.read_text("frontend", "server_config.json.bak-" + STAMP), '{"old": true}\n')
        self.assertIn("storage_server", self.config())


# ====================================================================== --dry-run
class TestDryRun(SetupCase):
    def test_a_dry_run_writes_and_runs_nothing(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        before = snapshot(self.root)
        for answers in (LOCAL_UPLOAD_TERMINAL, dict(LOCAL_UPLOAD_TERMINAL, **SERVER, **self.web("lab")),
                        dict(CLUSTER, **SERVER, interface="terminal")):
            code, out, err = self.setup(dict(answers, test_data=self.tmp), "--dry-run", "--force")
            self.assertEqual(code, 0, out + err)
            self.assertIn("would write", out)
            self.assertEqual(snapshot(self.root), before)
            self.assertEqual(sorted(os.listdir(self.path("frontend"))),
                             ["server_config.example.json", "server_config.json"], "no backup, no new file")
        self.assertEqual((self.calls, self.checks, self.test_runs), ([], [], []))

    def test_a_dry_run_says_what_it_would_run(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web("lab"), test_data=self.tmp), "--dry-run")
        self.assertIn("would run", out)
        self.assertIn("setup_igblast_databases.sh", out)
        self.assertIn("hash-password", out)
        self.assertIn("TBRC_BIND=0.0.0.0", out)
        self.assertIn("TBRC_PASSWORD_HASH", out)
        self.assertIn("TCRb_mini", out)

    def test_a_dry_run_needs_no_password(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me"), "--dry-run")
        self.assertEqual(code, 0, out + err)

    def test_a_dry_run_stops_like_the_real_run_when_a_file_exists_and_force_is_missing(self):
        self.write_text("frontend/server_config.json", '{"old": true}\n')
        before = snapshot(self.root)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--dry-run")
        self.assertEqual(code, 2)
        self.assertIn("--force", out + err)
        self.assertEqual(snapshot(self.root), before)

    def test_a_dry_run_still_refuses_bad_answers(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **dict(SERVER, server_name="local")), "--dry-run")
        self.assertEqual(code, 2)


# ====================================================================== IgBLAST databases and the check
class TestLocalBuild(SetupCase):
    def test_the_build_is_the_one_in_the_handoff_step(self):
        self.setup(LOCAL_UPLOAD_TERMINAL)
        build = self.docker_builds()[0]
        self.assertEqual(build["cmd"], ["docker", "run", "--rm", "-v", self.path("backend") + ":/work", "-w",
                                        "/work/igblast", "tbrc-backend:arm64", "bash", "-c", STEP6])

    def mark_built(self):
        os.makedirs(self.path("backend", "igblast", "internal_data"))
        for species in ("mouse", "human"):
            for panel in ("ig_all", "tcr_all"):
                folder = self.path("backend", "igblast", "work", species, panel)
                os.makedirs(folder)
                with open(os.path.join(folder, "%s_%s_V.nhr" % (species, panel)), "w") as handle:
                    handle.write("x")

    def test_databases_that_are_already_built_are_not_built_again(self):
        self.mark_built()
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.docker_builds(), [])
        self.assertIn("already built", out)
        self.assertEqual(len(self.checks), 1)

    def test_one_missing_panel_means_it_builds(self):
        self.mark_built()
        os.remove(self.path("backend", "igblast", "work", "human", "tcr_all", "human_tcr_all_V.nhr"))
        self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(len(self.docker_builds()), 1)

    def test_a_failed_build_is_reported_and_setup_does_not_claim_success(self):
        def run(cmd, **kwargs):
            self.calls.append({"cmd": list(cmd), "stdin": kwargs.get("stdin"), "stream": kwargs.get("stream")})
            return (1, "", "") if "bash" in cmd else (0, "", "")
        self.cli.run = run
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 1)
        self.assertIn("[FAIL]", out)

    def test_a_failing_check_gives_exit_status_1_but_the_files_are_written(self):
        self.check_status = 1
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 1)
        self.assertTrue(os.path.exists(self.path("frontend", "server_config.json")))


# ====================================================================== --test-data
class TestTestData(SetupCase):
    def library(self):
        folder = os.path.join(self.tmp, "mini")
        os.makedirs(folder)
        return folder

    def test_the_mini_library_is_run_with_the_settings_in_the_brief(self):
        folder = self.library()
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--test-data", folder)
        self.assertEqual(code, 0, out + err)
        run = self.test_runs[0]
        self.assertEqual((run["user"], run["run_name"], run["primer_set"], run["pipeline"], run["species"], run["panel"],
                          run["clonality"], run["fastq_dir"]),
                         ("setup_test", "TCRb_mini", "T", "tiago_vsearch", "mouse", "trb", True, folder))
        self.assertTrue(run["local"])
        self.assertIn("TCRb_mini", out)
        self.assertIn("finished", out)

    def test_on_a_cluster_it_is_not_run_locally(self):
        folder = self.library()
        self.setup(dict(CLUSTER, **SERVER, interface="terminal", test_data=folder), "--force")
        self.assertFalse(self.test_runs[0]["local"])

    def test_a_run_that_did_not_finish_is_reported(self):
        self.test_run_status = 2
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--test-data", self.library())
        self.assertEqual(code, 1)
        self.assertIn("[FAIL]", out)
        self.assertIn("TCRb_mini", out)

    def test_a_run_that_raises_a_refusal_is_reported_with_its_message(self):
        def refusing(a):
            raise self.cli.Fail("sample.json refused, 1 problem", "a hint")
        self.cli.cmd_run = refusing
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--test-data", self.library())
        self.assertEqual(code, 1)
        self.assertIn("sample.json refused", out)

    def test_blank_skips_it(self):
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, test_data=""))
        self.setup(LOCAL_UPLOAD_TERMINAL, "--force")
        self.assertEqual(self.test_runs, [])

    def test_the_folder_may_be_asked_for_in_a_terminal(self):
        folder = self.library()
        asked = self.terminal(["1", "1", "1", folder])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.test_runs[0]["fastq_dir"], folder)
        self.assertTrue(any("test" in q.lower() for q in asked))


# ====================================================================== the questions in a terminal
class TestTerminal(SetupCase):
    def test_every_question_is_asked_one_at_a_time_in_plain_words(self):
        asked = self.terminal(["2", "2", "labnas", "me@nas.example.org:/raw/", "me@nas.example.org:/res/", "",
                               "1", "me@login.example.org", "", "/home/me/TBRC/backend", "", "", "a1", "compute",
                               "",       # no test library
                               "y"],     # yes, change backend/config/cluster.json
                              hidden=[])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        config = self.config()
        self.assertEqual(config["storage_server"]["labnas"]["ssh_port"], 22, "the port defaults to 22")
        self.assertEqual(config["storage_server"]["labnas"]["results_path"], "me@nas.example.org:/res/")
        self.assertEqual(config["cluster"]["ssh_host"], "me@login.example.org")
        self.assertEqual(config["cluster"]["snakemake_env"], "snakemake", "the conda env defaults to snakemake")
        text = " ".join(asked).lower()
        for word in ("computer", "cluster", "ssh", "upload", "terminal", "web app"):
            self.assertIn(word, out.lower() + text)

    def test_a_bad_answer_in_a_terminal_is_asked_again(self):
        asked = self.terminal(["maybe", "3", "1",                 # where: two bad tries, then this computer
                               "2", "local", "labnas", "me@nas:/raw/", "me@nas:/my res/", "me@nas:/res/", "70000", "22",
                               "1", ""])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        config = self.config()
        self.assertEqual(list(config["storage_server"]), ["labnas"], "the refused name 'local' is gone")
        self.assertEqual(config["storage_server"]["labnas"]["results_path"], "me@nas:/res/")
        self.assertEqual(config["storage_server"]["labnas"]["ssh_port"], 22)
        self.assertGreaterEqual(out.count("Not accepted"), 5)

    def test_the_password_is_asked_twice_and_must_match(self):
        self.terminal(["1", "1", "2", "1", "", "n"],      # ... no test library; do not start the web app now
                      hidden=["a-long-password", "something else", PASSWORD, PASSWORD])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertIn("differ", out)
        self.assertEqual(self.hash_calls()[0]["stdin"], PASSWORD + "\n")
        self.assertNotIn(PASSWORD, out + err)

    def test_input_that_ends_early_is_a_clean_refusal(self):
        self.terminal(["1"])
        before = snapshot(self.root)
        code, out, err = self.setup(None)
        self.assertEqual(code, 2)
        self.assertEqual(snapshot(self.root), before)


# ====================================================================== the Docker images
class TestImages(SetupCase):
    BUILD_BACKEND = ["docker", "build", "--target", "backend", "-t", BACKEND_IMAGE]
    BUILD_UNIFIED = ["docker", "build", "--target", "unified", "-t", UNIFIED_IMAGE]

    def context(self):
        return self.path("deploy", "docker")

    def test_a_missing_image_without_a_terminal_stops_with_the_exact_command_and_changes_nothing(self):
        self.images.clear()
        before = snapshot(self.root)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 2, out + err)
        text = out + err
        self.assertIn(BACKEND_IMAGE, text)
        self.assertIn(" ".join(self.BUILD_BACKEND + [self.context()]), text)
        self.assertIn("--build-images", text)
        self.assertIn("15 minutes", text)
        self.assertEqual(snapshot(self.root), before)
        self.assertEqual(self.calls, [], "nothing is built or run without --build-images")

    def test_the_web_app_needs_the_unified_image_and_both_are_named_when_both_are_missing(self):
        self.images.clear()
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 2)
        text = out + err
        for image in (BACKEND_IMAGE, UNIFIED_IMAGE):
            self.assertIn(image, text)
        self.assertEqual(self.calls, [])

    def test_a_cluster_with_the_web_app_needs_only_the_unified_image(self):
        self.images.discard(BACKEND_IMAGE)
        code, out, err = self.setup(dict(CLUSTER, data="upload", **self.web()), "--force")
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.image_builds(), [])

    def test_build_images_builds_what_is_missing_with_the_dockerfiles_real_targets(self):
        self.images.clear()
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--build-images")
        self.assertEqual(code, 0, out + err)
        builds = self.image_builds()
        self.assertEqual([b["cmd"] for b in builds],
                         [self.BUILD_BACKEND + [self.context()], self.BUILD_UNIFIED + [self.context()]])
        self.assertTrue(all(b["stream"] for b in builds))
        # the images exist before the IgBLAST databases are built in one of them
        order = [("build" if c["cmd"][:2] == ["docker", "build"] else "igblast" if "bash" in c["cmd"] else "other")
                 for c in self.calls]
        self.assertLess(order.index("build"), order.index("igblast"))
        self.assertIn("15 minutes", out)

    def test_build_images_builds_only_the_missing_one(self):
        self.images.discard(UNIFIED_IMAGE)
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--build-images")
        self.assertEqual([b["cmd"] for b in self.image_builds()], [self.BUILD_UNIFIED + [self.context()]])

    def test_build_images_with_nothing_missing_builds_nothing(self):
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--build-images")
        self.assertEqual(self.image_builds(), [])

    def test_a_failed_build_is_reported_and_the_exit_status_is_1(self):
        self.images.clear()

        def run(cmd, **kwargs):
            self.calls.append({"cmd": list(cmd), "stdin": kwargs.get("stdin"), "stream": kwargs.get("stream"),
                               "cwd": kwargs.get("cwd")})
            return (1, "", "") if cmd[:2] == ["docker", "build"] else (0, "", "")
        self.cli.run = run
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--build-images")
        self.assertEqual(code, 1, out + err)
        self.assertIn("[FAIL]", out)
        self.assertTrue(os.path.exists(self.path("frontend", "server_config.json")), "the files were written first")

    def test_in_a_terminal_it_says_what_is_missing_and_asks_and_yes_builds_it(self):
        self.images.clear()
        asked = self.terminal(["1", "1", "1", "", "y"])      # this computer, upload, terminal, no test library, yes
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual([b["cmd"] for b in self.image_builds()], [self.BUILD_BACKEND + [self.context()]])
        self.assertTrue(any("build" in q.lower() and "now" in q.lower() for q in asked), asked)
        self.assertIn("15 minutes", out)

    def test_in_a_terminal_a_no_prints_the_commands_and_goes_on_without_the_image(self):
        self.images.clear()
        self.terminal(["1", "1", "1", "", "n"])
        code, out, err = self.setup(None)
        self.assertEqual(self.image_builds(), [])
        self.assertEqual(code, 1, "setup is not finished without the image")
        self.assertIn(" ".join(self.BUILD_BACKEND + [self.context()]), out)
        self.assertTrue(os.path.exists(self.path("frontend", "server_config.json")), "the configuration is kept")

    def test_a_dry_run_says_what_it_would_build_and_stops_like_the_real_run_without_the_flag(self):
        self.images.clear()
        before = snapshot(self.root)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--dry-run")
        self.assertEqual(code, 2)
        self.assertIn("would build", out)
        self.assertEqual(snapshot(self.root), before)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL, "--dry-run", "--build-images")
        self.assertEqual(code, 0, out + err)
        self.assertIn(" ".join(self.BUILD_BACKEND + [self.context()]), out)
        self.assertEqual(self.calls, [])

    def test_tools_on_the_path_mean_no_backend_image_is_needed(self):
        self.images.clear()
        self.native = set(self.cli.LOCAL_TOOLS)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.image_builds(), [])

    def test_without_docker_it_says_so(self):
        self.docker_missing = True
        before = snapshot(self.root)
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 2)
        self.assertIn("Docker", out + err)
        self.assertEqual(snapshot(self.root), before)

    def test_the_image_names_are_the_ones_tbrc_and_the_compose_file_look_for(self):
        self.assertIn(self.cli.SETUP_BACKEND_IMAGE, self.cli.DOCKER_IMAGES)
        self.assertEqual((self.cli.SETUP_BACKEND_IMAGE, self.cli.SETUP_UNIFIED_IMAGE), (BACKEND_IMAGE, UNIFIED_IMAGE))
        with open(os.path.join(REPO, "deploy", "docker", "env.example"), encoding="utf-8") as handle:
            self.assertIn("TBRC_IMAGE=" + UNIFIED_IMAGE, handle.read().splitlines())
        with open(os.path.join(REPO, "deploy", "docker", "docker-compose.yml"), encoding="utf-8") as handle:
            self.assertIn("${TBRC_IMAGE:-%s}" % UNIFIED_IMAGE, handle.read())


# ====================================================================== starting the web app
class TestStart(SetupCase):
    COMPOSE = ["docker", "compose", "--env-file", ".env", "up", "-d"]

    def test_start_runs_compose_with_the_env_file_then_seeds_the_barcode_sets(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--start")
        self.assertEqual(code, 0, out + err)
        compose, seed = self.compose_calls(), self.seed_calls()
        self.assertEqual([c["cmd"] for c in compose], [self.COMPOSE])
        self.assertEqual([c["cmd"] for c in seed], [["bash", "seed-barcodes.sh"]])
        folder = self.path("deploy", "docker")
        self.assertEqual((compose[0]["cwd"], seed[0]["cwd"]), (folder, folder))
        order = [c["cmd"][:2] for c in self.calls]
        self.assertLess(order.index(["docker", "compose"]), order.index(["bash", "seed-barcodes.sh"]))
        self.assertIn("http://localhost:8080", out)

    def test_without_start_a_run_without_a_terminal_only_prints_the_commands(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 0, out + err)
        self.assertEqual((self.compose_calls(), self.seed_calls()), ([], []))
        self.assertIn("docker compose --env-file .env up -d", out)
        self.assertIn("bash seed-barcodes.sh", out)
        self.assertIn("cd deploy/docker", out)

    def test_the_terminal_only_choice_starts_nothing(self):
        self.setup(LOCAL_UPLOAD_TERMINAL, "--start")
        self.assertEqual((self.compose_calls(), self.seed_calls()), ([], []))

    def test_in_a_terminal_it_asks_and_yes_starts_it(self):
        asked = self.terminal(["1", "1", "2", "1", "", "y"], hidden=[PASSWORD, PASSWORD])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual([c["cmd"] for c in self.compose_calls()], [self.COMPOSE])
        self.assertEqual(len(self.seed_calls()), 1)
        self.assertTrue(any("start" in q.lower() for q in asked), asked)

    def test_in_a_terminal_a_no_prints_the_commands(self):
        self.terminal(["1", "1", "2", "1", "", "n"], hidden=[PASSWORD, PASSWORD])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertEqual((self.compose_calls(), self.seed_calls()), ([], []))
        self.assertIn("docker compose --env-file .env up -d", out)

    def test_it_does_not_start_without_the_unified_image(self):
        self.images.discard(UNIFIED_IMAGE)
        self.terminal(["1", "1", "2", "1", "", "n", "y"], hidden=[PASSWORD, PASSWORD])
        code, out, err = self.setup(None)
        self.assertEqual(self.compose_calls(), [])
        self.assertEqual(code, 1)
        self.assertIn(" ".join(TestImages.BUILD_UNIFIED), out)

    def test_a_failed_start_is_reported(self):
        def run(cmd, **kwargs):
            self.calls.append({"cmd": list(cmd), "stdin": kwargs.get("stdin"), "stream": kwargs.get("stream"),
                               "cwd": kwargs.get("cwd")})
            if cmd[:2] == ["docker", "compose"]:
                return 1, "", ""
            return (0, FAKE_HASH + "\n", "") if "hash-password" in cmd else (0, "", "")
        self.cli.run = run
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--start")
        self.assertEqual(code, 1)
        self.assertIn("[FAIL]", out)
        self.assertEqual(self.seed_calls(), [], "no seeding when the app did not start")

    def test_a_dry_run_starts_nothing_and_says_so(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--dry-run", "--start")
        self.assertEqual(code, 0, out + err)
        self.assertIn("would run: docker compose --env-file .env up -d", out)
        self.assertEqual((self.compose_calls(), self.seed_calls()), ([], []))


# ====================================================================== small things the independent review found
class TestReviewMinors(SetupCase):
    def test_a_password_over_72_bytes_is_refused_with_the_reason_and_asked_again(self):
        # bcrypt takes 72 bytes; caddy hash-password exits 1 beyond that, which used to end setup after every question
        too_long = "p" * 73
        self.terminal(["1", "1", "2", "1", "", "n"], hidden=[too_long, too_long, PASSWORD, PASSWORD])
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        self.assertIn("72 bytes", out)
        self.assertNotIn(too_long, out + err)
        self.assertEqual(self.hash_calls()[0]["stdin"], PASSWORD + "\n")

    def test_72_bytes_pass_and_73_do_not_counting_bytes_not_characters(self):
        check = self.cli.check_password
        self.assertEqual(check("p" * 72), "p" * 72)
        self.assertEqual(check("\u00e9" * 36), "\u00e9" * 36)          # 36 characters, 72 bytes
        for bad in ("p" * 73, "\u00e9" * 37, "\U0001f600" * 19):       # 73 bytes, 74 bytes, 76 bytes
            with self.assertRaises(ValueError) as caught:
                check(bad)
            self.assertIn("72 bytes", str(caught.exception))

    def test_an_answers_file_with_a_long_password_is_refused_by_name(self):
        before = snapshot(self.root)
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, interface="web", audience="me", password="p" * 73))
        self.assertEqual(code, 2)
        self.assertIn("password", out + err)
        self.assertIn("72 bytes", out + err)
        self.assertEqual(snapshot(self.root), before)

    def test_the_env_backup_never_exists_with_a_wider_mode_than_600(self):
        self.write_text("deploy/docker/.env", "TBRC_USER=old\n", mode=0o644)
        opened = []
        real_open = os.open

        def recording_open(path, flags, mode=0o777, **kwargs):
            opened.append((os.fspath(path), mode))
            return real_open(path, flags, mode, **kwargs)
        old_umask = os.umask(0o000)       # the worst case: a umask that removes nothing
        self.addCleanup(os.umask, old_umask)
        with mock.patch("os.open", recording_open), \
                mock.patch.object(shutil, "copy2", side_effect=AssertionError("copy2 makes the copy at the umask default")):
            code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--force")
        self.assertEqual(code, 0, out + err)
        backups = [(p, m) for p, m in opened if ".env.bak-" in p]
        self.assertEqual(len(backups), 1, opened)
        self.assertEqual(backups[0][1] & 0o077, 0, "the backup is created with %o" % backups[0][1])
        backup = self.path("deploy", "docker", ".env.bak-" + STAMP)
        self.assertEqual(stat.S_IMODE(os.stat(backup).st_mode), 0o600)
        with open(backup) as handle:
            self.assertEqual(handle.read(), "TBRC_USER=old\n")

    def test_a_write_error_after_earlier_writes_says_what_was_already_written(self):
        real = self.cli.write_file
        calls = []

        def failing(path, text, private, stamp):
            calls.append(path)
            if len(calls) == 2:
                raise OSError(28, "No space left on device")
            return real(path, text, private, stamp)
        self.cli.write_file = failing
        code, out, err = self.setup(dict(CLUSTER, **SERVER, interface="terminal"), "--force")
        self.assertEqual(code, 2)
        self.assertIn("wrote frontend/server_config.json", out)
        self.assertIn("No space left", out)
        flat = " ".join(out.split())
        self.assertIn("already written", flat)
        self.assertIn("server_config.json", flat.split("already written", 1)[1])

    def test_a_ctrl_c_after_writing_lists_what_was_written(self):
        def interrupted(a):
            raise KeyboardInterrupt
        self.cli.cmd_check = interrupted
        code, out, err = self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(code, 130)
        flat = " ".join(out.split())
        self.assertIn("interrupted", flat)
        self.assertIn("already written", flat)
        self.assertIn("frontend/server_config.json", flat)

    def test_the_exit_codes_are_described_truthfully_in_the_readme(self):
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = " ".join(handle.read().split())
        start = readme.index("Setup ends with exit status 0")
        sentence = readme[start:start + 900]
        self.assertNotIn("refused an answer or a file and changed nothing", sentence)
        self.assertIn("part-way", sentence)
        self.assertIn("changed nothing", sentence)       # said only of what is refused before anything is written

    def test_one_python_version_is_stated_and_the_code_does_not_need_a_newer_one(self):
        import ast
        with open(TBRC, encoding="utf-8") as handle:
            source = handle.read()
        ast.parse(source, feature_version=(3, 9))
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = handle.read()
        docstring = source.split('"""', 2)[1]
        mentioned = set(re.findall(r"Python (3\.\d+)", readme) + re.findall(r"Python (3\.\d+)", docstring))
        self.assertEqual(mentioned, {"3.9"}, "the README and tbrc must name one minimum version")

    def test_the_readmes_igblast_command_is_the_one_setup_runs(self):
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = handle.read()
        flat = " ".join(readme.replace("\\\n", " ").split())
        self.assertIn(self.cli.IGBLAST_BUILD, flat)
        self.assertIn("-w /work/igblast tbrc-backend:arm64 bash -c", flat)


# ====================================================================== "the whole lab" says what it means
class TestWholeLab(SetupCase):
    HOST = "labmac.example.org"

    def setUp(self):
        super().setUp()
        self.cli.socket.gethostname = lambda: self.HOST
        self.addCleanup(setattr, self.cli.socket, "gethostname", socket.gethostname)

    def lab(self, *flags):
        return self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web("lab")), *flags)

    def test_it_says_plain_http_unencrypted_one_account_and_how_to_get_https(self):
        code, out, err = self.lab()
        self.assertEqual(code, 0, out + err)
        flat = " ".join(out.split())
        for phrase in ("plain HTTP", "unencrypted", "share one account", "deploy/docker/Caddyfile", "HTTPS"):
            self.assertIn(phrase, flat)
        self.assertIn("any user id", flat, "anyone can submit under any user id")

    def test_it_prints_the_real_host_name_and_not_a_placeholder(self):
        code, out, err = self.lab()
        self.assertIn("http://%s:8080" % self.HOST, out)
        self.assertNotIn("<this computer", out)

    def test_the_note_comes_before_the_password_is_asked(self):
        events = []
        real_say = self.cli.say
        self.cli.say = lambda text="": (events.append(("say", text)), real_say(text))[1]
        self.terminal(["1", "1", "2", "2", "", "n"], hidden=[PASSWORD, PASSWORD])
        self.cli.getpass.getpass  # patched by terminal(); wrap to note the moment
        inner = self.cli.getpass.getpass
        self.cli.getpass.getpass = lambda *a, **k: (events.append(("password", "")), inner(*a, **k))[1]
        code, out, err = self.setup(None)
        self.assertEqual(code, 0, out + err)
        said = [i for i, (kind, text) in enumerate(events) if kind == "say" and "plain HTTP" in text]
        asked = [i for i, (kind, _) in enumerate(events) if kind == "password"]
        self.assertTrue(said and asked and said[0] < asked[0], events)

    def test_a_dry_run_says_it_too(self):
        code, out, err = self.lab("--dry-run")
        self.assertEqual(code, 0, out + err)
        self.assertIn("plain HTTP", " ".join(out.split()))

    def test_the_menu_says_what_the_choice_means(self):
        asked = self.terminal(["1", "1", "2", "1", "", "n"], hidden=[PASSWORD, PASSWORD])
        self.setup(None)
        menu = " ".join(" ".join(asked).split())
        self.assertIn("plain HTTP", menu)
        self.assertIn("one shared password", menu)

    def test_for_this_computer_only_none_of_it_is_said_and_the_address_is_localhost(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web("me")))
        self.assertEqual(code, 0, out + err)
        self.assertNotIn("plain HTTP", out)
        self.assertNotIn("unencrypted", out)
        self.assertIn("http://localhost:8080", out)
        self.assertNotIn(self.HOST, out)

    def test_a_computer_with_no_name_gets_a_clear_fallback(self):
        self.cli.socket.gethostname = lambda: ""
        code, out, err = self.lab()
        self.assertEqual(code, 0, out + err)
        self.assertIn("this computer's name", out)

    def test_the_readme_says_it_too(self):
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = " ".join(handle.read().split())
        for phrase in ("plain HTTP", "unencrypted", "one account", "deploy/docker/Caddyfile"):
            self.assertIn(phrase, readme)


# ====================================================================== the writable mounts of the web app
class TestWritableMounts(SetupCase):
    """The web container gets the workflow read-only (it is the part facing the network, and this checkout also
    feeds ./tbrc and deploy_all.sh), plus three writable folders. A mount point cannot be created inside a
    read-only mount, so setup makes the folders before the app starts."""
    FOLDERS = ("data/runs", "results", ".parallel")

    def exist(self):
        return [os.path.isdir(self.path("backend", *f.split("/"))) for f in self.FOLDERS]

    def test_the_web_app_gets_its_three_folders_made_and_before_it_starts(self):
        seen = []

        def run(cmd, **kwargs):
            self.calls.append({"cmd": list(cmd), "stdin": kwargs.get("stdin"), "stream": kwargs.get("stream"),
                               "cwd": kwargs.get("cwd")})
            if cmd[:2] == ["docker", "compose"]:
                seen.append(self.exist())
            return (0, FAKE_HASH + "\n", "") if "hash-password" in cmd else (0, "", "")
        self.cli.run = run
        self.assertEqual(self.exist(), [False, False, False])
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--start")
        self.assertEqual(code, 0, out + err)
        self.assertEqual(seen, [[True, True, True]], "the folders exist when compose is called")
        self.assertEqual(list(self.cli.WEB_WRITABLE_FOLDERS), list(self.FOLDERS))

    def test_they_are_made_when_the_commands_are_only_printed_too(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.exist(), [True, True, True])
        self.assertEqual(self.compose_calls(), [])

    def test_a_terminal_only_setup_makes_none(self):
        self.setup(LOCAL_UPLOAD_TERMINAL)
        self.assertEqual(self.exist(), [False, False, False])

    def test_a_dry_run_makes_none_and_says_it_would(self):
        code, out, err = self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()), "--dry-run")
        self.assertEqual(code, 0, out + err)
        self.assertEqual(self.exist(), [False, False, False])
        self.assertIn("would make backend/data/runs, backend/results, backend/.parallel", out)

    def test_existing_folders_and_their_contents_are_left_alone(self):
        os.makedirs(self.path("backend", "results", "u", "run1"))
        with open(self.path("backend", "results", "u", "run1", "keep.txt"), "w") as handle:
            handle.write("x")
        self.setup(dict(LOCAL_UPLOAD_TERMINAL, **self.web()))
        self.assertTrue(os.path.isfile(self.path("backend", "results", "u", "run1", "keep.txt")))


# ====================================================================== the files in deploy/docker
class TestDeployDockerFiles(unittest.TestCase):
    """The files a new user gets from git clone. They are read from the repository itself."""

    DOCKER = os.path.join(REPO, "deploy", "docker")

    def read(self, name):
        with open(os.path.join(self.DOCKER, name), encoding="utf-8") as handle:
            return handle.read()

    def test_the_files_are_there(self):
        for name in ("Caddyfile", "docker-compose.yml", "Dockerfile", "env.example", "seed-barcodes.sh",
                     os.path.join("ssh", "config"), ".dockerignore"):
            self.assertTrue(os.path.isfile(os.path.join(self.DOCKER, name)), name)
        self.assertTrue(os.access(os.path.join(self.DOCKER, "seed-barcodes.sh"), os.X_OK))

    def test_no_file_assumes_a_handoff_a_kit_a_saved_image_or_a_home_folder(self):
        for folder, _, files in os.walk(self.DOCKER):
            for name in files:
                path = os.path.join(folder, name)
                with open(path, encoding="utf-8") as handle:
                    text = handle.read().lower()
                for word in (OLD_KIT_FOLDER, "handoff", "docker load", HOME_FOLDERS, "tar.gz"):
                    self.assertNotIn(word, text, "%s mentions %r" % (os.path.relpath(path, REPO), word))
                self.assertIsNone(re.search(r"\bkit\b", text), os.path.relpath(path, REPO))

    def test_nothing_that_ships_names_the_old_kit_folder(self):
        paths = [os.path.join(REPO, "tbrc"), os.path.join(REPO, "README.md")]
        for top in ("tests", os.path.join("deploy", "docker")):
            for folder, _, files in os.walk(os.path.join(REPO, top)):
                paths += [os.path.join(folder, name) for name in files if not name.endswith(".pyc")]
        for path in paths:
            with open(path, encoding="utf-8") as handle:
                text = handle.read()
            self.assertNotIn(OLD_KIT_FOLDER, text, os.path.relpath(path, REPO))

    def test_the_dockerfile_has_the_two_targets_setup_builds(self):
        text = self.read("Dockerfile")
        self.assertRegex(text, r"(?m)^FROM \S+ AS base$")
        self.assertRegex(text, r"(?m)^FROM base AS backend$")
        self.assertRegex(text, r"(?m)^FROM backend AS unified$")

    def test_the_login_and_the_loopback_default_stay(self):
        self.assertIn("basic_auth", self.read("Caddyfile"))
        compose = self.read("docker-compose.yml")
        self.assertIn("${TBRC_BIND:-127.0.0.1}", compose)
        self.assertIn("TBRC_BIND=127.0.0.1", self.read("env.example").splitlines())
        self.assertIn("TBRC_PASSWORD_HASH=", self.read("env.example").splitlines())

    def test_the_workflow_is_mounted_read_only_with_three_writable_folders(self):
        text = self.read("docker-compose.yml")
        mounts = [m.group(1) for m in re.finditer(r"(?m)^\s*- (\.\./\.\./backend\S*)\s*$", text)]
        self.assertEqual(mounts, ["../../backend:/srv/backend:ro",
                                  "../../backend/data/runs:/srv/backend/data/runs",
                                  "../../backend/results:/srv/backend/results",
                                  "../../backend/.parallel:/srv/backend/.parallel"])
        self.assertNotIn("/srv/backend\n", text, "no writable mount of the whole backend")
        for mount in mounts[1:]:
            self.assertFalse(mount.endswith(":ro"), mount)

    def test_the_manual_route_makes_the_folders_first(self):
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = handle.read()
        self.assertIn("mkdir -p backend/data/runs backend/results backend/.parallel", readme)
        start = readme.index("### Starting the web app by hand")
        self.assertLess(readme.index("mkdir -p backend/data/runs backend/results backend/.parallel", start),
                        readme.index("docker compose --env-file .env up -d", start))

    def test_the_parallel_folder_is_ignored_by_git(self):
        with open(os.path.join(REPO, ".gitignore"), encoding="utf-8") as handle:
            self.assertIn("backend/.parallel/", handle.read().splitlines())

    def test_the_claim_about_where_it_was_tested_is_the_modest_one_everywhere(self):
        claim = "Tested on one Apple Silicon Mac; not tested on Linux or Intel"
        with open(os.path.join(REPO, "README.md"), encoding="utf-8") as handle:
            readme = handle.read()
        for name, text in (("README.md", readme), ("Dockerfile", self.read("Dockerfile")),
                           ("docker-compose.yml", self.read("docker-compose.yml"))):
            self.assertIn(claim, " ".join(text.replace("# ", " ").split()), name)
            for overclaim in ("Apple Silicon Macs only", "built and tested", "tested on linux/arm64"):
                self.assertNotIn(overclaim, text, "%s claims more than that: %r" % (name, overclaim))

    def test_the_context_sends_the_daemon_only_the_dockerfile(self):
        lines = [l for l in self.read(".dockerignore").splitlines() if l and not l.startswith("#")]
        self.assertEqual(lines, ["*", "!Dockerfile"], "the login hash in .env must not enter a build context")

    def test_the_env_file_is_ignored_by_git(self):
        with open(os.path.join(REPO, ".gitignore"), encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        for pattern in ("deploy/docker/.env", "deploy/docker/.env.*"):
            self.assertIn(pattern, lines)



if __name__ == "__main__":
    unittest.main()
