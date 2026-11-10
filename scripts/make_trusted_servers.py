#!/usr/bin/env python3
"""Write the list of servers and program paths the cluster trusts.

    make_trusted_servers.py            # rewrite backend/config/trusted_servers.json
    make_trusted_servers.py --check    # exit 1 if that file is not what the config would give

The cluster reads backend/config/trusted_servers.json to refuse a sample.json whose destination or
program paths were not set by this installation: a value that is not on the list is not run. The list
is generated from frontend/server_config.json, which stays the one place to edit. After changing a
server, a path prefix or an IgBLAST path there, run this script and commit both files; deploy_all.sh
runs it with --check and stops when they have drifted apart.

Only the fields a sample.json carries are copied: web ports, e-mail addresses and the rest of the
config never leave the front end. Each IgBLAST setting is a list of accepted strings: today's value,
then the older install locations that finished runs still carry and rerun_final.sh reuses.

Python standard library only.
"""
import argparse
import json
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DEFAULT_CONFIG = os.path.join(REPO, "frontend", "server_config.json")
DEFAULT_OUT = os.path.join(REPO, "backend", "config", "trusted_servers.json")

# (field, type); every storage_server entry must be an object with all of these
SERVER_FIELDS = (("rawdata_path", str), ("results_path", str), ("ssh_port", int), ("transfer_mode", str))
CLUSTER_FIELDS = ("path_prefix", "snakemake_env")
IGBLAST_FIELDS = ("igblast_bin", "igblast_data", "igblast_auxiliary_data",
                  "igblast_db_v", "igblast_db_d", "igblast_db_j")
# no database prefix set is the empty string
IGBLAST_MAY_BE_EMPTY = ("igblast_db_v", "igblast_db_d", "igblast_db_j")

# Fields whose empty value is also accepted. A sample.json from `tbrc run --local`, `tbrc run
# --fastq-dir` or benchmark/run_benchmark.sh leaves these null, and run_igblast.sh then uses the
# copies bundled in backend/igblast, which is no outside path. Listed after the legacy paths.
IGBLAST_ALSO_EMPTY = ("igblast_data", "igblast_auxiliary_data")

# Install locations earlier than today's that real past runs still carry in their sample.json and
# that rerun_final.sh reuses. Appended after the configured value, never before it.
IGBLAST_LEGACY = {
    "igblast_data": [
        "/path/to/home/Downloads/ncbi-igblast-1.22.0/internal_data",
        "/path/to/scratch/pipelines/TBRC/igblast/internal_data",
    ],
    "igblast_auxiliary_data": [
        "/path/to/home/Downloads/ncbi-igblast-1.22.0/optional_file",
        "/path/to/scratch/pipelines/TBRC/igblast/optional_file",
    ],
    # The IgH D database prefix that two finished runs (2026-11-08 check of the stored sample.json
    # files) still carry; the current config sets none.
    "igblast_db_d": [
        "/path/to/scratch/pipelines/TBRC/igblast/work/mouse/igh/mouse_igh_D",
    ],
}


class ConfigError(Exception):
    """The config cannot be turned into a trusted list; the message says what to fix."""


def _require_str(container, field, where, allow_empty=False):
    value = container.get(field) if isinstance(container, dict) else None
    if not isinstance(value, str) or (not value and not allow_empty):
        raise ConfigError("%s.%s must be %s string" % (where, field, "a" if allow_empty else "a non-empty"))
    return value


def _server_entry(name, entry):
    where = "storage_server.%s" % name
    if not isinstance(entry, dict):
        raise ConfigError(
            "%s is %s, not an object. Use the object form with rawdata_path, results_path, "
            "ssh_port and transfer_mode" % (where, "a string (the older form)" if isinstance(entry, str)
                                            else type(entry).__name__))
    out = {}
    for field, kind in SERVER_FIELDS:
        if field not in entry:
            raise ConfigError("%s.%s is missing" % (where, field))
        value = entry[field]
        if kind is int:
            if isinstance(value, bool) or not isinstance(value, int) or not 0 < value < 65536:
                raise ConfigError("%s.%s must be a port number (an integer from 1 to 65535)" % (where, field))
        else:
            _require_str(entry, field, where)
        out[field] = value
    return out


def build(config):
    """The trusted list for a parsed server_config.json."""
    if not isinstance(config, dict):
        raise ConfigError("the config must be a JSON object")
    storage = config.get("storage_server")
    if not isinstance(storage, dict) or not storage:
        raise ConfigError("storage_server must list at least one server")
    cluster = config.get("cluster")
    if not isinstance(cluster, dict):
        raise ConfigError("cluster must be an object")

    servers = {name: _server_entry(name, entry) for name, entry in storage.items()}

    trusted_cluster = {f: _require_str(cluster, f, "cluster") for f in CLUSTER_FIELDS}

    igblast = {}
    for field in IGBLAST_FIELDS:
        current = _require_str(cluster, field, "cluster", allow_empty=field in IGBLAST_MAY_BE_EMPTY)
        accepted = [current]
        for legacy in IGBLAST_LEGACY.get(field, []):
            if legacy not in accepted:
                accepted.append(legacy)
        if field in IGBLAST_ALSO_EMPTY and "" not in accepted:
            accepted.append("")
        igblast[field] = accepted

    return {"servers": servers, "cluster": trusted_cluster, "igblast": igblast}


def render(trusted):
    return json.dumps(trusted, indent=2, sort_keys=True) + "\n"


def load_config(path):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except OSError as err:
        raise ConfigError("cannot read %s: %s" % (path, err.strerror or err))
    except ValueError as err:
        raise ConfigError("%s is not valid JSON: %s" % (path, err))


def write_atomically(path, text):
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    tmp = os.path.join(folder, ".%s.tmp%d" % (os.path.basename(path), os.getpid()))
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="write nothing; exit 1 if the trusted list is missing or out of step with the config")
    ap.add_argument("--config", default=DEFAULT_CONFIG, help="server config (default: frontend/server_config.json)")
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="trusted list (default: backend/config/trusted_servers.json)")
    args = ap.parse_args(argv)

    try:
        text = render(build(load_config(args.config)))
    except ConfigError as err:
        print("make_trusted_servers.py: %s" % err, file=sys.stderr)
        return 2

    if args.check:
        try:
            with open(args.out, encoding="utf-8") as f:
                current = f.read()
        except OSError:
            current = None
        if current == text:
            return 0
        print("make_trusted_servers.py: %s is %s %s. Run: python3 scripts/make_trusted_servers.py "
              "and commit the result." % (args.out, "missing, so it is not in step with" if current is None
                                          else "out of step with", args.config), file=sys.stderr)
        return 1

    write_atomically(args.out, text)
    print("wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
