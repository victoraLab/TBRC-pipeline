"""What a sample.json row may hold, checked when the Snakefile loads it and before any job starts.

The web container that writes sample.json might be hijacked, and it can still hand the cluster a
doctored file and launch it. So the Snakefile calls validate_sample() on every row, and a row with
one bad value stops the run. This is the cluster's copy of the form's rules (frontend/
validate_input.R, sections 1 and 2 of docs/superpowers/specs/2026-11-08-tbrc-shell-injection-
fix-design.md); the choice sets below are a second copy on purpose and backend/tests/
test_sample_rules.py keeps the two equal.

Two kinds of field need more than a good format:

  * where results go: server_results_path, server_ssh_port, server_transfer_mode;
  * which program runs: igblast_bin, igblast_data, igblast_auxiliary_data, igblast_db_v/d/j.

These must equal what backend/config/trusted_servers.json lists for the run's server_name, because
a compute node pushes with the cluster account's own key. A deployment with no such file keeps the
format checks only, until its turn comes.

A value is right as written or refused: nothing is trimmed or case-folded here, since the form
already cleaned it before writing. A refusal names the field and what is wrong, never the value.

Plain standard library, and nothing newer than Python 3.6, so it runs wherever the Snakefile does.
"""
import json
import math
import numbers
import os
import re

# ---------------------------------------------------------------- choice sets
# The form's lists (frontend/validate_input.R, igblast_panel_choices in frontend/global.R) and, where
# the pipeline accepts more than the form offers, the old value that earlier sample.json still carry.
METHODS = ("classic_luka", "tiago_vsearch", "tiago_nanopore", "tiago_sulkseq")
ARCHIVE_SCOPES = ("light", "full")
ARCHIVE_FORMATS = ("zip", "tar.gz")
SPECIES = ("human", "mouse")
RECEPTOR_CHECKS = ("auto", "warn", "off")
LUKA_LIGHT_PROTOCOLS = ("standard", "luka_light")
MIBS_GENE_GROUPINGS = ("strict", "transitive")
MIBS_CLONE_SCOPES = ("sample", "global")
INPUT_METHODS = ("Server", "Upload")
TRANSFER_MODES = ("push", "pull", "none")
# ig_light is no longer offered but older sample.json carry it and the backend still runs it.
IGBLAST_PANELS = ("igh", "ig_all", "tra", "trb", "tcr_all", "all_receptors", "ig_b6", "ig_light")

# Not in the trusted list: the CLI's local runs and the form's Upload mode. Nothing is pushed from
# either, so their destination must be empty. `tbrc run --fastq-dir` writes the empty string when
# there is no storage server, and that is a local run too.
BARE_SERVERS = ("local", "upload", "")

TRIM_ALPHABET = "ACGTRYSWKMBDHVN"
TRIM_MIN = 10
TRIM_MAX = 64

# sequence_number is twice the form's n_read (1 to 500); older files hold the unmultiplied count.
SEQUENCE_NUMBER_RANGE = (1, 1000)

# ---------------------------------------------------------------- patterns
# re.fullmatch only: "$" would let a final newline through.
_SEGMENT = re.compile(r"[A-Za-z0-9._-]+")
_TRIM = re.compile(r"[%s]{%d,%d}" % (TRIM_ALPHABET, TRIM_MIN, TRIM_MAX))
# The local part may not start with "-": an address that did would read as an option if it were ever
# put on a mail program's command line. Keep in step with is_valid_email() in frontend/global.R and the
# check in scripts/final.script.sh.
_EMAIL = re.compile(r"[A-Za-z0-9._%+][A-Za-z0-9._%+-]*@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
_PATH = re.compile(r"[A-Za-z0-9._/@:+-]+")        # user@host:/dir/ or /dir/
_PROGRAM = re.compile(r"[A-Za-z0-9._/+-]+")       # igblastn or /dir/igblastn
_PORT_TEXT = re.compile(r"[0-9]{1,5}")
_BOOL_TEXT = ("TRUE", "FALSE", "true", "false")

_SEGMENT_HELP = ("must be one path segment: letters, digits, dot, underscore and hyphen only, "
                 "not . or .. and not starting with -")
_PATH_HELP = ("must be a plain path: letters, digits and . _ / @ : + - only, with no .. and "
              "nothing starting with -")


class SampleRejected(ValueError):
    """sample.json holds a value the pipeline will not run. .problems has one line per field."""

    def __init__(self, problems):
        self.problems = list(problems)
        ValueError.__init__(self, "sample.json refused, %d problem%s:\n  %s" % (
            len(self.problems), "" if len(self.problems) == 1 else "s", "\n  ".join(self.problems)))


# ---------------------------------------------------------------- value helpers
def _is_nan(value):
    if isinstance(value, float):
        return value != value
    return type(value).__name__ in ("NAType", "NaTType")


def _is_blank(value, texts=("",)):
    """None, NaN, or one of the strings that stand for nothing."""
    if value is None or _is_nan(value):
        return True
    return isinstance(value, str) and value in texts


def _is_bool(value):
    return isinstance(value, bool) or type(value).__name__ == "bool_"


def _truth(value):
    """True or False for a real boolean or TRUE/FALSE/true/false; None for anything else."""
    if _is_bool(value):
        return bool(value)
    if isinstance(value, str) and value in _BOOL_TEXT:
        return value.upper() == "TRUE"
    return None


def _whole(value):
    """An integer for an integer or an integral float (a table with a gap turns ints into floats)."""
    if _is_bool(value):
        return None
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, float) and math.isfinite(value) and value == int(value):
        return int(value)
    return None


def _segment_ok(value):
    return (isinstance(value, str) and _SEGMENT.fullmatch(value) is not None
            and value not in (".", "..") and not value.startswith("-"))


def _path_ok(value, pattern=_PATH):
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        return False
    parts = re.split(r"[/:@]", value)
    return not any(part in (".", "..") or part.startswith("-") for part in parts)


def _port_number(value):
    if isinstance(value, str):
        return int(value) if _PORT_TEXT.fullmatch(value) else None
    return _whole(value)


# ---------------------------------------------------------------- the rules
def _rule(kind, present=False, blank=True, **more):
    """present: the Snakefile indexes this column unconditionally, so a missing one is a problem.
    blank: None, NaN or an empty string passes. blank_text: the strings that count as empty."""
    rule = {"kind": kind, "present": present, "blank": blank}
    rule.update(more)
    return rule


_MIBS_NUMBER = dict(kind="integer", present=False, blank=True, blank_text=("", "NA"))

FIELD_RULES = {
    "user_id": _rule("segment", present=True, blank=False),
    "user_email": _rule("email", present=True),
    "sample_name": _rule("segment"),
    "input_method": _rule("choice", blank=False, choices=INPUT_METHODS),
    "server_name": _rule("server_name", present=True, blank=False),
    "server_path": _rule("path", present=True, blank_text=("", "NA")),
    "server_results_path": _rule("results_path", present=True, blank=False),
    "server_ssh_port": _rule("port", present=True, blank=False),
    "server_transfer_mode": _rule("transfer_mode", present=True, blank=False),
    "primer_set": _rule("segment", present=True, blank=False),
    "trim_primers": _rule("boolean", present=True, blank=False),
    "sequence_number": _rule("integer", present=True, blank=False, range=SEQUENCE_NUMBER_RANGE),
    # An absent, null or NaN pipeline takes the default, an empty string does not.
    "pipeline": _rule("choice", choices=METHODS, blank_text=()),
    "luka_light_protocol": _rule("choice", present=True, blank=False, choices=LUKA_LIGHT_PROTOCOLS),
    "run_igblast": _rule("boolean", present=True, blank=False),
    "run_clonality": _rule("boolean", present=True, blank=False),
    "archive_format": _rule("choice", choices=ARCHIVE_FORMATS),
    "igblast_species": _rule("choice", present=True, blank=False, choices=SPECIES),
    "igblast_panel": _rule("choice", present=True, blank=False, choices=IGBLAST_PANELS),
    "mibs_junction_max_distance": dict(_MIBS_NUMBER, range=(0, 30)),
    "mibs_gene_grouping": _rule("choice", choices=MIBS_GENE_GROUPINGS),
    "mibs_clone_scope": _rule("choice", choices=MIBS_CLONE_SCOPES),
    "mibs_min_reads_per_well": dict(_MIBS_NUMBER, range=(1, 10000)),
    "mibs_min_wells_per_clone": dict(_MIBS_NUMBER, range=(1, 96)),
    "mibs_productive_only": _rule("boolean"),
    "mibs_min_reads_per_cluster": dict(_MIBS_NUMBER, range=(1, 1000)),
    "platemap_attached": _rule("boolean"),
    "sequencing_kit": _rule("segment"),
    "receptor_check": _rule("choice", choices=RECEPTOR_CHECKS),
    "igblast_bin": _rule("igblast_list", present=True, blank=False, pattern=_PROGRAM),
    "igblast_organism": _rule("choice", present=True, blank=False, choices=SPECIES),
    "igblast_db_v": _rule("igblast_list", present=True, pattern=_PATH),
    "igblast_db_d": _rule("igblast_list", present=True, pattern=_PATH),
    "igblast_db_j": _rule("igblast_list", present=True, pattern=_PATH),
    "igblast_data": _rule("igblast_list", pattern=_PATH),
    "igblast_auxiliary_data": _rule("igblast_list", present=True, pattern=_PATH),
    "trim_seq": _rule("trim_sequence", present=True, blank=False),
    "keep_raw": _rule("boolean", present=True, blank=False),
    "archive_scope": _rule("choice", choices=ARCHIVE_SCOPES),
    "run_group": _rule("segment"),
    "source": _rule("source", present=True, blank=False),
    "dest": _rule("dest", present=True, blank=False),
}


# ---------------------------------------------------------------- one check per kind
# Each takes (field, rule, value, row, trusted) and returns None, or what is wrong, without the value.
def _plain(check):
    """The kinds that look at one value alone: a blank passes or not as the rule says."""
    def run(field, rule, value, row, trusted):
        if _is_blank(value, rule.get("blank_text", ("",))):
            return None if rule["blank"] else "is empty"
        return check(rule, value)
    return run


def _segment(rule, value):
    return None if _segment_ok(value) else _SEGMENT_HELP


def _choice(rule, value):
    if isinstance(value, str) and value in rule["choices"]:
        return None
    return "must be one of: " + ", ".join(rule["choices"])


def _boolean(rule, value):
    return None if _truth(value) is not None else "must be true or false (TRUE, FALSE, true, false)"


def _integer(rule, value):
    low, high = rule["range"]
    number = _whole(value)
    if number is not None and low <= number <= high:
        return None
    return "must be a whole number from %d to %d" % (low, high)


def _email(rule, value):
    if isinstance(value, str) and len(value) <= 254 and _EMAIL.fullmatch(value):
        return None
    return "must be an email address, or empty"


def _path(rule, value):
    return None if _path_ok(value) else _PATH_HELP


def _trim_sequence(field, rule, value, row, trusted):
    if _is_blank(value):
        if _truth(row.get("trim_primers")) is False:
            return None
        return "is empty, but Trim primers is on"
    if isinstance(value, str) and _TRIM.fullmatch(value):
        return None
    return "must be %d to %d letters from %s" % (TRIM_MIN, TRIM_MAX, " ".join(TRIM_ALPHABET))


def _run_parts(value, counts):
    """The parts of a results/<user>/<run> style path, or None if it is not one."""
    if not isinstance(value, str):
        return None
    parts = value.split("/")
    if len(parts) not in counts:
        return None
    return parts


def _user_of(row):
    user = row.get("user_id")
    return user if _segment_ok(user) else None


def _dest(field, rule, value, row, trusted):
    parts = _run_parts(value, (3,))
    if parts is None or parts[0] != "results" or not all(_segment_ok(p) for p in parts[1:]):
        return "must be results/<user_id>/<run>, each part one path segment"
    user = _user_of(row)
    if user is not None and parts[1] != user:
        return "must be results/<user_id>/<run>, and its user is not user_id"
    return None


def _source(field, rule, value, row, trusted):
    parts = _run_parts(value, (1, 2))
    if parts is None or not all(_segment_ok(p) for p in parts):
        return "must be <user_id>/<run>, or one path segment, each part a path segment"
    user = _user_of(row)
    if len(parts) == 2 and user is not None and parts[0] != user:
        return "must be <user_id>/<run>, and its user is not user_id"
    return None


# -- where the results go
def _server_kind(row, trusted):
    """('trusted', entry) for a server in the list, ('bare', None) for local and upload, ('free', None)
    for any plain name when there is no list, (None, None) for a name that is none of these."""
    name = row.get("server_name")
    if isinstance(name, str):
        if trusted is not None and name in trusted["servers"]:
            return "trusted", trusted["servers"][name]
        if name in BARE_SERVERS:
            return "bare", None
        if trusted is None and _segment_ok(name):
            return "free", None
    return None, None


def _empty_destination_ok(row):
    """An empty server_results_path is accepted only with transfer mode pull or none."""
    return row.get("server_transfer_mode") in ("pull", "none")


def _pushes(row):
    """True when final.script.sh would rsync to the destination: it does so for every mode but pull
    (including none) as soon as server_results_path is not empty or NA."""
    return (row.get("server_transfer_mode") != "pull"
            and not _is_blank(row.get("server_results_path"), ("", "NA")))


def _server_name(field, rule, value, row, trusted):
    if _server_kind(row, trusted)[0] is not None:
        return None
    if trusted is None:
        return "must be local, upload, empty, or one path segment"
    return "is not local, upload, empty or a server in the trusted list"


def _results_path(field, rule, value, row, trusted):
    kind, entry = _server_kind(row, trusted)
    empty = _is_blank(value, ("", "NA"))
    if kind == "bare":
        return None if empty else "must be empty for a local run or an upload, which push nothing"
    if kind == "trusted":
        if empty:
            if _empty_destination_ok(row):
                return None
            return "is empty, which is accepted only with transfer mode pull or none"
        if isinstance(value, str) and value == entry["results_path"]:
            return None
        return "is not the destination trusted for this server"
    return None if empty or _path_ok(value) else _PATH_HELP


def _port(field, rule, value, row, trusted):
    kind, entry = _server_kind(row, trusted)
    if _is_blank(value, ("", "NA")):
        if kind == "trusted" and _pushes(row):
            return "is empty, but a push needs the trusted port"
        return None
    number = _port_number(value)
    if number is None or not 0 < number < 65536:
        return "must be a port number from 1 to 65535"
    if kind == "trusted" and _pushes(row) and number != _whole(entry["ssh_port"]):
        return "is not the port trusted for this server"
    return None


def _transfer_mode(field, rule, value, row, trusted):
    if not isinstance(value, str) or value not in TRANSFER_MODES:
        return "must be one of: " + ", ".join(TRANSFER_MODES)
    kind, entry = _server_kind(row, trusted)
    if kind == "trusted" and value not in ("pull", "none") and value != entry["transfer_mode"]:
        return "is not the transfer mode trusted for this server"
    return None


# -- which program runs
def _igblast_list(field, rule, value, row, trusted):
    if trusted is not None:
        accepted = trusted["igblast"].get(field) or []
        text = "" if _is_blank(value) else value
        if isinstance(text, str) and text in accepted:
            return None
        return "is not one of the values trusted on this cluster"
    if _is_blank(value):
        return None if rule["blank"] else "is empty"
    return None if _path_ok(value, rule["pattern"]) else _PATH_HELP


_CHECKS = {
    "segment": _plain(_segment),
    "choice": _plain(_choice),
    "boolean": _plain(_boolean),
    "integer": _plain(_integer),
    "email": _plain(_email),
    "path": _plain(_path),
    "trim_sequence": _trim_sequence,
    "dest": _dest,
    "source": _source,
    "server_name": _server_name,
    "results_path": _results_path,
    "port": _port,
    "transfer_mode": _transfer_mode,
    "igblast_list": _igblast_list,
}


def validate_sample(row, trusted):
    """Raise SampleRejected, listing every problem at once, unless `row` may run.

    `row` is one record of the sample table (samples_df.to_dict("records")). `trusted` is what
    load_trusted() returned: None means the format checks only."""
    if not isinstance(row, dict):
        raise SampleRejected(["sample: a row of fields was expected"])
    problems = []
    for field, rule in FIELD_RULES.items():
        if field not in row:
            if rule["present"]:
                problems.append("%s: is missing, and the workflow reads it" % field)
            continue
        message = _CHECKS[rule["kind"]](field, rule, row[field], row, trusted)
        if message:
            problems.append("%s: %s" % (field, message))
    if problems:
        raise SampleRejected(problems)


def validate_sample_name(name):
    """Raise SampleRejected unless `name` is one path segment.

    The Snakefile takes the run's name from the sample_name field, which validate_sample() checks,
    but when that is blank it derives it from the FASTQ file names. Those come from disk, not from
    sample.json, and the name ends up in every output path and shell command, so it gets the same
    segment rule here. The message does not repeat the name."""
    if not _segment_ok(name):
        raise SampleRejected(["sample_name: the name the run is called by (taken from the FASTQ file "
                              "names when sample.json has none) " + _SEGMENT_HELP])


def _run_folder_of(table_path):
    """(user, run) when `table_path` is a file directly inside data/runs/<user>/<run>/, else None.

    Looks at the end of the path only, so it matches however the pipeline is handed the table: relative
    to the root or to .parallel/<key> (whose data/ is a symlink to the root's), absolute on the cluster
    or inside a container. It reads the text; nothing on disk is touched."""
    if not isinstance(table_path, str) or not table_path:
        return None
    parts = os.path.normpath(table_path).split(os.sep)
    if len(parts) < 5 or parts[-5:-3] != ["data", "runs"]:
        return None
    return parts[-3], parts[-2]


def validate_run_folder(rows, table_path):
    """Raise SampleRejected unless every row belongs to the run folder its table sits in.

    The launchers put a run's sample.json at data/runs/<user>/<run>/, and the form builds the row from
    the same two names: user_id is <user>, dest is results/<user>/<run>, source is <user>/<run>. dest
    drives the cleanup (rm -rf) in final.script.sh, and the web app's key can put any file in any run
    folder, so a row that names another run is refused. dest and user_id are exact. source may also be
    the older one-segment <run>, which sample.json files from before the form added the user still
    hold; no rule reads it, and a re-run of such a run must keep working. A table anywhere else (a
    test fixture, a run outside data/runs) is left alone. The message names the fields and never
    repeats a value."""
    folder = _run_folder_of(table_path)
    if folder is None:
        return
    user, run = folder
    wanted = (("user_id", (user,), "must be the <user> of"),
              ("dest", ("results/%s/%s" % (user, run),), "must be results/<user>/<run> for"),
              ("source", ("%s/%s" % (user, run), run), "must be <user>/<run>, or just <run>, for"))
    problems = []
    for row in rows:
        if not isinstance(row, dict):
            problems.append("sample: a row of fields was expected")
            continue
        for field, accepted, what in wanted:
            value = row.get(field)
            if not (isinstance(value, str) and value in accepted):
                problems.append("%s: %s the run folder (data/runs/<user>/<run>/) that holds this sample.json"
                                % (field, what))
    if problems:
        raise SampleRejected(problems)


# ---------------------------------------------------------------- the trusted list
def load_trusted(path):
    """The cluster's trusted list, or None when the file is absent (the format checks then stand
    alone). A file that is there but cannot be used stops the run: silently dropping to format
    checks would switch the destination check off."""
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError) as err:
        raise SampleRejected(["%s: cannot be read as JSON (%s)" % (path, err)])
    problems = []
    if not isinstance(data, dict):
        raise SampleRejected(["%s: expected a JSON object" % path])

    servers = data.get("servers")
    if not isinstance(servers, dict):
        problems.append("%s: servers must be an object" % path)
    else:
        for name, entry in sorted(servers.items()):
            where = "%s: servers.%s" % (path, name)
            if not isinstance(entry, dict):
                problems.append(where + " must be an object")
                continue
            if not _path_ok(entry.get("results_path")):
                problems.append(where + ".results_path must be a plain non-empty path")
            port = entry.get("ssh_port")
            if _is_bool(port) or not isinstance(port, int) or not 0 < port < 65536:
                problems.append(where + ".ssh_port must be an integer port number")
            if entry.get("transfer_mode") not in TRANSFER_MODES:
                problems.append(where + ".transfer_mode must be one of " + ", ".join(TRANSFER_MODES))

    igblast = data.get("igblast")
    if not isinstance(igblast, dict):
        problems.append("%s: igblast must be an object" % path)
    else:
        for field, rule in sorted(FIELD_RULES.items()):
            if rule["kind"] != "igblast_list":
                continue
            values = igblast.get(field)
            if not (isinstance(values, list) and values and all(isinstance(v, str) for v in values)):
                problems.append("%s: igblast.%s must be a non-empty list of strings" % (path, field))
    if problems:
        raise SampleRejected(problems)
    return data
