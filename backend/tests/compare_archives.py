#!/usr/bin/env python3
"""List what differs between two delivered archives (or two folders).

    python3 backend/tests/compare_archives.py A B

A and B are each a .zip, a .tar / .tar.gz / .tgz, or a folder. Files are matched by their path
inside the archive and compared byte for byte (SHA-256). It prints

  * the files that are in only one of the two,
  * the files that are in both but whose contents differ,

and nothing else: file dates, owners and the order of members are ignored, since a rerun always
changes them.

Exit status 0: the two hold the same files with the same contents. 1: they differ. 2: A or B could
not be read.

Plain standard library; runs on Python 3.6 and later.
"""
import hashlib
import os
import sys
import tarfile
import zipfile


class Unreadable(Exception):
    pass


def _digest(handle):
    sha = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1 << 20), b""):
        sha.update(chunk)
    return sha.hexdigest()


def _clean(name):
    """The member name without a leading ./, so 'a/x' and './a/x' are the same file."""
    while name.startswith("./"):
        name = name[2:]
    return name


def digests(path):
    """{path inside the archive: sha256 of its contents}, directories left out."""
    found = {}
    if os.path.isdir(path):
        for folder, _, names in os.walk(path):
            for name in names:
                full = os.path.join(folder, name)
                with open(full, "rb") as handle:
                    found[_clean(os.path.relpath(full, path))] = _digest(handle)
        return found
    if not os.path.isfile(path):
        raise Unreadable("%s does not exist" % path)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if not info.filename.endswith("/"):
                    with archive.open(info) as handle:
                        found[_clean(info.filename)] = _digest(handle)
        return found
    try:
        with tarfile.open(path) as archive:
            for member in archive:
                if member.isfile():
                    found[_clean(member.name)] = _digest(archive.extractfile(member))
    except (tarfile.TarError, OSError) as err:
        raise Unreadable("%s is not a zip, a tar archive or a folder (%s)" % (path, err))
    return found


def compare(first, second):
    """(only in first, only in second, in both but different), each sorted."""
    only_first = sorted(set(first) - set(second))
    only_second = sorted(set(second) - set(first))
    different = sorted(name for name in set(first) & set(second) if first[name] != second[name])
    return only_first, only_second, different


def main(argv):
    if len(argv) != 3:
        print("Usage: compare_archives.py A B", file=sys.stderr)
        return 2
    try:
        first = digests(argv[1])
        second = digests(argv[2])
    except Unreadable as err:
        print("compare_archives: %s" % err, file=sys.stderr)
        return 2
    only_first, only_second, different = compare(first, second)
    for title, names in (("Only in A", only_first), ("Only in B", only_second),
                         ("Contents differ", different)):
        if names:
            print("%s (%d):" % (title, len(names)))
            for name in names:
                print("  " + name)
    if only_first or only_second or different:
        return 1
    print("Same %d files, same contents." % len(first))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
