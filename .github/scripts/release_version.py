"""Manage two-component package versions and validate release distributions."""

import argparse
import hashlib
import json
import math
import re
import sys
import tarfile
import time
import tomllib
import urllib.error
import urllib.request
import zipfile
from email.parser import BytesParser
from pathlib import Path


VERSION = re.compile(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)")
INDEX_URLS = {
    "pypi": "https://pypi.org/pypi/svzerod/{version}/json",
    "testpypi": "https://test.pypi.org/pypi/svzerod/{version}/json",
}


def parse_version(version):
    match = VERSION.fullmatch(version) if isinstance(version, str) else None
    if match is None:
        raise ValueError(f"Expected a MAJOR.MINOR version without leading zeros: {version!r}")
    return tuple(int(part) for part in match.groups())


def read_version(path):
    document = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    project = document.get("project")
    if not isinstance(project, dict) or "version" not in project:
        raise ValueError("Missing [project].version in pyproject.toml")
    version = project["version"]
    parse_version(version)
    return version


def next_version(version, bump):
    major, minor = parse_version(version)
    if bump == "major":
        return f"{major + 1}.0"
    if bump == "minor":
        return f"{major}.{minor + 1}"
    raise ValueError(f"Unsupported version increment: {bump}")


def set_version(path, version, expected_current):
    parse_version(version)
    parse_version(expected_current)
    path = Path(path)
    current = read_version(path)
    if current != expected_current:
        raise ValueError(f"Expected current version {expected_current}, found {current}")
    if version == current:
        return
    text = path.read_bytes().decode("utf-8")
    project = re.search(r"(?m)^[ \t]*\[project\][ \t]*(?:#[^\r\n]*)?\r?$", text)
    if project is None:
        raise ValueError("Cannot locate the [project] table for a version update")
    following_table = re.search(r"(?m)^[ \t]*\[", text[project.end():])
    end = project.end() + following_table.start() if following_table else len(text)
    section = text[project.end():end]
    assignments = list(re.finditer(
        r"(?m)^[ \t]*version[ \t]*=[ \t]*([\"'])([^\r\n]*?)\1[ \t]*(?:#[^\r\n]*)?\r?$",
        section,
    ))
    if len(assignments) != 1 or assignments[0].group(2) != current:
        raise ValueError("Expected one single-line version assignment in [project]")
    start, stop = assignments[0].span(2)
    start += project.end()
    stop += project.end()
    updated = text[:start] + version + text[stop:]
    expected = tomllib.loads(text)
    expected["project"]["version"] = version
    if tomllib.loads(updated) != expected:
        raise ValueError("Version update would modify other TOML values")
    path.write_bytes(updated.encode("utf-8"))


def check_metadata(content, version, filename):
    metadata = BytesParser().parsebytes(content)
    for field, expected in (("Name", "svzerod"), ("Version", version)):
        if metadata.get_all(field, []) != [expected]:
            raise ValueError(f"{filename}: expected exactly one {field}: {expected}")


def inspect_distributions(dist, version):
    parse_version(version)
    paths = sorted(Path(dist).iterdir())
    kinds = set()
    prefix = f"svzerod-{version}"
    wheel_pattern = re.compile(
        re.escape(prefix) + r"-(?:[0-9][A-Za-z0-9_.]*-)?[A-Za-z0-9_.]+-[A-Za-z0-9_.]+-[A-Za-z0-9_.]+\.whl"
    )
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Unexpected distribution entry: {path.name}")
        if path.name.endswith(".whl"):
            if wheel_pattern.fullmatch(path.name) is None:
                raise ValueError(f"Unexpected wheel filename: {path.name}")
            with zipfile.ZipFile(path) as archive:
                members = [entry for entry in archive.infolist() if entry.filename.endswith(".dist-info/METADATA")]
                if len(members) != 1 or members[0].filename != f"{prefix}.dist-info/METADATA":
                    raise ValueError(f"{path.name}: expected one matching wheel METADATA file")
                check_metadata(archive.read(members[0]), version, path.name)
            kinds.add("wheel")
        elif path.name.endswith(".tar.gz"):
            if path.name != f"{prefix}.tar.gz":
                raise ValueError(f"Unexpected source filename: {path.name}")
            with tarfile.open(path, "r:gz") as archive:
                members = [entry for entry in archive.getmembers()
                           if entry.name.endswith("/PKG-INFO") and len(entry.name.split("/")) == 2]
                if len(members) != 1 or members[0].name != f"{prefix}/PKG-INFO" or not members[0].isfile():
                    raise ValueError(f"{path.name}: expected one matching source PKG-INFO file")
                check_metadata(archive.extractfile(members[0]).read(), version, path.name)
            kinds.add("source")
        else:
            raise ValueError(f"Unexpected distribution file: {path.name}")
    if kinds != {"source", "wheel"}:
        raise ValueError("Release distributions must include both a source archive and a wheel")
    return paths


def distribution_hashes(dist, version):
    hashes = {}
    for path in inspect_distributions(dist, version):
        with path.open("rb") as stream:
            hashes[path.name] = hashlib.file_digest(stream, "sha256").hexdigest()
    return hashes


def fetch_index(index, version):
    with urllib.request.urlopen(INDEX_URLS[index].format(version=version), timeout=30) as response:
        return json.load(response)


def compare_index(payload, hashes):
    if not isinstance(payload, dict) or not isinstance(payload.get("urls"), list):
        raise ValueError("Package index response has no distribution list")
    present = set()
    for entry in payload["urls"]:
        if not isinstance(entry, dict) or not isinstance(entry.get("filename"), str):
            raise ValueError("Package index response has an invalid distribution entry")
        filename = entry["filename"]
        if filename not in hashes:
            raise ValueError(f"Unexpected distribution already on the package index: {filename}")
        if filename in present:
            raise ValueError(f"Duplicate distribution on the package index: {filename}")
        digests = entry.get("digests")
        if not isinstance(digests, dict) or digests.get("sha256") != hashes[filename]:
            raise ValueError(f"SHA256 mismatch for distribution on the package index: {filename}")
        present.add(filename)
    return set(hashes) - present


def check_index(index, dist, version):
    hashes = distribution_hashes(dist, version)
    try:
        payload = fetch_index(index, version)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return
        raise
    compare_index(payload, hashes)


def verify_index(index, dist, version, attempts=12, delay=5):
    if attempts < 1 or not math.isfinite(delay) or delay < 0:
        raise ValueError("Index verification requires positive attempts and a finite nonnegative delay")
    hashes = distribution_hashes(dist, version)
    for attempt in range(attempts):
        try:
            payload = fetch_index(index, version)
        except urllib.error.HTTPError as error:
            if error.code != 404 and error.code < 500:
                raise
            problem = str(error)
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            problem = str(error)
        else:
            missing = compare_index(payload, hashes)
            if not missing:
                return
            problem = "Missing distributions on the package index: " + ", ".join(sorted(missing))
        if attempt + 1 < attempts:
            time.sleep(delay)
    raise ValueError(f"Could not verify the package index after {attempts} attempts: {problem}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("current", "next", "set"):
        subparser = commands.add_parser(command)
        subparser.add_argument("--path", type=Path, default=Path("pyproject.toml"))
        if command == "next":
            subparser.add_argument("--bump", choices=("minor", "major"), required=True)
        if command == "set":
            subparser.add_argument("--version", required=True)
            subparser.add_argument("--expected-current", required=True)
    for command in ("verify-dist", "check-index", "verify-index"):
        subparser = commands.add_parser(command)
        subparser.add_argument("--dist", type=Path, default=Path("dist"))
        subparser.add_argument("--version", required=True)
        if command != "verify-dist":
            subparser.add_argument("--index", choices=tuple(INDEX_URLS), required=True)
        if command == "verify-index":
            subparser.add_argument("--attempts", type=int, default=12)
            subparser.add_argument("--delay", type=float, default=5)
    args = parser.parse_args()
    try:
        if args.command == "current":
            print(read_version(args.path))
        elif args.command == "next":
            print(next_version(read_version(args.path), args.bump))
        elif args.command == "set":
            set_version(args.path, args.version, args.expected_current)
            print(args.version)
        elif args.command == "verify-dist":
            paths = inspect_distributions(args.dist, args.version)
            print(f"Verified {len(paths)} distributions for svzerod {args.version}")
        elif args.command == "check-index":
            check_index(args.index, args.dist, args.version)
            print(f"Existing {args.index} files match the candidate distributions")
        else:
            verify_index(args.index, args.dist, args.version, args.attempts, args.delay)
            print(f"Verified all distributions on {args.index} for svzerod {args.version}")
    except (OSError, ValueError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
