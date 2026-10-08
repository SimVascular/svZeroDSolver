import importlib.util
import io
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
import urllib.error
import zipfile
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "release_version.py"


class ReleaseVersionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.project = self.root / "pyproject.toml"
        self.project.write_text(
            '[project]\nname = "pysvzerod"\nversion = "3.1" # current\n'
            '\n[tool.example]\nversion = "9.8.7"\n', encoding="utf-8"
        )

    def cli(self, *arguments):
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            cwd=self.root, capture_output=True, text=True, check=False,
        )

    def helper(self):
        spec = importlib.util.spec_from_file_location("release_version", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def distributions(self, *, version="3.2", name="pysvzerod"):
        dist = self.root / "dist"
        dist.mkdir(exist_ok=True)
        metadata = f"Metadata-Version: 2.4\nName: {name}\nVersion: {version}\n\n".encode()
        wheel = dist / "pysvzerod-3.2-cp311-cp311-linux_x86_64.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("pysvzerod-3.2.dist-info/METADATA", metadata)
        source = dist / "pysvzerod-3.2.tar.gz"
        with tarfile.open(source, "w:gz") as archive:
            for member_name in (
                "pysvzerod-3.2/PKG-INFO",
                "pysvzerod-3.2/pysvzerod.egg-info/PKG-INFO",
            ):
                member = tarfile.TarInfo(member_name)
                member.size = len(metadata)
                archive.addfile(member, io.BytesIO(metadata))
        return dist, wheel, source

    def remote_files(self, files):
        import hashlib

        return {"urls": [
            {"filename": path.name, "digests": {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()
            }} for path in files
        ]}

    def test_current_reads_project_version(self):
        result = self.cli("current")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "3.1")

    def test_next_supports_only_major_and_minor(self):
        for bump, expected in (("minor", "3.2"), ("major", "4.0")):
            with self.subTest(bump=bump):
                result = self.cli("next", "--bump", bump, "--path", str(self.project))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), expected)
        self.assertNotEqual(self.cli("next", "--bump", "patch").returncode, 0)

    def test_current_rejects_malformed_versions_and_toml(self):
        for version in ("3.1.2", "03.1", "3.01", "v3.1", "3", "3.2rc1"):
            with self.subTest(version=version):
                self.project.write_text(f'[project]\nversion = "{version}"\n')
                self.assertNotEqual(self.cli("current").returncode, 0)
        self.project.write_text('[project]\nversion = "3.1"\nversion = "3.2"\n')
        self.assertNotEqual(self.cli("current").returncode, 0)

    def test_set_preserves_comments_other_versions_and_newlines(self):
        original = self.project.read_bytes().replace(b'"3.1"', b"'3.1'").replace(b"\n", b"\r\n")
        self.project.write_bytes(original)
        result = self.cli("set", "--version", "4.0", "--expected-current", "3.1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.project.read_bytes(), original.replace(b"'3.1'", b"'4.0'"))

    def test_set_rejects_stale_version_without_changing_file(self):
        original = self.project.read_bytes()
        result = self.cli("set", "--version", "4.0", "--expected-current", "3.0")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.project.read_bytes(), original)

    def test_set_current_version_is_noop(self):
        original = self.project.read_bytes()
        result = self.cli("set", "--version", "3.1", "--expected-current", "3.1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.project.read_bytes(), original)

    def test_set_rejects_patch_version(self):
        original = self.project.read_bytes()
        self.assertNotEqual(self.cli("set", "--version", "3.1.1", "--expected-current", "3.1").returncode, 0)
        self.assertEqual(self.project.read_bytes(), original)

    def test_set_does_not_change_embedded_table_text(self):
        original = (
            'description = """\n[project]\nversion = "3.1"\n"""\n'
            '[project]\nname = "pysvzerod"\nversion = "3.1"\n'
        ).encode()
        self.project.write_bytes(original)
        result = self.cli("set", "--version", "3.2", "--expected-current", "3.1")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.project.read_bytes(), original)

    def test_verify_distribution_metadata(self):
        dist, wheel, _ = self.distributions()
        wheel.rename(wheel.with_name("pysvzerod-3.2-1build-cp311-cp311-linux_x86_64.whl"))
        result = self.cli("verify-dist", "--dist", str(dist), "--version", "3.2")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_requires_source_and_wheel(self):
        dist, wheel, source = self.distributions()
        source.unlink()
        self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)
        wheel.unlink()
        self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_rejects_wrong_metadata_name_and_version(self):
        for version, name in (("3.1", "pysvzerod"), ("3.2", "another-project")):
            with self.subTest(version=version, name=name):
                dist, _, _ = self.distributions(version=version, name=name)
                self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_rejects_mismatched_filename(self):
        dist, wheel, _ = self.distributions()
        wheel.rename(wheel.with_name("pysvzerod-3.3-cp311-cp311-linux_x86_64.whl"))
        self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_rejects_missing_or_duplicate_wheel_metadata(self):
        dist, wheel, _ = self.distributions()
        for members in ([], ["pysvzerod-3.2.dist-info/METADATA", "other-3.2.dist-info/METADATA"]):
            with self.subTest(members=members):
                with zipfile.ZipFile(wheel, "w") as archive:
                    for name in members:
                        archive.writestr(name, "Name: pysvzerod\nVersion: 3.2\n")
                self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_rejects_missing_or_duplicate_source_metadata(self):
        dist, _, source = self.distributions()
        for count in (0, 2):
            with self.subTest(count=count):
                with tarfile.open(source, "w:gz") as archive:
                    for _ in range(count):
                        content = b"Name: pysvzerod\nVersion: 3.2\n"
                        member = tarfile.TarInfo("pysvzerod-3.2/PKG-INFO")
                        member.size = len(content)
                        archive.addfile(member, io.BytesIO(content))
                self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_rejects_source_version_different_from_wheel(self):
        dist, _, source = self.distributions()
        with tarfile.open(source, "w:gz") as archive:
            content = b"Name: pysvzerod\nVersion: 3.1\n"
            member = tarfile.TarInfo("pysvzerod-3.2/PKG-INFO")
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
        result = self.cli("verify-dist", "--dist", str(dist), "--version", "3.2")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(source.name, result.stderr)

    def test_rejects_duplicate_metadata_fields(self):
        dist, wheel, _ = self.distributions()
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("pysvzerod-3.2.dist-info/METADATA", "Name: pysvzerod\nVersion: 3.2\nVersion: 3.2\n")
        self.assertNotEqual(self.cli("verify-dist", "--dist", str(dist), "--version", "3.2").returncode, 0)

    def test_index_accepts_complete_matching_files(self):
        helper = self.helper()
        dist, wheel, source = self.distributions()
        response = io.BytesIO(json.dumps(self.remote_files([wheel, source])).encode())
        with patch.object(helper.urllib.request, "urlopen", return_value=response) as request:
            helper.verify_index("testpypi", dist, "3.2", attempts=1, delay=0)
        self.assertEqual(request.call_args.args[0], "https://test.pypi.org/pypi/pysvzerod/3.2/json")

    def test_index_rejects_matching_filename_with_different_hash(self):
        helper = self.helper()
        dist, wheel, source = self.distributions()
        data = self.remote_files([wheel, source])
        data["urls"][0]["digests"]["sha256"] = "0" * 64
        for method in (helper.check_index, helper.verify_index):
            with self.subTest(method=method.__name__):
                with patch.object(helper.urllib.request, "urlopen", return_value=io.BytesIO(json.dumps(data).encode())):
                    with self.assertRaisesRegex(ValueError, "SHA256"):
                        method("pypi", dist, "3.2")

    def test_index_rejects_unexpected_remote_files(self):
        helper = self.helper()
        dist, wheel, source = self.distributions()
        data = self.remote_files([wheel, source])
        data["urls"].append({"filename": "unexpected.whl", "digests": {"sha256": "0" * 64}})
        for method in (helper.check_index, helper.verify_index):
            with self.subTest(method=method.__name__):
                with patch.object(helper.urllib.request, "urlopen", return_value=io.BytesIO(json.dumps(data).encode())):
                    with self.assertRaisesRegex(ValueError, "Unexpected"):
                        method("pypi", dist, "3.2")

    def test_check_index_accepts_partial_matching_upload(self):
        helper = self.helper()
        dist, wheel, _ = self.distributions()
        with patch.object(helper.urllib.request, "urlopen", return_value=io.BytesIO(json.dumps(self.remote_files([wheel])).encode())):
            helper.check_index("pypi", dist, "3.2")

    def test_check_index_accepts_absent_version(self):
        helper = self.helper()
        dist, _, _ = self.distributions()
        with patch.object(helper.urllib.request, "urlopen", side_effect=urllib.error.HTTPError("url", 404, "missing", {}, None)):
            helper.check_index("pypi", dist, "3.2")

    def test_verify_index_retries_missing_files_and_network_errors(self):
        helper = self.helper()
        dist, wheel, source = self.distributions()
        responses = [
            urllib.error.HTTPError("url", 404, "missing", {}, None),
            urllib.error.URLError("connection interrupted"),
            io.BytesIO(json.dumps(self.remote_files([wheel])).encode()),
            io.BytesIO(json.dumps(self.remote_files([wheel, source])).encode()),
        ]
        with patch.object(helper.urllib.request, "urlopen", side_effect=responses) as request:
            with patch.object(helper.time, "sleep") as sleep:
                helper.verify_index("pypi", dist, "3.2", attempts=4, delay=5)
        self.assertEqual(request.call_count, 4)
        self.assertEqual(sleep.call_count, 3)

    def test_verify_index_exhausts_retry_limit(self):
        helper = self.helper()
        dist, wheel, _ = self.distributions()
        data = json.dumps(self.remote_files([wheel])).encode()
        with patch.object(helper.urllib.request, "urlopen", side_effect=lambda *args, **kwargs: io.BytesIO(data)) as request:
            with patch.object(helper.time, "sleep"):
                with self.assertRaisesRegex(ValueError, "Missing"):
                    helper.verify_index("pypi", dist, "3.2", attempts=2, delay=0)
        self.assertEqual(request.call_count, 2)


if __name__ == "__main__":
    unittest.main()
