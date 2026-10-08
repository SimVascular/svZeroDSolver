"""Exercise release recording against isolated local Git repositories."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "record_release.sh"
VERSION_HELPER = """\
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("command", choices=["set"])
parser.add_argument("--expected-current", required=True)
parser.add_argument("--version", required=True)
args = parser.parse_args()
path = Path("pyproject.toml")
text = path.read_text()
old = f'version = "{args.expected_current}"'
assert text.count(old) == 1
path.write_text(text.replace(old, f'version = "{args.version}"'))
"""


class RecordReleaseTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.remote = self.root / "remote.git"
        self.seed = self.root / "seed"
        self.environment = os.environ.copy()
        for key in list(self.environment):
            if key.startswith("GIT_"):
                del self.environment[key]
        self.environment.update(
            GIT_CONFIG_GLOBAL=os.devnull,
            GIT_CONFIG_NOSYSTEM="1",
            GIT_AUTHOR_NAME="Fixture Author",
            GIT_AUTHOR_EMAIL="fixture@example.com",
            GIT_COMMITTER_NAME="Fixture Author",
            GIT_COMMITTER_EMAIL="fixture@example.com",
        )
        self.git(self.root, "init", "--bare", "--initial-branch=master", str(self.remote))
        self.git(self.root, "clone", str(self.remote), str(self.seed))
        (self.seed / ".github/scripts").mkdir(parents=True)
        (self.seed / ".github/scripts/release_version.py").write_text(VERSION_HELPER)
        (self.seed / "pyproject.toml").write_text(
            '[project]\nname = "pysvzerod"\nversion = "3.1"\n'
        )
        (self.seed / "README.md").write_text("A solver package.\n")
        self.git(self.seed, "add", ".")
        self.git(self.seed, "commit", "-m", "Initial fixture")
        self.source = self.git(self.seed, "rev-parse", "HEAD")
        self.git(self.seed, "push", "origin", "master")
        self.checkout_number = 0

    def git(self, directory, *arguments):
        result = subprocess.run(
            ["git", *arguments],
            cwd=directory,
            env=self.environment,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    def checkout(self):
        self.checkout_number += 1
        directory = self.root / f"checkout-{self.checkout_number}"
        self.git(self.root, "clone", str(self.remote), str(directory))
        self.git(directory, "checkout", "--detach", self.source)
        return directory

    def record(self, directory):
        output = directory / "release-output"
        environment = self.environment | {
            "DEFAULT_BRANCH": "master",
            "SOURCE_SHA": self.source,
            "CURRENT_VERSION": "3.1",
            "RELEASE_VERSION": "3.2",
            "GITHUB_OUTPUT": str(output),
        }
        result = subprocess.run(
            ["bash", str(SCRIPT)],
            cwd=directory,
            env=environment,
            capture_output=True,
            text=True,
        )
        return result, output.read_text() if output.exists() else ""

    def refs(self):
        return self.git(self.remote, "show-ref")

    def release_commit(self):
        return self.git(self.remote, "rev-parse", "refs/tags/v3.2^{commit}")

    def test_records_version_and_tag_in_one_commit(self):
        checkout = self.checkout()
        (checkout / "dist").mkdir()
        (checkout / "dist/package.whl").write_text("Untracked distribution")
        result, output = self.record(checkout)
        self.assertEqual(result.returncode, 0, result.stderr)
        commit = self.release_commit()
        self.assertEqual(output, f"commit_sha={commit}\n")
        self.assertEqual(self.git(self.remote, "rev-parse", "master"), commit)
        self.assertEqual(self.git(self.remote, "rev-parse", f"{commit}^"), self.source)
        self.assertEqual(
            self.git(self.remote, "diff", "--name-only", self.source, commit),
            "pyproject.toml",
        )
        self.assertEqual(
            self.git(self.remote, "show", "-s", "--format=%an <%ae>|%cn <%ce>|%s", commit),
            "github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>|"
            "github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>|"
            "Release pysvzerod 3.2",
        )
        self.assertIn('version = "3.2"', self.git(self.remote, "show", f"{commit}:pyproject.toml"))

    def test_matching_tag_retry_reuses_the_recorded_commit(self):
        first, _ = self.record(self.checkout())
        self.assertEqual(first.returncode, 0, first.stderr)
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output, f"commit_sha={self.release_commit()}\n")
        self.assertEqual(self.refs(), before)

    def test_matching_retry_accepts_later_default_branch_commits(self):
        first, _ = self.record(self.checkout())
        self.assertEqual(first.returncode, 0, first.stderr)
        self.git(self.seed, "pull", "--ff-only", "origin", "master")
        (self.seed / "README.md").write_text("A later change.\n")
        self.git(self.seed, "commit", "-am", "Update fixture documentation")
        self.git(self.seed, "push", "origin", "master")
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output, f"commit_sha={self.release_commit()}\n")
        self.assertEqual(self.refs(), before)

    def test_advanced_branch_without_a_tag_is_rejected(self):
        checkout = self.checkout()
        (self.seed / "README.md").write_text("A later change.\n")
        self.git(self.seed, "commit", "-am", "Advance fixture branch")
        self.git(self.seed, "push", "origin", "master")
        before = self.refs()
        result, output = self.record(checkout)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("advanced", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)

    def test_rejected_branch_update_does_not_publish_the_tag(self):
        hook = self.remote / "hooks/update"
        hook.write_text('#!/bin/sh\n[ "$1" != "refs/heads/master" ]\n')
        hook.chmod(0o755)
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("could not be recorded", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)

    def test_existing_tag_with_wrong_parent_is_rejected(self):
        self.git(self.seed, "tag", "v3.2", self.source)
        self.git(self.seed, "push", "origin", "v3.2")
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)

    def test_existing_tag_with_wrong_tree_is_rejected(self):
        (self.seed / "README.md").write_text("Unrelated tagged change.\n")
        self.git(self.seed, "commit", "-am", "Create an unrelated tagged commit")
        self.git(self.seed, "tag", "v3.2")
        self.git(self.seed, "push", "origin", "master", "v3.2")
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)

    def test_matching_tag_outside_default_branch_is_rejected(self):
        metadata = self.seed / "pyproject.toml"
        metadata.write_text(metadata.read_text().replace('"3.1"', '"3.2"'))
        self.git(self.seed, "commit", "-am", "Create an unmerged release")
        self.git(self.seed, "tag", "v3.2")
        self.git(self.seed, "push", "origin", "v3.2")
        before = self.refs()
        result, output = self.record(self.checkout())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)

    def test_tracked_changes_are_rejected_before_editing_the_version(self):
        for staged in (False, True):
            with self.subTest(staged=staged):
                checkout = self.checkout()
                metadata = (checkout / "pyproject.toml").read_text()
                (checkout / "README.md").write_text("Local edits.\n")
                if staged:
                    self.git(checkout, "add", "README.md")
                before = self.refs()
                result, output = self.record(checkout)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("clean", result.stderr)
                self.assertEqual((checkout / "pyproject.toml").read_text(), metadata)
                self.assertEqual(output, "")
                self.assertEqual(self.refs(), before)

    def test_wrong_source_checkout_is_rejected(self):
        checkout = self.checkout()
        self.git(checkout, "commit", "--allow-empty", "-m", "Move checkout")
        before = self.refs()
        result, output = self.record(checkout)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("SOURCE_SHA", result.stderr)
        self.assertEqual(output, "")
        self.assertEqual(self.refs(), before)


if __name__ == "__main__":
    unittest.main()
