# Publishing svzerod

## Run a release

1. Open **Actions → Build and publish PyPI package → Run workflow**.
2. Select `master` as the source branch, with your changes already merged.
3. Choose a `destination`:

   | Destination | Result |
   | --- | --- |
   | `testpypi` (default) | Upload the next version to TestPyPI for testing. |
   | `pypi` | Upload the next version to PyPI and create a GitHub release. |

4. Choose `version_bump`: **minor** (`3.1 → 3.2`) or **major** (`3.1 → 4.0`).
5. Click **Run workflow**. For PyPI, have a configured reviewer approve the
   waiting deployment.

Use `testpypi` to check a release before selecting `pypi` for production.

The workflow builds, tests, and validates the packages before uploading. A
successful PyPI release updates `pyproject.toml`, creates a `vX.Y` tag, and
creates a GitHub release with notes and the packages attached. Do not update
versions or create release tags manually.

The uploaded package reports its version through `svzerod.__version__`.

TestPyPI uploads leave the repository version unchanged and do not create a
GitHub release. Pull requests and tags never publish packages automatically.

## If a release fails

- If the default branch changed before uploading, start a new run from its
  latest commit.
- Use **Re-run failed jobs** to retry an interrupted upload or GitHub release
  with the same packages. Published versions cannot be overwritten; avoid
  rebuilding an already uploaded version.
- If PyPI upload succeeds but recording the version is blocked by branch changes
  or permissions, reconcile the version commit, tag, and GitHub release with
  that run's source and packages before starting another production release.
