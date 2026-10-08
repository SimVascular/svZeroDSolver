# Publishing pysvzerod

## One-time setup

1. Create GitHub environments named `pypi` and `testpypi`. Restrict `pypi` to
   the default branch (`master`); allow testing branches in `testpypi`.
2. Configure a trusted publisher for `pysvzerod` on each package index:
   owner `SimVascular`, repository `svZeroDSolver`, workflow `pypi.yml`, and
   the matching environment (`pypi` or `testpypi`).
3. Allow the workflow to commit version updates to the default branch, create
   tags, and create GitHub releases. Branch protection must permit these updates.

## Run a release

1. Open **Actions → Build and publish PyPI package → Run workflow**.
2. Select the source branch. For production, merge your changes into `master`
   first and select that branch.
3. Choose a `destination`:

   | Destination | Result |
   | --- | --- |
   | `testpypi` (default) | Upload the next version to TestPyPI for testing. |
   | `pypi` | Upload the next version to PyPI and create a GitHub release. |

4. Choose `version_bump`: **minor** (`3.1 → 3.2`) or **major** (`3.1 → 4.0`).
   Click **Run workflow** and approve the deployment if reviewers are configured.

The workflow builds, tests, and validates the packages before uploading. A
successful PyPI release updates `pyproject.toml`, creates a `vX.Y` tag, and
creates a GitHub release with notes and the packages attached. Do not update
versions or create release tags manually.

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
