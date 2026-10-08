# Publishing pysvzerod to PyPI

The `.github/workflows/pypi.yml` workflow builds a source distribution and
CPython wheels for Linux x86_64 and aarch64, Windows x86_64, and macOS Intel and
Apple Silicon, and tests each wheel. Linux wheels use native runners, including
`ubuntu-24.04-arm` for aarch64. Start releases manually from GitHub Actions;
tags do not trigger publication. Pull requests that change packaging files
build, test, and validate the current package version without publishing.

Each wheel runs `tests/test_solver.py`,
`tests/test_calibrator.py::test_steady_flow_calibration`, and the `svzerodsolver`
console command. The full I/O and calibration suites run in normal CI. CMake
package builds are limited to two parallel jobs.

## First-time setup

1. In the `SimVascular/svZeroDSolver` repository, create GitHub environments
   named `pypi` and `testpypi`. Restrict `pypi` deployments to the default
   branch, currently `master`; allow the feature branches you want to test in
   `testpypi`. Add required reviewers if release approval is desired.
2. On PyPI, create `pysvzerod` in the SimVascular organization (Your
   organizations, Manage, Projects). In the project's Publishing settings, add
   a GitHub trusted publisher with owner `SimVascular`, repository
   `svZeroDSolver`, workflow `pypi.yml`, and environment `pypi`.
3. Create a separate `pysvzerod` project and trusted publisher on TestPyPI,
   using the same owner, repository, and workflow with environment `testpypi`.
   Both destinations authenticate through OIDC trusted publishing.
4. Enable repository Actions write permissions and permit the workflow to
   commit the version update directly to the default branch, create tags, and
   create GitHub Releases. If branch protection blocks the version commit,
   configure designated release permissions or reconcile the release manually
   as described below.

A pending publisher in a personal account also works, but its first upload
creates the project under that account, and an organization owner must then
transfer it. Pending publishers do not reserve the project name until the first
upload.

## Choose a run

In Actions, open **Build and publish PyPI package**, select **Run workflow**,
and choose the source branch and `destination`:

| Destination | Behavior |
| --- | --- |
| `build-only` (default) | Build and validate the current `pyproject.toml` version; upload nothing to a package index. |
| `testpypi` | Publish the selected next version to TestPyPI; feature branches are allowed. |
| `pypi` | Publish the selected next version to PyPI from the current default branch and record the release. |

For publication, choose `version_bump`: `minor` (default) changes the current
`3.1` to `3.2`; `major` changes it to `4.0`. Build-only and pull-request runs
keep `3.1` while that is the repository version. Do not edit the version before
a release: the workflow calculates it. Patch and prerelease versions are not
supported by this workflow. Publication requires explicitly selecting
`testpypi` or `pypi`.

The workflow pins the source commit and injects the same selected version into
every package build. After builds and tests pass, a separate `validate` job
checks every distribution's package name and version and runs
`twine check --strict` on the source distribution and all wheels. It saves the
checked files as `python-distributions`; publishing downloads that same
artifact.

## Test the artifacts

For a build-only run, download the validated artifact and install its wheel in
a fresh environment outside the repository:

```sh
gh run download RUN_ID -R SimVascular/svZeroDSolver --name python-distributions --dir pysvz-wheels
python -m pip install numpy pandas
python -m pip install --no-index --find-links pysvz-wheels --only-binary=:all: --no-deps pysvzerod==3.1
```

Replace `RUN_ID` with the workflow run ID and `3.1` with its artifact version.
Then run a model and the `svzerodsolver` command.

After publishing a candidate to TestPyPI, test its exact version in another
fresh environment. Install dependencies from PyPI first, then install only the
package from TestPyPI:

```sh
python -m pip install numpy pandas
python -m pip install --index-url https://test.pypi.org/simple/ --only-binary=:all: --no-deps pysvzerod==3.2
```

Replace `3.2` with the candidate version. TestPyPI publication verifies the
uploaded files and hashes, but does not commit the version or create a stable
GitHub release. A new build cannot replace an existing TestPyPI version with
different files. Use build-only runs for further development until an unused
major or minor candidate is available.

## Each release

1. Merge the intended source changes into the default branch. Run with
   `destination: build-only` and test the downloaded artifacts. TestPyPI is
   available for checking a publishing candidate before production.
2. Run the workflow on the default branch with `destination: pypi` and the
   intended `version_bump`. Approve the `pypi` deployment if a reviewer is
   configured. After uploading, the workflow verifies the files and hashes on
   PyPI, commits the new `pyproject.toml` version as
   `Zachary Sexton <zsexton@stanford.edu>`, tags it as `vX.Y`, and creates a
   GitHub Release with generated notes and the validated distribution files.
3. Install the released package in a fresh environment and check its version:

   ```sh
   python -m pip install pysvzerod==3.2
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

   Replace `3.2` with the released version.

## Retry and recovery

Production runs are serialized from preparation through release recording;
TestPyPI runs are serialized separately. Both use `cancel-in-progress: false`.
Package versions and uploaded filenames are immutable. If an upload stops
partway through, rerun the failed jobs to reuse the same validated artifact.
Do not rerun all jobs or start a fresh build to retry that upload: rebuilding
can produce different bytes. Before uploading, the workflow rejects existing
files that conflict with the artifact. It skips matching existing files, then
verifies the complete set of filenames and hashes on the selected index.

The version commit and tag are pushed together atomically. If the default
branch moves before that update, automated recording stops instead of
force-pushing. An upload may already have succeeded.
Use the run's exact source commit, selected version, and validated artifact to
reconcile the version commit, tag, and GitHub Release before starting another
production release. The same reconciliation may be needed if repository
permissions block recording after upload. If a published release is broken,
yank it on PyPI and publish a corrected release with a new minor or major
version.

If the version commit and tag were recorded but creating the GitHub Release
was interrupted, rerun the failed job. It recognizes the exact recorded commit
and resumes without another version increment. Existing release assets must
match the validated files; completed releases are left as they are.

The wheels use pinned CMake dependencies. Update those pins intentionally when
changing the supported Python versions. The source distribution includes the
CMake project and C++ sources so users without a matching wheel can compile it.
