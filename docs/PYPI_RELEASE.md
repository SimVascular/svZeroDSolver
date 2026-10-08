# Publishing pysvzerod to PyPI

The `.github/workflows/pypi.yml` workflow builds a source distribution and
CPython wheels for Linux x86_64 and aarch64, Windows x86_64, and macOS Intel and
Apple Silicon, and tests each wheel. Linux wheels use native runners, including
`ubuntu-24.04-arm` for aarch64. Pull requests that change packaging files and
manual workflow runs build, test, and validate the artifacts without publishing
them. A `v<version>` tag publishes to PyPI after validation passes.

Each wheel runs `tests/test_solver.py`,
`tests/test_calibrator.py::test_steady_flow_calibration`, and the `svzerodsolver`
console command. The full I/O and calibration suites run in normal CI. CMake
package builds are limited to two parallel jobs.

## First-time setup

1. In the `SimVascular/svZeroDSolver` repository, create a GitHub environment
   named `pypi`. Limit its deployments to `v*` tags and add a required
   reviewer, because an upload cannot be undone.
2. On PyPI, create `pysvzerod` in the SimVascular organization (Your
   organizations, Manage, Projects). In the project's Publishing settings, add
   a GitHub trusted publisher with owner `SimVascular`, repository
   `svZeroDSolver`, workflow `pypi.yml`, and environment `pypi`.

A pending publisher in a personal account also works, but its first upload
creates the project under that account, and an organization owner must then
transfer it. Pending publishers do not reserve the project name until the first
upload.

## Testing wheels before a release

On pull requests, manual runs, and tags, a separate `validate` job runs
`twine check --strict` on the source distribution and all wheels after their
builds and tests pass. It saves the checked files as the `python-distributions`
artifact. To try a wheel on another machine, download that artifact and install
it in a fresh environment outside the repository:

```sh
gh run download -R SimVascular/svZeroDSolver --name python-distributions --dir pysvz-wheels
python -m pip install numpy pandas
python -m pip install --no-index --find-links pysvz-wheels --only-binary=:all: --no-deps pysvzerod==3.1
```

Replace `3.1` with the version in the artifact being tested.
Without a run ID, `gh run download` takes the latest artifact with that name;
pass a run ID to pick a specific run. Then run a model and the `svzerodsolver`
command.

## Each release

Repository releases and the package share one version number. Every new `v*`
tag starts the workflow and must equal `v` plus the `pyproject.toml` version,
including a tag created for a GitHub release.

1. Set the version in `pyproject.toml` in a pull request and merge it. Run the
   workflow manually on that commit and test its wheels as described above.
2. Create a GitHub release on that commit with a tag matching the version, for
   example `v3.1`. The tag build verifies the version, builds and tests the
   source distribution and wheels, then validates their metadata. After
   validation, the publish job downloads and uploads that same
   `python-distributions` artifact through the trusted publisher once the
   `pypi` deployment is approved.
3. Install the released package in a fresh environment and check its version:

   ```sh
   python -m pip install pysvzerod==3.1
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

Uploads are serialized with `cancel-in-progress: false`, so a newer run does
not cancel an upload in progress. If an upload fails partway through, inspect
the files already on PyPI before retrying; existing filenames cannot be
overwritten.

A version number can be uploaded only once. If a release turns out broken, yank
it on PyPI, which makes pip skip it, and publish a fixed version such as
`3.1.1`. For releases with large packaging changes, publish a pre-release
first, for example version `3.2rc1` with tag `v3.2rc1`. Once a final release
exists, pip skips pre-releases unless asked for them.

The wheels use pinned CMake dependencies. Update those pins intentionally when
changing the supported Python versions. The source distribution includes the
CMake project and C++ sources so users without a matching wheel can compile it.
