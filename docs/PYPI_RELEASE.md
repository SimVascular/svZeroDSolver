# Publishing pysvzerod to PyPI

The `.github/workflows/pypi.yml` workflow builds a source distribution and
CPython wheels for Linux x86_64, Windows x86_64, and macOS Intel and Apple
Silicon, and tests each wheel. Pull requests that change packaging files and
manual workflow runs build and test the artifacts without publishing them. A
`v<version>` tag publishes to PyPI after all build jobs pass.

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

Each workflow run keeps its wheels as artifacts named after the runner, such as
`wheels-macos-15` or `wheels-ubuntu-latest`. To try them on another machine,
download one and install it in a fresh environment outside the repository:

```sh
gh run download -R SimVascular/svZeroDSolver --name wheels-macos-15 --dir pysvz-wheels
python -m pip install --find-links pysvz-wheels --only-binary pysvzerod pysvzerod
```

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
   source distribution and wheels, and uploads them to PyPI through the
   trusted publisher once the `pypi` deployment is approved.
3. Install the released package in a fresh environment and check its version:

   ```sh
   python -m pip install pysvzerod==3.1
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

A version number can be uploaded only once. If a release turns out broken, yank
it on PyPI, which makes pip skip it, and publish a fixed version such as
`3.1.1`. For releases with large packaging changes, publish a pre-release
first, for example version `3.2rc1` with tag `v3.2rc1`. Once a final release
exists, pip skips pre-releases unless asked for them.

The wheels use pinned CMake dependencies. Update those pins intentionally when
changing the supported Python versions. The source distribution includes the
CMake project and C++ sources so users without a matching wheel can compile it.
