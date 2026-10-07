# Publishing pysvzerod to PyPI

The `.github/workflows/pypi.yml` workflow builds a source distribution and
CPython wheels for Linux x86_64, Windows x86_64, and macOS Intel and Apple
Silicon. A manual workflow run builds and checks the artifacts without
publishing. A `v<version>` tag publishes them after all build jobs pass.

## First-time setup

1. In the GitHub repository that will own releases, create an environment named
   `pypi`. Configure required reviewers if releases need manual approval.
2. In PyPI account settings, create a pending trusted publisher for the
   `pysvzerod` project. Enter the exact GitHub owner and repository that will
   run the workflow, the workflow filename `pypi.yml`, and environment `pypi`.
   A pending publisher does not reserve the name; publish promptly after setup.
3. Merge the packaging changes into that repository's default branch and run
   the workflow manually once. Inspect all wheel and source build results.

## Each release

1. Set the version in `pyproject.toml`. Run the project tests and the
   PyPI workflow manually from the release commit.
2. After the changes are merged, create and push a tag matching that version,
   for example `v2.0`. The tag build verifies the version, builds the source
   distribution and wheels, checks metadata, and uploads to PyPI through the
   trusted publisher.
3. Install the released package in a fresh environment and check its version:

   ```sh
   python -m pip install pysvzerod==2.0
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

The wheels use pinned CMake dependencies. Update those pins intentionally when
changing the supported Python versions. The source distribution includes the
CMake project and C++ sources so users without a matching wheel can compile it.
