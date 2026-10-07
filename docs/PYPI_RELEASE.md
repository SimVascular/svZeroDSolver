# Publishing pysvzerod to PyPI

The `.github/workflows/pypi.yml` workflow builds a source distribution and
CPython wheels for Linux x86_64, Windows x86_64, and macOS Intel and Apple
Silicon. A manual workflow run builds and checks the artifacts. Select the
`publish_testpypi` input to upload them to TestPyPI; its default is off.
A `v<version>` tag publishes to production PyPI after all build jobs pass.

## First-time setup

1. In the GitHub repository that will own releases, create an environment
   named `testpypi`. Configure required reviewers if uploads need manual
   approval.
2. [Create a TestPyPI account](https://test.pypi.org/account/register/). It is
   separate from a production PyPI account. In its account Publishing settings,
   create a pending trusted publisher for `pysvzerod` with the exact GitHub
   owner and repository that will run this workflow, filename `pypi.yml`, and
   environment `testpypi`.
3. Merge the packaging changes into that repository's default branch. Run
   the workflow manually with `publish_testpypi` selected. The workflow builds
   and tests the packages before uploading them to TestPyPI. Inspect the
   [TestPyPI project page](https://test.pypi.org/project/pysvzerod/) and install
   the exact version in a fresh environment:

   ```sh
   python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pysvzerod==2.0.1
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

Before publishing to production PyPI, create a separate `pypi` GitHub
environment and a pending publisher in the production PyPI account with the
same owner, repository, and workflow filename, but environment `pypi`.
Pending publishers do not reserve the project name until the first upload.

## Local TestPyPI upload before merging

For a quick test from this branch, create a TestPyPI account and an API token
there. Build and inspect the source distribution from a clean checkout:

```sh
uv build --sdist
uvx twine check dist/*
uvx twine upload --repository testpypi dist/*
```

Twine prompts for a username and password. Enter `__token__` as the username
and the TestPyPI API token as the password. Upload only the source
distribution: a local wheel is tagged for this machine and is not repaired, and
TestPyPI rejects Linux `linux_x86_64` wheels. The GitHub workflow builds the
full wheel matrix. Each filename can be uploaded only once, even after
deletion, so set a dev version such as `2.0.2.dev1` in `pyproject.toml` for
repeated local tests.

## Each release

1. Set the version in `pyproject.toml`. Run the project tests and the
   PyPI workflow manually from the release commit.
2. After the changes are merged, create and push a tag matching that version,
   for example `v2.0.1`. The tag build verifies the version, builds the source
   distribution and wheels, checks metadata, and uploads to PyPI through the
   trusted publisher.
3. Install the released package in a fresh environment and check its version:

   ```sh
   python -m pip install pysvzerod==2.0.1
   python -c "import importlib.metadata; print(importlib.metadata.version('pysvzerod'))"
   ```

The wheels use pinned CMake dependencies. Update those pins intentionally when
changing the supported Python versions. The source distribution includes the
CMake project and C++ sources so users without a matching wheel can compile it.
