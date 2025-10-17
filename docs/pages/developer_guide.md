@page developer_guide Developer Guide

[TOC]

# Architecture

svZeroDSolver is written in a highly modular manner to enable reuse of code
for many different applications. It is divided into a header based library
in the `src` directory and a collection of different applications in the
`applications` folder. Each application is written for a different use-case
of svZeroDSolver, namely:

* svZerodSolver in `svzerodsolver.cpp`
* Python API in `pysvzerod.cpp`
* svZeroDCalibrator in `svzerodcalibrator.cpp`
* svZeroDVisualization for visualizing 0D models and results
* svZeroDGUI for creating new 0D models grahically.

[Architecture for svZeroDVisualization](@ref visualization).

[Architecture for svZeroDGUI](@ref GUI).


# Build in debug mode

For debug purposes it is recommended to build svZeroDSolver in Debug mode.

```bash
mkdir Debug
cd Debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build .
```

# Install with pip

Execute this command in the root folder to install the current source:
```bash
pip install -e ".[dev]"
```
This is useful when continuously running the integration tests during development.

# Contributing to svZeroDSolver

**NOTE: To contribute new developments to the main branch of svZeroDSolver, developers must first open an issue on the svZeroDSolver Github repository to describe the planned changes.** 

* The changes should be implemented in a feature branch of the developer's fork of svZeroDSolver. 
* Once the changes are implemented, the developer should make sure the build, documentation, and code format tests are passing on the user's feature branch. 

# Testing code changes
* Tests are automatically run as Github workflows when pushing changes to the developer's remote branch on Github. These tests are defined in the `.github/workflows` directory.
* It may be helpful and more efficient to run the build tests locally with each incremental change the developer makes. This can be done by following these steps:
  1. If you do not have `pytest` installed in your virtual environment, install it with `pip install pytest`
  2. We first need to build the svZeroDSolver project. This can be done by running `pip install -e .` from the svZeroDSolver directory
  3. From the main svZeroDSolver directory, `cd tests`.
  4. Now that you are in the `tests` directory, run `pytest`. The terminal output will show which tests are passing.
    * To run a specific subset of tests matching a certain pattern, e.g., only the calibrator, you can use `pytest -k calibrator`. 
    * To see all simulation output, run `pytest -sv`.
    * If you want to run the tests using the binaries in the `Release` directory instead of the python module, you can call `pytest --coverage`
* The tests for the C++ interface require the `CMake` install and can be run by building the tests in `svZeroDSolver/tests/test_interface`.
* Code formatting can be performed following the instructions in the [Formatting](#formatting) section. 
* The documentation can be built following the instructions in the [Documentation](#documentation) section. 
* Once all the tests are passing, the developer should open a pull request from the feature branch and link the relevant issue.

# Adding new blocks

The modular architecture of svZeroDSolver relies on "blocks", such as blood vessels, junctions, valves, boundary conditions, etc. These blocks are assembled in a manner specified by the `.json` configuration file, which dictates the assembled governing equations for the model. We are always interested in adding new blocks to expand the funcitonality of svZeroDSolver.

Detailed steps required to implement a new block in svZeroDSolver are available [here](@ref add_block).

Steps required to visualize a new block with svZeroDSolver Visualization application are available [here](@ref visualization).

# Code Style

We follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

## Formatting {#formatting}

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to automatically 
format our code accoring to the [Google Style](https://google.github.io/styleguide/cppguide.html), 
as specified in the `.clang-format` file. This increases readability and maintainability of the code 
while enabling you to focus on coding.

There are tools for your favorite IDE to automatically format your code. Examples are:
- [Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)
- [vim](https://github.com/rhysd/vim-clang-format)
- [and many more](https://clang.llvm.org/docs/ClangFormat.html)

Before formatting the code with clang-format, ensure that you have latest clang-format version by running `clang-format --version`. To upgrade `clang-format` use the following commands
  * on MacOS: `brew upgrade clang-format`
  * on Ubuntu: `sudo apt upgrade clang-format`

Before committing any changes, run the following
command **from your build folder** to format all your files:

```bash
make codeformat
```

You can also just check **if** a file would be formatted without actually formatting
it with:

```bash
make codecheck
```

If the above commands do not work on your platform (it does not work on Sherlock at Stanford)
you can run the following command **from the svZeroDSolver folder** to format all your files:

```bash
find src/**/*.h src/**/*.cpp | xargs clang-format -style=Google -i
```

The latter check is also performed in the GitHub CI/CD (a.k.a. Actions) and
indicates on merge requests when the code doesnt yet meet all style
requirements.

On Sherlock at Stanford, clang-format is included in the `llvm` module.

# Documentation {#documentation}

We use [Doxygen](https://doxygen.nl) to automatically build an html documentation
from source code. Please have at Doxygen's [Documentation Guide](https://www.doxygen.nl/manual/docblocks.html)
for an introduction into the syntax of writing documentation in C++. For more
inspiration, you can look at the existing source files and how they use
documentation. 

**NOTE: Undocumented code will fail the automated code checks on Github 
and cannot be merged.**

In the following you can find a short recap of the most important
commands:

## Latex equations
For inline equations use `\f$a+b=c\f$` and for block equations use:
```
\f[
a+b=c
\f]
```

## Citations
If you want to cite a piece literature in your documentation, add
a respective BibTeX citation to `docs/references.bib` and use `\cite name_of_citation` to
cite the document.

## Drawing circuits
As the elements of the svZeroDSolver are often represented
in the form of electrical circuits, we use [CircuiTikZ](https://ctan.org/pkg/circuitikz?lang=en)
to draw circuits in the documentation (see blocks in Block for examples). 
To start a CircuitTikZ drawing use the following command:
```
\f[
\begin{circuitikz}
...
\end{circuitikz}
\f]
```

## Build
The documentation is automatically built in the GitHub CI/CD and published
on GitHub pages. If you want to build the documentation locally, you can use:

```
doxygen docs/Doxyfile
```
You can then view the documentation locally in your browser by opening `docs/build/html/index.html`.

If you do not have Doxygen install you can do that with `brew install doxygen`
on macOS or with `sudo apt-get install doxygen` on Linux.

# Profiling

Profiling helps to easily identify bottlenecks in the code. A profiling report
lists the executation time spend on certain parts of the code. If you experience
performance issue with svZeroDSolver, you can follow this guide
to create a profiling report. The generation of profiling reports requires
[Docker](https://docs.docker.com/get-docker/) to be installed on your machine.

```docker
docker build -t profile_svzerodsolver -f container/profiling/Dockerfile .
docker run -it -v $(PWD):/opt/data --rm profile_svzerodsolver path/to/simulation_config.json
```

This will generate a file called `profiling_report.pdf` in your current working directory.
