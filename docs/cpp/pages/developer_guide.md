# Developer Guide

[TOC]

## Code Style

We follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

## Formatting

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to automatically 
format our code accoring to the [Google Style](https://google.github.io/styleguide/cppguide.html).
This increases readability and maintainability of the code while enabling you
to focus on coding. Before committing any changes, you can run the following
comment to format all your files:

```bash
cd src
find **/*.hpp **/*.cpp | xargs clang-format -i --style=Google
```

You can also just check **if** a file would be formatted without actually formatting
it with:

```bash
cd src
find **/*.hpp **/*.cpp | xargs clang-format --dry-run --style=Google --Werror
```

The latter check is also performed in the GitHub CI/CD (a.k.a. Actions) and
indicates on merge requests when the code doesnt't yet meet all style
requirements.

## Documentation

We use [Doxygen](https://doxygen.nl) to automatically build an html documentation
from source code. Please have at Doxygen's [Documentation Guide](https://www.doxygen.nl/manual/docblocks.html)
for an introduction into the syntax of writing documentation in C++. For more
inspiration, you can look at the existing source files and how they use
documentation. In the following you can find a short recap of the most important
commands:

* **Latex equations**: For inline equations use `\f$a+b=c\f$` and for block equations use:
```
\f[
a+b=c
\f]
```
* **Citation**: If you want to cite a piece literature in your documentation, add
    a respective BibTeX citation to `docs/cpp/references.bib` and use `\cite name_of_citation` to
    cite the document.
* **Drawing circuits**: As the elements of the svZeroDSolver are often represented
    in the form of electrical circuits, we make use of [CircuiTikZ](https://ctan.org/pkg/circuitikz?lang=en)
    to draw circuits in the documentation. The start a CircuitTikZ drawing use
    the following command:
```
\f[
\begin{circuitikz}
...
\end{circuitikz}
\f]
```

The documentation is automatically build in the GitHub CI/CD and published
on GitHub pages. If you want to build the documentation locally, you can use:

```
doxygen docs/cpp/Doxyfile
```

If you do not have Doxygen install you can do that with `brew install doxygen`
on macOS or with `sudo apt-get install doxygen` on Linux.

## Profiling

Profiling helps to easily identify bottlenecks in the code. A profiling report
lists the executation time spend on certain parts of the code. If you experience
performance issue with svZeroDSolver, you can follow this guide
to create a profiling report.

### macOS

Install Google Performance Tools:

```bash
brew install gperftools
```

Go to the folder of the build you want to profile (i.e. `cd Debug` or `cd Release`)
and run:

```bash
sudo DYLD_INSERT_LIBRARIES=/usr/local/lib/libprofiler.dylib CPUPROFILE=main.prof CPUPROFILE_FREQUENCY=100000 ./svzerodsolver path/to/input_file.json path/to/output_file.csv
pprof --pdf ./svzerodsolver main.prof > profiling_report.pdf
```

This will generate a profiling report called `profiling_report.pdf` in the respective
folder. Make sure to build svZeroDSolver before starting profiling.

### Linux

Install Google Performance Tools:

```bash
sudo apt install google-perftools libgoogle-perftools-dev
```

Go to the folder of the build you want to profile (i.e. `cd Debug` or `cd Release`)
and run:

```bash
LD_PRELOAD=/usr/local/lib/libprofiler.so CPUPROFILE=main.prof CPUPROFILE_FREQUENCY=100000 ./svzerodsolver ../solver_0d.in ../solver_output.csv
google-pprof --pdf ./svzerodsolver main.prof > profiling_report.pdf
```

If the first line of the above command doesn't work, you can try
exchanging the `LD_PRELOAD` path with `/usr/lib/x86_64-linux-gnu/libprofiler.so`.

This will generate a profiling report called `profiling_report.pdf` in the respective
folder. Make sure to build svZeroDSolver before starting profiling.
