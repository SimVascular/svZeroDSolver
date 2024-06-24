@page developer_guide Developer Guide

[TOC]

# Architecture

svZeroDSolver is written in a highly modular manner to enable reuse of code
for many different applications. It is divided into a header based library
in the `src` directory and a collection of different applications in the
`applications` folder. Each application is written for a different use-case
of svZeroDSolver, namely:

* svZeroDCalibrator in `svzerodcalibrator.cpp`
* svZerodSolver in `svzerodsolver.cpp`
* Python API in `pysvzerod.cpp`


## Build in debug mode

For debug purposes it is recommended to build svZeroDSolver in Debug mode.

```bash
mkdir Debug
cd Debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build .
```

## Install with pip

Execute this command in the root folder to install the current source:
```bash
pip install -e ".[dev]"
```
This is useful when continuously running the integration tests during development.

## Adding new blocks

**1. Add the new block to the following lists/dictionaries:**
- `BlockType` in src/model/BlockType.h
- `block_factory_map` in src/model/Model.cpp. Here the dictionary key should match the string that specifies the type of block read from the JSON config file.
- If required, `BlockClass` if it's a new kind of block that might require special handling (most blocks don't).

**2. Create a new class inherited from `Block` for the new block. Specific things to pay attention to:**
- Define a constructor of the form:
`GenericBlock(int id, Model *model) : Block(id, model, BlockType::block_type, BlockClass::block_class, {{"Param1", InputParameter()}, {"Param2", InputParameter()}, ..., {"ParamN", InputParameter()}}) {}`
- In the constructor, `GenericBlock` is the name of the new class
- `block_type` and `block_class` are the same as what was added in Step 1 above.
- The input parameters of the block are `Param1`, ... , `ParamN`. The properties of each parameter is defined by `InputParameter`, which specifies whether it is optional, an array, a scalar, and its default value. See `InputParameter` in `Parameter.h` for details.
- The names `Param1`, ... , `ParamN` should be the same as the strings read from the JSON config file.
- The class can have an  `enum ParamId` object that relates the parameter indices to their names. This makes it easier to reference the parameters while implementing the governing equations of the block (discussed below). The order of parameters in the `ParamId` object should match the order in the constructor.
- The class should have `setup_dofs` and `update_constant` functions. It should also have a `TripletsContributions num_triplets{?, ?, ?}` object that specifies how many elements the governing equations of the block contribute to the `F`, `E` and `dCdy` matrices respectively. The class may also contain `update_time` and `update_solution` functions. Details are below.

**3. Now implement the governing equations for the block.**
- The local state vector for each block is always arranged as `[P_in, Q_in, P_out, Q_out, InternalVariable1, ..., InternalVariableN]`. Here, `InternalVariable*` refers to any variable in the governing equations that are not in the inlet and outlet flow and pressure.
- The equations should be written in the form `E(t)*ydot + F(t)*y + c(y,t) = 0`. Here, `y` is the local state vector mentioned above, `ydot` is the time-derivative of the state vector, `E` and `F` are matrices and `c` is a vector containing all non-linear and constant terms in the equation. `E` and `F` are size (number_of_equations*size_of_state_vector), while `c` is length number_of_equations.
- All matrix elements that are constant are specified in `update_constant`. Matrix elements that depend only on time (not the solution of the problem itself) are specified in `update_time`. Matrix elements that change with the solution (i.e. depend on the state variables themselves) are specified in `update_solution`. Not all blocks will require the latter two functions.
- The elements of the `E` and `F` matrix are populated using the syntax `system.F.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_ids]) = A`. Here, `current_block_equation_id` goes from 0 to number_of_equations (for the current block) and `current_block_variable_ids` goes from 0 to size_of_state_vector for the current block. The indices correspond to the block's local state vector mentioned above.
- If the block contains non-linear equations, these terms must be specified in `update_solution` as `system.C(global_eqn_ids[current_block_equation_id]) = non_linear_term`.
- For non-linear equations, the derivative of the term specified above with respect to each state variable must also be provided. This goes into a `dC_dy` matrix using the following syntax `system.dC_dy.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_id]) = A`, where `A` is the derivative of the non-linear term in the equation with ID `current_block_equation_id` with respect to the local state variable with ID `current_block_variable_id'. For example, if the non-linear term is in the first equation, `current_block_equation_id = 0`. For the derivative of this term with respect to `P_in`, `current_block_variable_id = 0` and for the derivative of this term with respect to `P_out`, `current_block_variable_id = 2`.

**4. Add the new files (`GenericBlock.h` and `GenericBlock.cpp`) to src/model/CMakeLists.txt**

This looks quite daunting when listed out like this. We need a way to explain this better before putting it into the documentation.


## Code Style

We follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

## Formatting

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to automatically 
format our code accoring to the [Google Style](https://google.github.io/styleguide/cppguide.html), 
as specified in the `.clang-format` file. This increases readability and maintainability of the code 
while enabling you to focus on coding.

There are tools for your favorite IDE to automatically format your code. Examples are:
- [Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)
- [vim](https://github.com/rhysd/vim-clang-format)
- [and many more](https://clang.llvm.org/docs/ClangFormat.html)

Before committing any changes, you can run the following
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

## Documentation

We use [Doxygen](https://doxygen.nl) to automatically build an html documentation
from source code. Please have at Doxygen's [Documentation Guide](https://www.doxygen.nl/manual/docblocks.html)
for an introduction into the syntax of writing documentation in C++. For more
inspiration, you can look at the existing source files and how they use
documentation. 

**NOTE: Undocumented code will fail the automated code checks on Github 
and cannot be merged.**

In the following you can find a short recap of the most important
commands:

### Latex equations
For inline equations use `\f$a+b=c\f$` and for block equations use:
```
\f[
a+b=c
\f]
```

### Citations
If you want to cite a piece literature in your documentation, add
a respective BibTeX citation to `docs/references.bib` and use `\cite name_of_citation` to
cite the document.

### Drawing circuits
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

### Build
The documentation is automatically built in the GitHub CI/CD and published
on GitHub pages. If you want to build the documentation locally, you can use:

```
doxygen docs/Doxyfile
```
You can then view the documentation locally in your browser by opening `docs/build/html/index.html`.

If you do not have Doxygen install you can do that with `brew install doxygen`
on macOS or with `sudo apt-get install doxygen` on Linux.

## Profiling

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
