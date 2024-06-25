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

The modular architecture of svZeroDSolver relies on "blocks", such as blood vessels, junctions, valves, boundary conditions, etc. These blocks are assembled in a manner specified by the `.json` configuration file, which dictates the assembled governing equations for the model. We are always interested in adding new blocks to expand the funcitonality of svZeroDSolver.

To add a new block, developers must first open an issue to describe the planned block on the svZeroDSolver Github repository. The new block should be implemented in a feature branch of the developer's fork of svZeroDSolver. Once the new block is implemented, the developer should open a pull request from the feature branch and link the relevant issue.

Below are details on the steps required to implement a new block in svZeroDSolver.

**1. Add the new block to the following lists/dictionaries:**
* `BlockType` in src/model/BlockType.h
* `block_factory_map` in src/model/Model.cpp
  * *Note:* In `block_factory_map`, the dictionary key should match the string specifying the type of block in the `.json` configuration/input file, and the dictionary value should match the class constructor name for the block.
* If the new block requires special handling that is different from the current blocks (most new blocks do not), add a new category to `BlockClass` in src/model/BlockType.h

**2. Create a new class inherited from `Block` for the new block.**
* *Note:* The best way to implement a new block is to look at examples of existing block classes. See the `ValveTanh` class for an example.
* Define a constructor of the form:
```
    GenericBlock(int id, Model *model) 
    : Block(id, model, BlockType::block_type, BlockClass::block_class, i
    {{Param_1, InputParameter()}, 
    {Param_2, InputParameter()}, ..., 
    {Param_N, InputParameter()}}) {}
```
  * `GenericBlock` is the name of the new class
  * `block_type` and `block_class` are the same as what was added in Step 1 above.
  * The names of the input parameters of the block are `Param_1`, ... , `Param_N`. 
  * The properties of each parameter are defined by `InputParameter`, which specifies whether it is optional, an array, a scalar, a function, and its default value.
  * The names `Param_1`, ... , `Param_N` should be the same as the parameter names within the block definition in the `.json` configuration/input file.

* The class should have a `setup_dofs(DOFHandler &dofhandler)` function.
  * This function typically only includes a call to the following function:
  ```
  Block::setup_dofs_(DOFHandler &dofhandler, int num_equations, const std::list<std::string> &internal_var_names)
  ```
  * In the above function, `num_equations` is the number of governing equations for the new block. 
  * `internal_var_names` is a list of strings that specify names for variables that are internal to the block, i.e. all variables for the block apart from the flow and pressure at the block's inlets and outlets.

* The class should have a `TripletsContributions num_triplets{*, *, *}` object. 
  * This specifies how many elements the governing equations of the block contribute to the global `F`, `E` and `dCdy` matrices respectively. Details are in Step 3 below. 

* The class should have an `update_constant` function and may also contain `update_time` and `update_solution` functions. These funtions implement the governing equations for the block. Details are below.

* (Optional) The class can have an  `enum ParamId` object that relates the parameter indices to their names. 
  * This makes it easier to reference the parameters while implementing the governing equations of the block (discussed below). 
  * The order of parameters in the `ParamId` object should match the order in the constructor. 

**3. Now implement the governing equations for the block.**
* The local state vector for each block is always arranged as `[P_in, Q_in, P_out, Q_out, InternalVariable_1, ..., InternalVariable_N]`. 
  * Here, `InternalVariable*` refers to any variable in the governing equations that are not the inlet and outlet flow and pressure. These are the same as those discussed above in the function `setup_dofs`.
  * The length of the state vector is typically four (inlet and outlet pressure and flow) plus the number of internal variables.

* The equations should be written in the form `E(t)*ydot + F(t)*y + C(y,t) = 0`. 
  * `y` is the local state vector mentioned above 
  * `ydot` is the time-derivative of the local state vector 
  * `E` and `F` are matrices of size `number_of_equations*size_of_state_vector`. 
  * `c` is a vector of length `number_of_equations`. 
  * `E` and `F` contain terms of the governing equation that multiply the respective components of `ydot` and `y` respectively.
  * `C` is a vector containing all non-linear and constant terms in the equation. 
  * If the equation contains non-linear terms, the developer should also derive the derivative of `C` with respect to `y` and `ydot`. These will be stored in the block's `dC_dy` and `dC_dydot` matrices, both of which are size `number_of_equations*size_of_state_vector`.

* For example, let's assume a block has the following non-linear governing equations:
```
a*dQ_in/dt + b*P_in + c*(dP_in/dt)*P_in + d = 0
```
```
e*dP_out/dt + f*Q_out*Q_out + g*P_out + h*I_1 = 0
```
  * For this block, the `P_in` and `Q_in` are the pressure and flow at the inlet respectively, `P_out` and `Q_out` are the pressure and flow at the outlet, and `I_1` is an internal variable.
  * The state vector is `[P_in, Q_in, P_out, Q_out, I_1]`.
  * The contributions to the local `E` matrix are `E[0,1] = a` and `E[1,2] = e`.
  * The contributions to the local `F` matrix are `F[0,0] = b`, `F[1,2] = g` and `F[1,4] = h`.
  * The contributions to the local `C` vector are `C[0] = c*(dP_in/dt)*P_in + d` and `C[1] = f*Q_out*Q_out`.
  * The contributions to the local `dC_dy` matrix are `dC_dy[0,0] = c*(dP_in/dt)` and `dC_dy[1,3] = 2*f*Q_out`.
  * The contributions to the local `dC_dydot` matrix are `dC_dydot[0,0] = c*dP_in`.

* All matrix elements that are constant are specified in `update_constant`. Matrix elements that depend only on time (not the solution of the problem itself) are specified in `update_time`. Matrix elements that change with the solution (i.e. depend on the state variables themselves) are specified in `update_solution`. Not all blocks will require the latter two functions.

* The elements of the `E`, `F`, `dC_dy` and `dC_dydot` matrices are populated using the syntax 
```
system.F.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_ids]) = a
```
  * Here, `current_block_equation_id` goes from 0 to `number_of_equations-1` (for the current block) and `current_block_variable_ids` goes from 0 to `size_of_state_vector-1` for the current block. 

* If the block contains non-linear equations, these terms must be specified in `update_solution` as `system.C(global_eqn_ids[current_block_equation_id]) = non_linear_term`.

* For non-linear equations, the derivative of the term specified above with respect to each state variable must also be provided. This goes into a `dC_dy` matrix using the following syntax 
```
system.dC_dy.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_id]) = a
```
  * Here, `a` is the derivative of the non-linear term in the equation with ID `current_block_equation_id` with respect to the local state variable with ID `current_block_variable_id`. 
  * For example, if the non-linear term is in the first equation, `current_block_equation_id = 0`. 
  * For the derivative of this term with respect to `P_in`, `current_block_variable_id = 0` and for the derivative of this term with respect to `P_out`, `current_block_variable_id = 2`.

**4. Add the new files (`GenericBlock.h` and `GenericBlock.cpp`) to src/model/CMakeLists.txt**

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
