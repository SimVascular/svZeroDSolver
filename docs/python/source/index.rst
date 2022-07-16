.. _mainpage:

******************************
svZeroDSolver documentation
******************************

**Version**: |release|

**Useful links**:
`Source Repository <https://github.com/SimVascular/svZeroDSolver>`_ |
`Issue Tracker <https://github.com/SimVascular/svZeroDSolver/issues>`_ |
`SimVascular <https://simvascular.github.io>`_

svZeroDSolver is a Python code that simulates the hemodynamics in zero-dimensional
(0D) lumped parameter models of vascular networks. These 0D models are governed
by differential algebraic equations (DAEs).

The solver uses a highly modular framework to model the vascular anatomy, using
individual 0D elements to represent different parts of the vascular anatomy
(and boundary conditions). The individual 0D elements and their associated
governing equations defined in `blocks.py`. 
In `solver.py`, the blocks are assembled and simulated using the generalized-alpha
time-stepping method defined in `time_integration.py`.

The svZeroDSolver Python files are in the `svzerodsolver` Python package directory. 

svZeroDSolver currently supports the following vascular 0D modeling options and
boundary conditions:

**Vascular 0D elements**:

* Resistor
* Resistor-capacitor
* Resistor-inductor
* Resistor-capacitor-inductor

**Boundary conditions**:

* Pressure
* Resistor
* RCR
* Coronary
* Flow

Installation
************

svZeroDSolver and all its dependencies can be installed easily via pip.

.. code-block:: bash

   pip install git+https://github.com/SimVascular/svZeroDSolver.git


For Contributers
================

The following guide provides all necessary steps to install your local
svZeroDSolver repository via pip in editable mode to allow for local code changes
to reflect in the package. 

If you are contributing to svZeroDSolver, it is highly recommended to use a virtual
environment like [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
After installing Miniconda you can create a new environment and enter it using:

.. code-block:: bash

   conda create -n zerodsolver python=3.9
   conda activate zerodsolver

After that, enter the repository folder and install the svZeroDSolver
**with development related dependencies** using:

.. code-block:: bash

   pip install -e .[dev]


*If you are using the `zsh` shell, enter: `pip install -e ".[dev]"`*

Usage
*****

Command line
============

To run svZeroDSolver form the command line, run:

.. code-block:: bash

   zerod SOLVER_INPUT_FILE 


For more information about command line options, enter:

.. code-block:: bash

   zerod --help

As a python module
==================

.. code-block:: python

   import svzerodsolver
   svzerodsolver.solver.set_up_and_run_0d_simulation('input.json')

This variant enables running svZeroDSolver within a user-defined Python code
(e.g. parameter optimization, uncertainty quantification)

More information
****************

.. toctree::
   :maxdepth: 1

   Documentation <modules>
   References <rst/bibliography>
