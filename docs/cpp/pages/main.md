@mainpage svZeroDPlus

[TOC]

svZeroDPlus is an application for performing different tasks with 0D hemodynamic
computer models. It is written in C++ to enable the
highest performance. It also offers a Python API for easily integrating it into
other routines.

* <a href="https://github.com/StanfordCBCL/svZeroDPlus">Source
repository</a>
* <a href="https://simvascular.github.io">About SimVascular</a>

Zero-dimensional (0D) models
are lightweight methods to simulate bulk hemodynamic quantities in the
cardiovascular system. Unlike 3D and 1D models, 0D models are purely
time-dependent; they are unable to simulate spatial patterns in the
hemodynamics. 0D models are analogous to electrical circuits. The flow rate
simulated by 0D models represents electrical current, while the pressure
represents voltage. Three primary building blocks of 0D models are resistors,
capacitors, and inductors Resistance captures the viscous effects of blood
flow, capacitance represents the compliance and distensibility of the vessel
wall, and inductance represents the inertia of the blood flow. Different
combinations of these building blocks, as well as others, can be formed to
reflect the hemodynamics and physiology of different cardiovascular
anatomies. These 0D models are governed by differential algebraic equations
(DAEs).

For more background information on 0D models, have a look at SimVascular's
[ROM Simulation Guide](http://simvascular.github.io/docsROMSimulation.html).

# Installation

svZeroDPlus can be installed in two different ways. For using the Python
API, an installation via pip is recommended.

## Using pip

For a pip installation, simply run the following command
(cloning of the repository is not required):

```bash
pip install git+https://github.com/StanfordCBCL/svZeroDPlus.git
```

## Using CMake

If you want to build svZeroDPlus manually from source, clone the repository
and run the following commands from the top directory of the project:

```bash
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```


@remark <details>
  <summary>**Building on Sherlock**</summary>

```bash
module load cmake/3.23.1 gcc/12.1.0 binutils/2.38
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/share/software/user/open/gcc/12.1.0/bin/g++ -DCMAKE_C_COMPILER=/share/software/user/open/gcc/12.1.0/bin/gcc ..
cmake --build .
```

</details>

# Getting started

* [Running simulations with svZeroDSolver](@ref svzerodsolver)
* [Calibrating 0D models with svZeroDCalibrator](@ref svzerodcalibrator)
* [Embedding svZeroDPlus in your code using the Python API](@ref svzerodplus)

# Developer Guide

If you are a developer and want to contribute to svZeroDSolver, you can find
more helpful information in our [Developer Guide](@ref developer_guide).
