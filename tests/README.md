# Testing Guide

[Integration testing](https://en.wikipedia.org/wiki/Integration_testing) is an essential part of software development. It is performed when integrating code changes into the main development branch to verify that the code works as expected. The following sections describe how to run and add integration tests used to the svMultiPhysics program.


# Build svZeroDSolver

Running a test case requires you to [build svZeroDSolver](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-install) and install with pip in order to run test cases using `pytest`

# Running tests using pytest

to run the full set of svZeroDSolver test cases, navigate to the `tests/` directory and run `pytest` from the terminal

## Create a new solver test

In order to create a new test for the solver, complete the following steps:

1. Create a svZeroDSolver json file for the test (e.g. `pulsatileFlow_R_RCR.json`) in the `tests/cases` directory
2. Verify the results. If possible, verify against an analystical solution. If this is not feasible, find some way to verify that the test results are correct. For example, you might use reference solutions from other codes, manufactured solutions, convergence analysis, or something of that variety.
3. Compute a reference solution for the test by navigating to the `tests` directory and running `python compute_ref_sol.py <your_test_case_name..json>` from the terminal. This should create a resultant `result_your_test_case_name.json` file in the `tests/cases/results` directory. Verify the result is there.
4. Add your test case to the `pytest` parametrization in `test_solver.py`. simply add the name of your test `json` to the list of strings after the `@pytest.mark.parametrize` decorator.
5. Verify that the test runs properly by running `pytest` from the terminal while in the `tests` directory.



