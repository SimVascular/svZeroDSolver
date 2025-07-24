@page jacobian Jacobian Generator for svZeroDSolver

This tool generates C++ code for block implementations in the svZeroDSolver framework using symbolic mathematics.

## Overview

The script `script/jacobian.py` reads block definitions from YAML files and generates C++ code for implementing the mathematical models in the solver. It uses symbolic differentiation (via SymPy) to automatically derive the necessary Jacobian matrices and system contributions.

## Usage

```bash
python jacobian.py <yaml_file>
```

## YAML File Format

The YAML files define the mathematical model of a block with the following structure:

### Required Sections

- `variables`: List of variable names in the model
- `derivatives`: List of derivative names (must match variables with `_dt` suffix)
- `constants`: List of parameter names
- `residuals`: List of residual equations that define the system

### Optional Sections

- `time_dependent`: List of parameters that depend on time (e.g., activation functions)
- `helper_functions`: Python code defining helper functions used in the residuals

### Example

```yaml
variables:
  - Pin
  - Qin
  - Pout
  - Qout

derivatives:
  - dPin_dt
  - dQin_dt
  - dPout_dt
  - dQout_dt

constants:
  - R
  - C
  - L
  - S

residuals:
  - Pin - Pout - (R + S * abs(Qin)) * Qin - L * dQout_dt
  - Qin - Qout - C * dPin_dt + C * (R + 2 * S * abs(Qin)) * dQin_dt
```

## Output

The script generates three C++ function implementations:

1. `update_constant` - Sets up constant matrix coefficients for the system
2. `update_time` - Updates time-dependent parameters
3. `update_solution` - Computes solution-dependent terms and Jacobians

## Workflow

1. Create a YAML file defining your block's mathematical model
2. Run `jacobian.py` on this file to generate C++ code
3. Copy the generated code to your block implementation file
4. Complete the implementation with necessary boilerplate code

## Tips

- Ensure the number of variables equals the number of derivatives
- The number of residuals should be equal to the number of variables minus 2
- Use helper functions for complex expressions to improve readability
- Define time-dependent constants separately

## Examples

See the provided YAML examples:
- `ChamberSphere.yaml` - Spherical heart chamber model
- `BloodVessel.yaml` - Blood vessel model with optional stenosis
- `ClosedLoopCoronaryBC.yaml` - Coronary boundary condition model