# Activation Functions for Chamber Models

## Overview

Chamber models in svZeroDSolver support different activation functions to model the time-varying elastance of cardiac chambers. The activation function determines how the chamber contracts and relaxes over the cardiac cycle.

## Available Activation Functions

### 1. Half Cosine Activation

A half cosine wave activation function to model contraction:

```
A(t) = -0.5 * cos(2π * t_contract / t_twitch) + 0.5  if t_contract ≤ t_twitch
A(t) = 0  otherwise
```

where `t_contract = max(0, t - t_active)`

**Parameters:**
- `t_active`: Time when activation begins within the cardiac cycle
- `t_twitch`: Duration of the contraction twitch

**Example usage in JSON:**
```json
{
  "type": "ChamberElastanceInductor",
  "name": "ventricle",
  "values": {
    "Emax": 1.057,
    "Emin": 0.091,
    "Vrd": 26.1,
    "Vrs": 18.0,
    "Impedance": 0.000351787,
    "activation_function_type": "half_cosine",
    "activation_function_values": {
      "t_active": 0.2,
      "t_twitch": 0.3
    }
  }
}
```

### 2. Piecewise Cosine Activation

A piecewise activation function that models separate contraction and relaxation phases:

```
φ(t) = 0.5 * [1 - cos(π * (t - t_C) / T_C)]  during contraction
φ(t) = 0.5 * [1 + cos(π * (t - t_R) / T_R)]  during relaxation
φ(t) = 0  otherwise
```

**Parameters:**
- `contract_start`: Time when contraction starts
- `relax_start`: Time when relaxation starts
- `contract_duration`: Duration of contraction phase
- `relax_duration`: Duration of relaxation phase

**Example usage in JSON:**
```json
{
  "type": "LinearElastanceChamber",
  "name": "left_atrium",
  "values": {
    "Emax": 199.95,
    "Epass": 260.59,
    "Vrest": 26.24,
    "activation_function_type": "piecewise_cosine",
    "activation_function_values": {
      "contract_start": 0.025,
      "relax_start": 0.08625,
      "contract_duration": 0.06125,
      "relax_duration": 0.18375
    }
  }
}
```

### 3. Two Hill Activation

A more flexible and physiologically realistic activation function based on the two-hill model. This allows for more precise control over the activation waveform shape.

The activation is computed as:
```
A(t) = C * [g₁(t) / (1 + g₁(t))] * [1 / (1 + g₂(t))]
```

where:
```
g₁(t) = (t_shifted / τ₁)^m₁
g₂(t) = (t_shifted / τ₂)^m₂
t_shifted = (t - t_shift) mod T_cardiac
```

C is a normalization constant ensuring maximum activation is 1.

**Parameters:**
- `t_shift`: Time shift parameter to align activation with cardiac cycle
- `tau_1`: Time constant for first hill (controls rise)
- `tau_2`: Time constant for second hill (controls fall)
- `m1`: Exponent for first hill (controls steepness of rise)
- `m2`: Exponent for second hill (controls steepness of fall)

**Example usage in JSON:**
```json
{
  "type": "ChamberElastanceInductor",
  "name": "ventricle",
  "values": {
    "Emax": 1.057,
    "Emin": 0.091,
    "Vrd": 26.1,
    "Vrs": 18.0,
    "Impedance": 0.000351787,
    "activation_function_type": "two_hill",
    "activation_function_values": {
      "t_shift": 0.15,
      "tau_1": 0.25,
      "tau_2": 0.35,
      "m1": 1.5,
      "m2": 10.0
    }
  }
}
```

## Switching Between Activation Functions

**Both `ChamberElastanceInductor` and `LinearElastanceChamber` require the `activation_function_type` parameter to be explicitly specified as a string:**

- `activation_function_type: "half_cosine"` - Half Cosine activation
- `activation_function_type: "piecewise_cosine"` - Piecewise Cosine activation
- `activation_function_type: "two_hill"` - Two Hill activation

All activation function-specific parameters must be grouped under `activation_function_values` as a nested dictionary.

### Important Notes

- The `activation_function_type` parameter is **required** for both chamber types
- Only string values are accepted (e.g., "half_cosine", not numeric values)
- Parameters must be specified within the `activation_function_values` nested dictionary

## Elastance Calculation

For all activation functions, the time-varying elastance E(t) is computed as:

**For ChamberElastanceInductor:**
```
E(t) = (Emax - Emin) * A(t) + Emin
Vrest(t) = (1 - A(t)) * (Vrd - Vrs) + Vrs
```

**For LinearElastanceChamber:**
```
E(t) = Epass + Emax * φ(t)
```

## References

The two-hill activation function is described in:
- Kaiser, A. D., et al. (2022). "A design-based model of the aortic valve for fluid-structure interaction." Biomechanics and Modeling in Mechanobiology. https://link.springer.com/article/10.1007/s10439-022-03047-3
