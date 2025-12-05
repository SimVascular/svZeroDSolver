"""
Compare results between constant and time-varying coronary boundary conditions.

This script loads results from two test cases:
1. pulsatileFlow_R_coronary.json - constant microvascular resistance
2. pulsatileFlow_R_coronary_varres.json - time-varying microvascular resistance

It visualizes the effect of time-varying resistance on pressure and flow.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def compute_Ram(t_cycle, Ram_min, Ram_max, T_vc, T_vr):
    """Compute time-varying resistance Ram(t_cycle)."""
    if t_cycle <= T_vc:
        # Contraction phase
        e_t = 0.5 * (1.0 - np.cos(np.pi * t_cycle / T_vc))
    elif t_cycle <= T_vc + T_vr:
        # Relaxation phase
        e_t = 0.5 * (1.0 + np.cos(np.pi * (t_cycle - T_vc) / T_vr))
    else:
        # Rest phase
        e_t = 0.0

    sqrt_Ram_min = np.sqrt(Ram_min)
    sqrt_Ram_max = np.sqrt(Ram_max)
    Ram_t = sqrt_Ram_min + (sqrt_Ram_max - sqrt_Ram_min) * e_t
    return Ram_t * Ram_t

# Load results
results_dir = os.path.join(os.path.dirname(__file__), 'cases', 'results')
ref_constant = pd.read_json(os.path.join(results_dir, 'result_pulsatileFlow_R_coronary.json'))
ref_varres = pd.read_json(os.path.join(results_dir, 'result_pulsatileFlow_R_coronary_varres.json'))

# Extract last cardiac cycle (steady state)
# Constant resistance case
t_constant = ref_constant['time'].values
# Find last cycle
t_cycle_start = t_constant[-1] - 1.0
mask_constant = t_constant >= t_cycle_start
t_constant_last = t_constant[mask_constant] - t_cycle_start
p_in_constant = ref_constant['pressure_in'].values[mask_constant]
p_out_constant = ref_constant['pressure_out'].values[mask_constant]
q_in_constant = ref_constant['flow_in'].values[mask_constant]
q_out_constant = ref_constant['flow_out'].values[mask_constant]

# Time-varying resistance case
t_varres = ref_varres['time'].values
t_cycle_start = t_varres[-1] - 1.0
mask_varres = t_varres >= t_cycle_start
t_varres_last = t_varres[mask_varres] - t_cycle_start
p_in_varres = ref_varres['pressure_in'].values[mask_varres]
p_out_varres = ref_varres['pressure_out'].values[mask_varres]
q_in_varres = ref_varres['flow_in'].values[mask_varres]
q_out_varres = ref_varres['flow_out'].values[mask_varres]

# Compute time-varying resistance for the variable resistance case
Ram_min = 50.0
Ram_max = 200.0
T_vc = 0.3
T_vr = 0.4
Ram_varres = np.array([compute_Ram(t, Ram_min, Ram_max, T_vc, T_vr) for t in t_varres_last])

# For constant case, Ram is constant at 100 (typical value for standard coronary BC)
Ram_constant = np.ones_like(t_constant_last) * 100.0

# Create comparison plots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Effect of Time-Varying Microvascular Resistance on Coronary Flow', fontsize=14, fontweight='bold')

# Plot 1: Microvascular Resistance
ax = axes[0, 0]
ax.plot(t_constant_last, Ram_constant, 'b-', linewidth=2, label='Constant Ram')
ax.plot(t_varres_last, Ram_varres, 'r-', linewidth=2, label='Time-varying Ram')
ax.axvspan(0, T_vc, alpha=0.2, color='red', label='Contraction')
ax.axvspan(T_vc, T_vc + T_vr, alpha=0.2, color='yellow', label='Relaxation')
ax.axvspan(T_vc + T_vr, 1.0, alpha=0.2, color='green', label='Rest')
ax.set_xlabel('Time in cardiac cycle (s)', fontsize=11)
ax.set_ylabel('Ram (microvascular resistance)', fontsize=11)
ax.set_title('Microvascular Resistance Variation', fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right')
ax.set_xlim([0, 1.0])

# Plot 2: Outlet Pressure
ax = axes[0, 1]
ax.plot(t_constant_last, p_out_constant, 'b-', linewidth=2, label='Constant Ram')
ax.plot(t_varres_last, p_out_varres, 'r-', linewidth=2, label='Time-varying Ram')
ax.axvspan(0, T_vc, alpha=0.1, color='red')
ax.axvspan(T_vc, T_vc + T_vr, alpha=0.1, color='yellow')
ax.axvspan(T_vc + T_vr, 1.0, alpha=0.1, color='green')
ax.set_xlabel('Time in cardiac cycle (s)', fontsize=11)
ax.set_ylabel('Outlet Pressure (Pa)', fontsize=11)
ax.set_title('Outlet Pressure Comparison', fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim([0, 1.0])

# Plot 3: Flow Rate
ax = axes[1, 0]
ax.plot(t_constant_last, q_out_constant, 'b-', linewidth=2, label='Constant Ram')
ax.plot(t_varres_last, q_out_varres, 'r-', linewidth=2, label='Time-varying Ram')
ax.axvspan(0, T_vc, alpha=0.1, color='red')
ax.axvspan(T_vc, T_vc + T_vr, alpha=0.1, color='yellow')
ax.axvspan(T_vc + T_vr, 1.0, alpha=0.1, color='green')
ax.set_xlabel('Time in cardiac cycle (s)', fontsize=11)
ax.set_ylabel('Outlet Flow Rate (mL/s)', fontsize=11)
ax.set_title('Outlet Flow Rate Comparison', fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim([0, 1.0])

# Plot 4: Pressure-Flow Relationship
ax = axes[1, 1]
ax.plot(q_out_constant, p_out_constant, 'b-', linewidth=2, label='Constant Ram', alpha=0.7)
ax.plot(q_out_varres, p_out_varres, 'r-', linewidth=2, label='Time-varying Ram', alpha=0.7)
ax.set_xlabel('Flow Rate (mL/s)', fontsize=11)
ax.set_ylabel('Outlet Pressure (Pa)', fontsize=11)
ax.set_title('Pressure-Flow Loop', fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()

# Print summary statistics
print("=" * 70)
print("Comparison of Constant vs. Time-Varying Coronary Resistance")
print("=" * 70)
print()
print("CONSTANT RESISTANCE (Ram = 100):")
print(f"  Outlet Pressure: {p_out_constant.mean():.2f} ± {p_out_constant.std():.2f} Pa")
print(f"  Outlet Flow:     {q_out_constant.mean():.3f} ± {q_out_constant.std():.3f} mL/s")
print(f"  Pressure range:  [{p_out_constant.min():.2f}, {p_out_constant.max():.2f}] Pa")
print(f"  Flow range:      [{q_out_constant.min():.3f}, {q_out_constant.max():.3f}] mL/s")
print()
print("TIME-VARYING RESISTANCE (Ram: 50 → 200 → 50):")
print(f"  Outlet Pressure: {p_out_varres.mean():.2f} ± {p_out_varres.std():.2f} Pa")
print(f"  Outlet Flow:     {q_out_varres.mean():.3f} ± {q_out_varres.std():.3f} mL/s")
print(f"  Pressure range:  [{p_out_varres.min():.2f}, {p_out_varres.max():.2f}] Pa")
print(f"  Flow range:      [{q_out_varres.min():.3f}, {q_out_varres.max():.3f}] mL/s")
print()
print("DIFFERENCES (Time-varying - Constant):")
print(f"  Mean pressure:   {p_out_varres.mean() - p_out_constant.mean():+.2f} Pa ({(p_out_varres.mean()/p_out_constant.mean()-1)*100:+.1f}%)")
print(f"  Mean flow:       {q_out_varres.mean() - q_out_constant.mean():+.3f} mL/s ({(q_out_varres.mean()/q_out_constant.mean()-1)*100:+.1f}%)")
print(f"  Pressure std:    {p_out_varres.std() - p_out_constant.std():+.2f} Pa")
print(f"  Flow std:        {q_out_varres.std() - q_out_constant.std():+.3f} mL/s")
print()

# Phase analysis
print("PHASE-SPECIFIC ANALYSIS:")
print()

# Contraction phase (0 to T_vc)
mask_contract_const = t_constant_last <= T_vc
mask_contract_var = t_varres_last <= T_vc
print(f"CONTRACTION PHASE (0 → {T_vc}s, Ram increases):")
print(f"  Constant:     P_out = {p_out_constant[mask_contract_const].mean():.2f} Pa, Q = {q_out_constant[mask_contract_const].mean():.3f} mL/s")
print(f"  Time-varying: P_out = {p_out_varres[mask_contract_var].mean():.2f} Pa, Q = {q_out_varres[mask_contract_var].mean():.3f} mL/s")
print(f"  Effect: Pressure {(p_out_varres[mask_contract_var].mean()/p_out_constant[mask_contract_const].mean()-1)*100:+.1f}%, Flow {(q_out_varres[mask_contract_var].mean()/q_out_constant[mask_contract_const].mean()-1)*100:+.1f}%")
print()

# Relaxation phase (T_vc to T_vc + T_vr)
mask_relax_const = (t_constant_last > T_vc) & (t_constant_last <= T_vc + T_vr)
mask_relax_var = (t_varres_last > T_vc) & (t_varres_last <= T_vc + T_vr)
print(f"RELAXATION PHASE ({T_vc} → {T_vc + T_vr}s, Ram decreases):")
print(f"  Constant:     P_out = {p_out_constant[mask_relax_const].mean():.2f} Pa, Q = {q_out_constant[mask_relax_const].mean():.3f} mL/s")
print(f"  Time-varying: P_out = {p_out_varres[mask_relax_var].mean():.2f} Pa, Q = {q_out_varres[mask_relax_var].mean():.3f} mL/s")
print(f"  Effect: Pressure {(p_out_varres[mask_relax_var].mean()/p_out_constant[mask_relax_const].mean()-1)*100:+.1f}%, Flow {(q_out_varres[mask_relax_var].mean()/q_out_constant[mask_relax_const].mean()-1)*100:+.1f}%")
print()

# Rest phase (T_vc + T_vr to 1.0)
mask_rest_const = t_constant_last > T_vc + T_vr
mask_rest_var = t_varres_last > T_vc + T_vr
print(f"REST PHASE ({T_vc + T_vr} → 1.0s, Ram = {Ram_min}):")
print(f"  Constant:     P_out = {p_out_constant[mask_rest_const].mean():.2f} Pa, Q = {q_out_constant[mask_rest_const].mean():.3f} mL/s")
print(f"  Time-varying: P_out = {p_out_varres[mask_rest_var].mean():.2f} Pa, Q = {q_out_varres[mask_rest_var].mean():.3f} mL/s")
print(f"  Effect: Pressure {(p_out_varres[mask_rest_var].mean()/p_out_constant[mask_rest_const].mean()-1)*100:+.1f}%, Flow {(q_out_varres[mask_rest_var].mean()/q_out_constant[mask_rest_const].mean()-1)*100:+.1f}%")
print()
print("=" * 70)

output_path = os.path.join(os.path.dirname(__file__), 'coronary_bc_comparison.png')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print()
print(f"Plot saved to: {output_path}")
print()

# Uncomment the line below to display the plot
# plt.show()
