import pandas as pd
import pdb
import matplotlib.pyplot as plt
import numpy as np

out_name = "out.csv"
res = pd.read_csv(out_name)

variables = ['flow:inlet_valve:ventricle',
             'flow:ventricle:outlet_valve',
             'pressure:inlet_valve:ventricle',
             'pressure:ventricle:outlet_valve',
             'r:ventricle',
             'v:ventricle',
             'S:ventricle',
             'tau:ventricle',
             'V:ventricle']

# Conversion factors for each variable
# Flows: m^3/s -> mL/s (multiply by 1e6)
# Pressures: Pa -> mmHg (divide by 133.322)
# Radius: m -> cm (multiply by 100)
# Velocity: m/s -> cm/s (multiply by 100)
# Stresses: Pa -> mmHg (divide by 133.322)
conversion_factors = [1e6, 1e6, 1/133.322, 1/133.322, 100, 100, 1/133.322, 1/133.322, 1e6]

y_axis_titles = ['Inlet Flow (mL/s)', 'Outlet Flow (mL/s)', 'Inlet Pressure (mmHg)',
                 'Outlet Pressure (mmHg)', 'Radius (cm)', 'Velocity (cm/s)',
                 'Stress (mmHg)', 'Active Stress (mmHg)', 'Volume (mL)']

fig, axs = plt.subplots(6, 2, figsize=(12, 24))
axs = axs.ravel()

for idx, var in enumerate(variables):
    name = f"{var}"
    ids = res.name == name
    time = np.array(res[ids].time)
    y = np.array(res[ids].y) * conversion_factors[idx]  # Apply conversion factor
    axs[idx].plot(time, y, label=var)
    axs[idx].set_xlabel('Time (s)')
    axs[idx].set_ylabel(y_axis_titles[idx])
    axs[idx].legend(y_axis_titles[idx])

# --- Net flow plot (outlet - inlet) ---
inlet_ids = res.name == 'flow:inlet_valve:ventricle'
outlet_ids = res.name == 'flow:ventricle:outlet_valve'

time = np.array(res[inlet_ids].time)
inlet_flow = np.array(res[inlet_ids].y) * 1e6  # Convert to mL/s
outlet_flow = np.array(res[outlet_ids].y) * 1e6  # Convert to mL/s
net_flow = outlet_flow - inlet_flow

axs[9].plot(time, net_flow, label='net flow (outlet - inlet)')
axs[9].set_xlabel('Time (s)')
axs[9].set_ylabel('Net Flow')
axs[9].legend()

# --- PV loop plot ---
pressure_ids = res.name == 'pressure:ventricle:outlet_valve'
volume_ids = res.name == 'V:ventricle'

volume = np.array(res[volume_ids].y) * 1e6  # Convert to mL
pressure = np.array(res[pressure_ids].y) / 133.322  # Convert to mmHg

axs[10].plot(volume, pressure, color='red')
axs[10].set_xlabel('Volume (mL)')
axs[10].set_ylabel('Pressure (mmHg)')
axs[10].set_title('PV Loop')
#axs[10].grid(True)


plt.tight_layout()
plt.savefig('output_plot.pdf')




# name = "tau:ventricle"
# ids = res.name == name
# out = np.array(res[ids])
# plt.plot(np.array(res[ids].time), np.array(res[ids].y), label=name)
# plt.show()

# # # Load out.csv for time values
# out_name = "out.csv"
# res_out = pd.read_csv(out_name)

# # Load screen_out.csv for y values (assume it has no titles)
# screen_out_name = "screen_out.csv"
# res_screen = pd.read_csv(screen_out_name, header=None)  # Load without headers

# # Extract time values from out.csv
# time_name = "tau:ventricle"
# time_ids = res_out.name == time_name
# time_values = np.array(res_out[time_ids].time)

# # Extract y values directly from the first column of screen_out.csv
# y_values = np.array(res_screen.iloc[0::2,0])  # Use the first column and every other row

# # Plot time (x-axis) vs y (y-axis)
# plt.plot(time_values, y_values, label="Screen Output")
# plt.xlabel("Time")
# plt.ylabel("f")
# plt.legend()
# plt.show()