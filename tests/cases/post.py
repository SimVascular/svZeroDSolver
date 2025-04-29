import pandas as pd
import pdb
import matplotlib.pyplot as plt
import numpy as np

out_name = "out.csv"
res = pd.read_csv(out_name)

variables = ['flow:INFLOW:ventricle',
             'flow:ventricle:J0',
             'pressure:INFLOW:ventricle',
             'pressure:ventricle:J0',
             'r:ventricle',
             'v:ventricle',
             'S:ventricle',
             'tau:ventricle',
             'V:ventricle']
fig, axs = plt.subplots(5, 2, figsize=(12, 15))
axs = axs.ravel()

for idx, var in enumerate(variables):
    name = f"{var}"
    ids = res.name == name
    axs[idx].plot(np.array(res[ids].time), np.array(res[ids].y), label=var)
    axs[idx].set_xlabel('Time')
    axs[idx].set_ylabel(var)
    axs[idx].legend()

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