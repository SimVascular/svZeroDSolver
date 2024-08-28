@page visualization Visualization Guide

[TOC]

# About

svZeroDVisualization is a web application designed for visualizing 0D simulation results and the 0D network. It allows users to interactively explore and analyze their simulation data through an intuitive interface. 
This application is available in the  `applications` folder.


# Architecture
svZeroDVisualization is built using a robust architecture that includes:
- Frontend: Utilizes HTML, CSS, Dash, and Plotly for creating a dynamic and interactive user interface. This setup allows for effective visualization and interaction with the 0D network and simulation results.
- Backend: Powered by a Flask application that handles the server-side logic. It leverages NetworkX for managing and visualizing the network graph and a Python code to determine network connections.

# Installing Dependencies 
1. We recommend using a virtual environment to help manage project-specific 
dependencies and avoid conflicts with other projects.
- Using venv:
```bash
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```
- Using Conda:
```bash
conda create --name myenv python=3.12  # Replace with your desired Python version
conda activate myenv
```

2. Install the necessary packages:
```bash
pysvzerod
pandas
matplotlib
networkx
dash
plotly
numpy
argparse
```


# How to Use
Note: Files related to this application are in the `applications`folder, within the `dirgraph_visualization` subdirectory.


1. Command Line Execution: Pass the filepath to your input JSON file and the output_directory where you want the visualization to be saved as command line arguments. 
Pass a third argument `export_csv` optionally if you want to save svZeroDSolver raw output.
   - The program will execute svZeroDSolver, generate a directed graph visualization of your network, parse simulation results, 
      and display the results along with the corresponding nodes on a local Flask server.

```bash
python applications/svZeroDVisualization/visualize_simulation.py 'tests/cases/chamber_elastance_inductor.json' './output/circuit_img/dir_graph'
```

2. Once the server is open, you can click on a node to inspect further. 
The data for that node will be displayed, including the simulation parameters input for that node, pressure/flow data, and any internal variables if present. 
   - Additional features include the ability to download figures and use the trace function 
   for more detailed inspection of network elements. The trace feature allows users to filter the 
   view by specific element types, such as isolating and examining only the blood vessels or 
   identifying the locations of the chambers within the network. This functionality enhances the 
   ability to focus on and analyze particular components of the network with precision.

   
# How to Visualize a New Block
1. **Update JSON Parsing**:
   - When parsing the JSON file, ensure that the new block type is included. 
   - Update the `load_json_input_file` function to create a new pandas DataFrame for the new block type.

2. **Create a New Function for the Block Type**:
   - Develop a function that processes the new block type. This function should take in:
     - `d` (the dictionary of parameters loaded from the JSON file)
     - Any relevant DataFrames (e.g., `df_vessels`, `junctions_expanded`)
   - Within this function:
     - Initialize a dictionary to hold all blocks of the new type.
     - Iterate through potential connectors for the new block type.
     - For each connector, update the `connecting_block_list` and determine the `flow_directions` (use +1 for downstream and -1 for upstream).

3. **Update the Parameter Dictionary**:
   - Add the newly created blocks to the `d["blocks"]` dictionary.
   
4. **Update Existing Functions**:
- Ensure that the new block type is integrated into the existing block structure, allowing it to interact with other components in the visualization.
- For instance, if the new block type can be connected to vessels, make sure to update functions like `create_vessel_blocks` to handle connections involving the new block type. 
This includes updating any relevant functions that create or manage connections between blocks.

By following these steps, you will ensure that the new block type is properly parsed, processed, and integrated into the visualization system.
