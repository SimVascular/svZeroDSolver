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

A user guide is available on the [SimVascular website](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-visualization).
   
# How to Visualize a New Block
1. **Update JSON Parsing**:
   - When parsing the JSON file, ensure that the new block type is included. 
   - Update the `load_json_input_file` function to create a new pandas DataFrame for the new block type.

2. **Create a New Function for the %Block Type**:
   - Develop a function that processes the new block type. This function should take in:
     - `d` (the dictionary of parameters loaded from the JSON file)
     - Any relevant DataFrames (e.g., `df_vessels`, `junctions_expanded`)
   - Within this function:
     - Initialize a dictionary to hold all blocks of the new type.
     - Iterate through potential connectors for the new block type.
     - For each connector, update the `connecting_block_list` and determine the `flow_directions` (use +1 for downstream and -1 for upstream).

3. **Update the %Parameter Dictionary**:
   - Add the newly created blocks to the `d["blocks"]` dictionary.
   
4. **Update Existing Functions**:
   - Ensure that the new block type is integrated into the existing block structure, allowing it to interact with other components in the visualization.
   - For instance, if the new block type can be connected to vessels, make sure to update functions like `create_vessel_blocks` to handle connections involving the new block type. 
This includes updating any relevant functions that create or manage connections between blocks.

By following these steps, you will ensure that the new block type is properly parsed, processed, and integrated into the visualization system.
