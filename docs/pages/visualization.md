@page visualization Visualization Guide

[TOC]

# About

svZeroDVisualization is a web application designed for visualizing 0D simulation results and the 0D network. It allows users to interactively explore and analyze their simulation data through an intuitive interface. 
This application is available in the  `applications` folder.


# Architecture
svZeroDVisualization is built using a robust architecture that includes:
- Frontend: Utilizes HTML, CSS, Dash, and Plotly for creating a dynamic and interactive user interface. This setup allows for effective visualization and interaction with the 0D network and simulation results.
- Backend: Powered by a Flask application that handles the server-side logic. It leverages NetworkX for managing and visualizing the network graph and a Python version of svZeroDSolver to draw connections.


# How to Use
1. Navigate to the `applications` folder and then into the `dirgraph_visualization`.  subdirectory.

2. Open the dirgraph_main.py file.

3. In the file, locate line 47 where the dirgraph function is called. 
Pass the filepath to your input JSON file and specify the output_directory where you want the visualization to be saved.

```bash
results, parameters, G = dirgraph( filepath='' , output_dir='')
```

4. Run the script. It will execute the svZeroDSolver, generate a directed graph visualization of your network, parse simulation results, 
and display the results along with the corresponding nodes on a local Flask server.

5. Once the server is open, you can select a node you want to inspect further. The data for that node will be displayed, including the simulation parameters input for that node, pressure/flow data, and any internal variables if present.

6. Additional actions include downloading figures and using the trace feature to further inspect elements by their type (e.g., view only blood vessels or junctions in the network).