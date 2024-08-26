@page GUI GUI Guide

[TOC]

# About
The svZeroDGUI application is designed to facilitate the creation of 0D model input files 
through an intuitive graphical user interface. Located in the `applications` folder, 
this tool allows users to generate input files for the svZeroDSolver by visually 
drawing and configuring the network. 

Unlike manual file creation, which can be 
cumbersome and error-prone, svZeroDGUI provides an easy-to-use interface that 
simplifies the process of defining network components such as vessels, junctions, and 
boundary conditions. This application is especially valuable for users who lack access to 
3D models or seek an efficient alternative to manual file generation, making the model creation 
process both faster and more user-friendly.


# Architecture

svZeroDGUI is built using a robust architecture that includes:
* Frontend: The frontend is developed with HTML, CSS, and JavaScript to create a 
responsive and user-friendly interface. It utilizes Cytoscape.js, a popular package for creating
interactive elements and graphical networks.

*  Backend: Flask app, Node.js for server-side logic, and Cypress for testing.
This architecture supports an intuitive user experience for 
generating and managing 0D input files through a graphical interface.


# How to Use
1. Navigate to the `applications` folder and then to the `create_0dmodel` subdirectory.
2. Launch the `app.py` file.
```bash
python applications/svZeroDGUI/app.py
```
3. Select a node type and name the node.
- For vessels, after drawing the node, click on it to open a form 
where you can enter details such as vessel length, diameter, and more.
- For junctions, click the node to specify if it’s a Normal Junction 
or a Blood Vessel junction.
4. To draw edges between nodes, toggle the `Draw on` button on the right. 
Once active, you can start connecting nodes by drawing edges between them.
5. When you wish to stop drawing edges and continue adding or moving nodes, 
click the `Draw off` button.
6. Once you’ve completed the network, click `Export to JSON` on the right. 
If there are any incorrect connections or patterns, an alert will prompt you 
to make necessary changes so the network can be processed by svZeroDSolver.
7. Open the downloaded JSON file and add any additional information, 
such as boundary condition data, before running it through svZeroDSolver.
