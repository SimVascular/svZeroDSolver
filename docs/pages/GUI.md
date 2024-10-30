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

# User guide

A user guide is available on the [SimVascular website](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-gui). 
