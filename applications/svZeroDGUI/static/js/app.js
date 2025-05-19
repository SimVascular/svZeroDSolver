// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause
/**
 * app.js
 *
 * This script handles the interactive functionality of the svZeroDGUI application,
 * including event listeners for UI elements such as collapsible sections,
 * form inputs, and dynamic content toggling.
 */

// Toggles visibility of the boundary condition type container based on selected node type.
function updateBoundaryConditionTypeVisibility() {
    const nodeType = document.getElementById('node-type').value;
    const boundaryConditionTypeContainer = document.getElementById('boundary-condition-type-container');
    if (nodeType === 'boundary_condition') {
        boundaryConditionTypeContainer.style.display = 'block';
        // Scroll to reveal the container
        requestAnimationFrame(() => {
          const dropdownHeight = boundaryConditionTypeContainer.offsetHeight;
          window.scrollBy(0, dropdownHeight);
        });
    } else {
        boundaryConditionTypeContainer.style.display = 'none';
    }
}

/*
 * Defines the file paths for the icons used in the application.
 * These icons represent different components within the model, such as vessels,
 * valves, chambers, junctions, and various types of boundary conditions.
 */
const vessel_icon = '/static/css/vessel.png';
const valve_icon = '/static/css/valve.png';
const chamber_icon = '/static/css/chamber.png';
const junction_icon = '/static/css/junction.png';
const resistance_icon = '/static/css/resistance.png';
const pressure_icon = '/static/css/pressure.png';
const RCR_icon = '/static/css/RCR.png';
const coronary_icon = '/static/css/coronary.png';
const flow_icon = '/static/css/flow.png';


document.addEventListener('DOMContentLoaded', function() {
    let simulation_parameters_dict = {}
    var deleteMode = false;

     document.querySelector('#submitJunctionButton').addEventListener('click', function() {
            submitJunctionInfo();
        });

     document.querySelector('#submitVesselButton').addEventListener('click', function() {
            submitVesselInfo();
        });

     document.querySelector('#submitSimParamButton').addEventListener('click', function() {
            submitSimParameters();
        });

    var node_count = {'boundary_condition': 0, 'vessel':0, 'valve': 0, 'chamber': 0, 'junction': 0};
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
        layout: {
            name: 'grid',
            rows: 2,
            cols: 2
        },
        style: [
            {
                selector: 'node[name]',
                style: {
                    'content': 'data(name)'
                }
            }
        ]
        });

    setSimParameters(); // Opens the Parameters form when you first load the page.
    cy.domNode(); // Register the domNode extension for making HTML classes for each node to assist Cypress testing

    cy.on('tap', function(event) {
       if (event.target === cy) {  // Filters clicks to those within the cytoscape drawing container.
            let nodeName = document.getElementById('node-name').value.trim();
            let nodeType = document.getElementById('node-type').value;
            let boundaryConditionType = nodeType === 'boundary_condition' ? document.getElementById('boundary-condition-type').value : null;

            if (nodeName !== "") {
                let pos = event.position || event.cyRenderedPosition;
                let nodeId = `${nodeType}-${node_count[nodeType]}`;
                node_count[nodeType] += 1;

                let color;
                if (nodeType === 'boundary_condition') {
                    // Sets the icons and edge color for each node
                    switch (boundaryConditionType) {
                        case 'FLOW':
                            color = '#FF00FF';  // Magenta
                            node_icon = flow_icon;
                            break;
                        case 'RESISTANCE':
                            color = 'purple';
                            node_icon = resistance_icon;
                            break;
                        case 'PRESSURE':
                            color = 'orange';
                            node_icon = pressure_icon;
                            break;
                        case 'RCR':
                            color = '#ADD8E6';  // Light Blue
                            node_icon = RCR_icon;
                            break;
                        case 'CORONARY':
                            color = '#800020';  // Burgundy
                            node_icon = coronary_icon;
                            break;
                        default:
                            color = 'grey';  // Default color for unknown types
                            break;
                    }
                } else {
                    switch (nodeType) {
                        case 'vessel':  // Cardinal Red
                            color = '#C41E3A';
                            node_icon = vessel_icon;
                            break;
                        case 'valve':
                            color = 'black';
                            node_icon = valve_icon; // https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.vectorstock.com%2Froyalty-free-vector%2Fheart-valve-disease-line-icon-vector-47167009&psig=AOvVaw39gwsTXwtIDM9KxgEODvmZ&ust=1721929537132000&source=images&cd=vfe&opi=89978449&ved=0CBEQjRxqFwoTCLi788edwIcDFQAAAAAdAAAAABAE
                            break;
                        case 'chamber':
                            color = 'pink';
                            node_icon = chamber_icon;
                            break;
                        case 'junction':
                            color = '#046791';  // Dark Blue
                            node_icon = junction_icon;
                            break;
                        default:
                            color = 'black';
                            node_icon = None;
                            break;
                    }
            }
            let div = document.createElement("div");
            div.innerHTML = "";
            div.classList.add('draggable'); // Adds custom class so the nodes can be tracked in Cypress E2E testing
            div.style.width = '30px';
            div.style.height = '30px';
            let newNode = cy.add({
                   group: 'nodes',
                   data: { id: nodeId, type: nodeType, name: nodeName,
                   cls_name: boundaryConditionType || nodeType, dom: div},
                   position: { x: pos.x, y: pos.y },
                   style: {'background-image': `url(${node_icon})`,
                    'background-fit': 'cover', // Ensure the image covers the node
                    'background-opacity': 1,   // Ensure the background is fully opaque
                    'border-color': color,     // Example border color
                    'border-width': 2,         // Example border width
                    'background-color': 'white' // Ensure background color does not override image
                    }
            });
            console.log('Node created:', newNode.json());
            }
            else {
                alert("Please enter a node name.");
            }
        }
    });


    var eh = cy.edgehandles({
        enabled: false // Initializes the edge handle as turned off until user clicks draw on button
    });

   // Enables edge drawing when the user clicks the draw-on button
   document.querySelector('#draw-on').addEventListener('click', function() {
        eh.enableDrawMode();
        eh.enable();
   });

    // Disables edge drawing when the user clicks the draw-off button
    document.querySelector('#draw-off').addEventListener('click', function() {
        eh.disableDrawMode();
        eh.disable();
    });

    // Function to calculate R, C, L
    function calculateRCL(additionalData) {
         console.log(additionalData);
         let vesselLength = additionalData.vessel_length;
         let vesselRadius = additionalData.vessel_radius;

         if (!vesselLength || !vesselRadius) {
            console.error("Invalid length or radius");
            return { R: null, C: null, L: null };
         }

        const R = (8 * additionalData.mu * vesselLength) / (Math.PI * Math.pow(vesselRadius, 4));
        console.log(R);
        const C = (3 * vesselLength * Math.PI * Math.pow(vesselRadius, 3)) / (2 * additionalData.young_mod * additionalData.h);
        const L = (additionalData.rho * vesselLength) / (Math.PI * Math.pow(vesselRadius, 2));
        return { vesselLength, R, C, L };
    }

    // Exports the 0d network to a json file
    function exportGraphData() {
        let detected_objects = {
            simulation_parameters: simulation_parameters_dict,
            boundary_conditions: [],
            junctions: [],
            vessels: [],
            valves: [],
            chambers: []
        };

        cy.nodes().forEach(node => {
        let data = node.data();
        let additionalData = data.additional_data || {};
        switch(data.type) {
            case 'boundary_condition':
                let boundaryCondition = {
                    bc_name: data.name,
                    bc_type: data.cls_name,
                    bc_values: {}
                }

                switch (data.cls_name) {  // Assuming cls_name is the field that specifies the boundary condition type
                    case "FLOW":
                        boundaryCondition.bc_values = {"Q": [], "t":[]};
                        break;
                    case "RESISTANCE":
                        boundaryCondition.bc_values = {"Pd": '', "R":''};
                        break;
                    case "PRESSURE":
                        boundaryCondition.bc_values = {"P": [], "t": []};
                        break;
                    case "RCR":
                        boundaryCondition.bc_values = {"C": '', "Pd":'', "Rd":'', "Rp":''};
                        break;
                    default:
                        boundaryCondition.bc_values = {
                            "Ca": '', "Cc":'' , "Pim": [], "P_v": '',
                            "Ra1": '', "Ra2": '' , "Rv1": '', "t": []
                        };
                        break;
                }

                detected_objects.boundary_conditions.push(boundaryCondition);
                break;
            case 'vessel':
                let inlet = getInlet(node);
                let outlet = getOutlet(node);

                // Calculate R, C, L
                let {vesselLength, R, C, L } = calculateRCL(additionalData);

                let vesselObject = {
                    boundary_conditions: {},
                    vessel_id: parseInt(data.id.replace('vessel-', '')),
                    vessel_length: vesselLength || "",
                    vessel_name: data.name,
                    zero_d_element_type: "BloodVessel",
                    zero_d_element_values: {
                        C: C || 0 ,
                        L: L || 0,
                        R_poiseuille: R || 0,
                        stenosis_coefficient: additionalData.stenosis_diameter || 0
                    }
                };

                if (inlet || outlet) {
                    if (inlet) {
                        vesselObject.boundary_conditions.inlet = inlet;
                    }
                    if (outlet) {
                        vesselObject.boundary_conditions.outlet = outlet;
                    }
                } else {
                    delete vesselObject.boundary_conditions;  // Remove the empty boundary_conditions object if neither inlet nor outlet is present
                }

                detected_objects.vessels.push(vesselObject);
                break;
            case 'valve':
                detected_objects.valves.push({
                    type: 'ValveTanh',
                    name: data.name,
                    params: {
                        Rmax: "",
                        Rmin: "",
                        Steepness: "",
                        upstream_block: getUpstream(node),
                        downstream_block: getDownstream(node)
                    }
                });
                break;
            case 'chamber':
                detected_objects.chambers.push({
                    name: data.name,
                    type: 'ChamberElastanceInductor',
                    values: {
                        Emax: '',
                        Emin: '',
                        Vrd: '',
                        Vrs: '',
                        t_active: '',
                        t_twitch: '',
                        Impedance: ''
                    }
                });
                break;
            case 'junction':
                 let junction_type = additionalData.junction_type;
                 let junctionObject = {
                    inlet_vessels: getInletVessels(node),
                    junction_name: data.name,
                    junction_type: junction_type,
                    junction_values: {
                        R_poiseuille: [],
                        L: [],
                        stenosis_coefficient: []
                    },
                    outlet_vessels: getOutletVessels(node)
                };

                if (junction_type == 'BloodVesselJunction') {
                    junctionObject.junction_values.R_poiseuille.push(additionalData.R1,additionalData.R2);
                    junctionObject.junction_values.L.push(additionalData.L1, additionalData.L2);
                    junctionObject.junction_values.stenosis_coefficient.push(additionalData.junction_stenosis_coefficient1,
                                                    additionalData.junction_stenosis_coefficient2);
                } else {
                    delete junctionObject.junction_values;  // Remove the empty junction_values object if it is a normal junction
                }
                detected_objects.junctions.push(junctionObject);
                break;
            }
        });
        // If any of the lists are empty, delete them from the json file.
        if (detected_objects.junctions.length === 0) delete detected_objects.junctions;
        if (detected_objects.vessels.length === 0) delete detected_objects.vessels;
        if (detected_objects.valves.length === 0) delete detected_objects.valves;
        if (detected_objects.chambers.length === 0) delete detected_objects.chambers;
        return detected_objects;
    }

    // Returns the boundary condition inlet if present.
    function getInlet(node) {
        let inboundEdges = node.incomers('edge');
        if (inboundEdges.length > 0) {
            if (inboundEdges[0].source().data('id').startsWith('boundary_condition')) {
                return inboundEdges[0].source().data('name');
            }
        }
        return;
    }

    // Returns the boundary condition outlet if present.
    function getOutlet(node) {
        let outboundEdges = node.outgoers('edge');
        if (outboundEdges.length > 0) {
            if (outboundEdges[0].target().data('id').startsWith('boundary_condition')) {
                return outboundEdges[0].target().data('name');
            }
        }
        return;
    }

    // Returns the upstream node if present.
    function getUpstream(node) {
        let inboundEdges = node.incomers('edge');
        if (inboundEdges.length > 0) {
            return inboundEdges[0].source().data('name');
        }
        return;
    }

    // Returns the downstream node if present.
    function getDownstream(node) {
        let outboundEdges = node.outgoers('edge');
        if (outboundEdges.length > 0) {
            return outboundEdges[0].target().data('name');
        }
        return;
    }

    // Returns the inlet vessels to a node.
    function getInletVessels(node) {
        let inboundEdges = node.incomers('edge');
        let inletVessels = [];

        inboundEdges.forEach(edge => {
            let sourceNode = edge.source();
            let type = sourceNode.data('type');
            console.log("id", sourceNode.data('id'));
            if (type === 'vessel') {
                let vessel_id = parseInt(sourceNode.data('id').replace('vessel-', ''));
                inletVessels.push(vessel_id);
            }
        });
        return inletVessels;
    }

    // Returns the outlet vessels to a node.
    function getOutletVessels(node) {
          let outboundEdges = node.outgoers('edge');
          let outletVessels = [];

          outboundEdges.forEach(edge => {
                let targetNode = edge.target();
                let type = targetNode.data('type');
                if (type === 'vessel') {
                    let vessel_id = parseInt(targetNode.data('id').replace('vessel-', ''));
                    outletVessels.push(vessel_id);
                }
            });
          return outletVessels;
    }

    // Handles the downloading logic.
    function downloadJSON(data, filename = 'graph_data.json') {
        let blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
        let url = URL.createObjectURL(blob);
        let a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
    }
    window.downloadJSON = downloadJSON;

    // Returns how many incoming connections there are to a node.
    function consistency_inlet(node) {
        let inboundEdges = node.incomers('edge');
        return inboundEdges.length;
    }

    // Returns how many outgoing connections there are to a node.
    function consistency_outlet(node) {
        let outboundEdges = node.outgoers('edge');
        return outboundEdges.length;
    }

    // Debugger for creating 0d models. Ensures all proper connections are present within the drawn graph.
    function Consistency_Check(data) {
        let alerts = [];
        let inletCounts = {};
        let outletCounts = {};
        let oneToOneJunctions = {};
        let junctionInletsMap = {};

        cy.nodes().forEach(node => {
            let data = node.data();
            let type = data['type'];
            let name = data['name'];
            let id = data['id'];

            switch (type) {
                        case 'vessel':
                            if (consistency_inlet(node) + consistency_outlet(node) != 2) {
                                alerts.push(`Vessel ${name} does not have exactly two connections`);
                            }
                            break;
                        case 'boundary_condition':
                            if (consistency_inlet(node) + consistency_outlet(node) != 1) {
                                alerts.push(`Boundary condition ${name} does not have exactly one connection`);
                            }
                            break;

                        case 'junction':
                            if (consistency_inlet(node) < 1) {
                                alerts.push(`Junction ${name} does not have an inlet`);
                            }
                            if (consistency_outlet(node) < 1) {
                                console.log(name);
                                alerts.push(`Junction ${name} does not have an outlet`);
                            }
                            let inlets = getInletVessels(node);
                            let outlets = getOutletVessels(node);

                            inlets.forEach((inlet, index) => {
                                if (!inletCounts[inlet]) {
                                    inletCounts[inlet] = 0;
                                }
                                inletCounts[inlet]++;
                                console.log(`Inlet [${index}]: ${inlet}, Count: ${inletCounts[inlet]}`);
                            });

                            outlets.forEach((outlet, index) => {
                                if (!outletCounts[outlet]) {
                                    outletCounts[outlet] = 0;
                                }
                                outletCounts[outlet]++;
                                console.log(`Outlet [${index}]: ${outlet}, Count: ${outletCounts[outlet]}`);
                            });

                            if (inlets.length === 1 && outlets.length === 1) {
                                oneToOneJunctions[outlets[0]] = name; // Store junction by outlet
                            }
                            junctionInletsMap[name] = inlets;
                            break;
                        default:
                            break;
                    }
        });

        console.log(oneToOneJunctions);

        // Initialize a dictionary to store junctions and their inlets that are 1:1 outlets
        let junctionToInletsFromOneToOne = {};

        // Populate the dictionary with all junctions receiving inlets that are 1:1 outlets
        Object.keys(junctionInletsMap).forEach(junctionId => {
            let inlets = junctionInletsMap[junctionId];
            let oneToOneInlets = inlets.filter(inlet => oneToOneJunctions.hasOwnProperty(inlet));

            // Store these 1:1 inlets for each junction if there are any
            if (oneToOneInlets.length > 0) {
                junctionToInletsFromOneToOne[junctionId] = oneToOneInlets;
            }
        });

        // Check for any junction that has multiple inlets originating from 1:1 junction outlets
        Object.keys(junctionToInletsFromOneToOne).forEach(junctionId => {
            let inlets = junctionToInletsFromOneToOne[junctionId];
            if (inlets.length > 1) { // More than one 1:1 junction outlet merging into the same junction
                alerts.push(`Junction ${junctionId} has multiple 1:1 junction outlets merging into it.`);
            }
        });

        // Check the counts to make sure a vessel is inlet to a junction max once.
        for (let vessel in inletCounts) {
            if (inletCounts[vessel] > 1) {
                console.log(`Checking vessel: ${vessel}, Count: ${inletCounts[vessel]}`);
                alerts.push(`Vessel ${vessel} is an inlet to more than one junction`);
            }
        }

        // Check the counts to make sure a vessel is outlet to a junction max once.
        for (let vessel in outletCounts) {
            if (outletCounts[vessel] > 1) {
                alerts.push(`Vessel ${vessel} is an outlet to more than one junction`);
            }
        }

        if (data['boundary_conditions'].length < 2) {
            alerts.push('The model needs at least two boundary conditions');
        }

        if (alerts.length == 0) {
            return 1;
        }

        if (alerts.length !== 0) {
            alert(alerts.join('\n'));
            return 0;
        }
    }

    window.Consistency_Check = Consistency_Check;

    // Event listener for the Export to JSON button
    document.getElementById('export-json').addEventListener('click', function() {
        let graphData = exportGraphData();
        if (Consistency_Check(graphData) == 1) {
            downloadJSON(graphData);
        }
    });

    let currentNode = null;

    // Sets styling for the simulation parameters form.
    function setSimParameters() {
        document.getElementById('nodeInfoModal').style.display = 'block';
        document.getElementById('SimParametersForm').style.display = 'block';
        document.getElementById('vesselForm').style.display = 'none';
        document.getElementById('junctionForm').style.display = 'none';
    }

    // Parses the results from the simulation parameters form once the user submits it.
    function submitSimParameters() {
        simulation_parameters_dict.number_of_cardiac_cycles = parseInt(document.getElementById('numcycles').value, 10);
        simulation_parameters_dict.number_of_time_pts_per_cardiac_cycle = parseInt(document.getElementById('numtimepts').value,10);
        simulation_parameters_dict.output_all_cycles = document.getElementById('output_all_cycles').value.toLowerCase() === 'true';
        simulation_parameters_dict.output_variable_based = document.getElementById('output_variable_based').value.toLowerCase() === 'true';
        document.getElementById('nodeInfoModal').style.display = 'none';
        hideNodeInfoModal()
    }


    // Function to show the modal form
    function showNodeInfoModal(nodeType) {
        document.getElementById('nodeInfoModal').style.display = 'block';
        document.getElementById('SimParametersForm').style.display = 'none';

        // Show specific form based on node type
        if (nodeType === 'vessel') {
            document.getElementById('vesselForm').style.display = 'block';
            document.getElementById('junctionForm').style.display = 'none';
        }
        else if (nodeType === 'junction') {
            document.getElementById('junctionForm').style.display = 'block';
            document.getElementById('vesselForm').style.display = 'none';
        }
    }

    // Show or hide the junction parameters form based on the selected junction type
    document.getElementById('junction-type').addEventListener('change', function() {
        if (this.value === 'BloodVesselJunction') {
            document.getElementById('Junction_Parameters_Form').style.display = 'block';
        }
    });

    // Function to hide the modal form
    function hideNodeInfoModal() {
        document.getElementById('nodeInfoModal').style.display = 'none';
        document.getElementById('vesselForm').style.display = 'none'; // Hide the form
    }

    // Show or hide the custom Young's modulus input field based on the selected value
    document.getElementById('young_mod_input').addEventListener('change', function() {
        if (this.value === 'custom') {
            document.getElementById('custom_young_mod').style.display = 'inline';
        } else {
            document.getElementById('custom_young_mod').style.display = 'none';
        }
    });

    // Show or hide the custom viscosity (mu) input field based on the selected value
    document.getElementById('mu_input').addEventListener('change', function() {
        if (this.value === 'custom') {
            document.getElementById('custom_mu').style.display = 'inline';
        } else {
            document.getElementById('custom_mu').style.display = 'none';
        }
    });

    // Show or hide the custom wall thickness (h) input field based on the selected value
    document.getElementById('h_input').addEventListener('change', function() {
        if (this.value === 'custom') {
            document.getElementById('custom_thickness').style.display = 'inline';
        } else {
            document.getElementById('custom_thickness').style.display = 'none';
        }
    });

    // Show or hide the custom density (rho) input field based on the selected value
    document.getElementById('rho_input').addEventListener('change', function() {
        if (this.value === 'custom') {
            document.getElementById('custom_rho').style.display = 'inline';
        } else {
            document.getElementById('custom_rho').style.display = 'none';
        }
    });

    // Parses the information the user inputted for a junction node.
    function submitJunctionInfo() {
        let additionalData = {};
        let junction_type = document.getElementById('junction-type').value;
        if (junction_type === 'BloodVesselJunction') {
            additionalData.junction_type = junction_type;
            additionalData.R1 = parseFloat(document.getElementById('R_poiseuille1').value);
            additionalData.R2 = parseFloat(document.getElementById('R_poiseuille2').value);
            additionalData.L1 = parseFloat(document.getElementById('L1').value);
            additionalData.L2 = parseFloat(document.getElementById('L2').value);
            additionalData.junction_stenosis_coefficient1 = parseFloat(document.getElementById('junction_stenosis_coefficient1').value);
            additionalData.junction_stenosis_coefficient2 = parseFloat(document.getElementById('junction_stenosis_coefficient2').value);
        }
        else {
            additionalData.junction_type = junction_type;
        }
        currentNode.data('additional_data', additionalData); // Update node with additional data
        hideNodeInfoModal(); // Hide the modal form

    }

    // Function to submit node info and update the node
    function submitVesselInfo() {
        let additionalData = {};

        if (currentType === 'vessel') {
            // add the form information to additionalData so it can be attached to the vessel
            additionalData.vessel_length = parseFloat(document.getElementById('vesselLengthInput').value);
            additionalData.vessel_radius = parseFloat(document.getElementById('vesselRadiusInput').value);
            additionalData.stenosis_diameter = parseFloat(document.getElementById('vesselStenosisDiameterInput').value);
            let youngModulus = document.getElementById('young_mod_input').value;
            if (youngModulus === 'custom') {
                youngModulus = parseFloat(document.getElementById('custom_young_mod').value);
            } else {
                youngModulus = parseFloat(youngModulus);
            }
            additionalData.young_mod = youngModulus;

            // Handle dynamic viscosity
            let mu = document.getElementById('mu_input').value;
            if (mu === 'custom') {
                mu = parseFloat(document.getElementById('custom_mu').value);
            } else {
                mu = parseFloat(mu);
            }
            additionalData.mu = mu;

            // Handle vessel thickness
            let h = document.getElementById('h_input').value;
            if (h === 'custom') {
                h = parseFloat(document.getElementById('custom_thickness').value);
            } else {
                h = parseFloat(h);
            }
            additionalData.h = h;

            // Handle density
            let rho = document.getElementById('rho_input').value;
            if (rho === 'custom') {
                rho = parseFloat(document.getElementById('custom_rho').value);
            } else {
                rho = parseFloat(rho);
            }
            additionalData.rho = rho;

            currentNode.data('additional_data', additionalData); // Update node with additional data

            // reset form inputs
            document.getElementById('vesselLengthInput').value = '';
            document.getElementById('vesselStenosisDiameterInput').value = '';
            document.getElementById('vesselRadiusInput').value = '';
        }
        hideNodeInfoModal(); // Hide the modal form
    }

    // Toggles the text on the delete-mode-button.
    function toggleDeleteMode() {
            deleteMode = !deleteMode;
            if (deleteMode) {
                document.getElementById('delete-mode-button').innerText = 'Exit Delete Mode';
            } else {
                document.getElementById('delete-mode-button').innerText = 'Enter Delete Mode';
            }
        }

    // Add event listener to the delete button
    document.getElementById('delete-mode-button').addEventListener('click', toggleDeleteMode);


    // Event listener to show the form when a node is clicked
    cy.on('tap', 'node', function(event) {
        currentNode = event.target;
        currentType = currentNode.data('type');
        if (deleteMode) {
                var node = event.target;
                node.remove();
        }
        else if (currentType == 'vessel') {
            showNodeInfoModal(currentType); // Show the modal form
        }
        else if (currentType == 'junction') {
            showNodeInfoModal(currentType); // Show the modal form
        }

    });

    // Get all elements with the class 'collapsible'
    var coll = document.getElementsByClassName("collapsible");
    var i;

    // Loop through all collapsible elements and add a click event listener to toggle visibility of their associated content
    for (i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = this.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
      });
    }
});
