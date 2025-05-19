# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

import sys
import argparse
import pysvzerod
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import json
import dash
from dash import html, dcc
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import os
from dirgraph_utils import set_up_0d_network


'''
This file enables the visualization of 0D simulation results in a web app that displays the 0D network 
as a directed graph. Users can interactively select nodes to view their parameters and simulation results.

Enter the filepath for your simulation input JSON file and the directory where 
you want to save the output directed graph as command line arguments to run the script. 
If you want to save the raw svZeroDSolver simulation results, add "export-csv" as the third command line argument.

'''

def dirgraph(filepath, output_dir, export_csv):
    solver = pysvzerod.Solver(filepath)
    solver.run()
    results = pd.DataFrame(solver.get_full_result())

    if export_csv:
        results.to_csv('results.csv', sep=',', index=False, encoding='utf-8')
        print(f"Results exported to results.csv")

    with open(filepath, 'r') as infile:
        parameters = json.load(infile)

    set_up_0d_network(
        filepath,
        name_type='id', # Options 'name' or 'id' specifies whether vessel names or ids should be used for each node
        draw_directed_graph= False,  # Enter True if you want to save the directed graph
        output_dir= output_dir
    )
    base_name = filepath.rsplit('/', 1)[-1]
    output_file = os.path.join(output_dir, os.path.splitext(base_name)[0] + "_directed_graph.dot")
    G = nx.DiGraph(nx.nx_pydot.read_dot(output_file))
    return results, parameters, G


def main():
    parser = argparse.ArgumentParser(
        description="Generate a directed graph visualization of your 0d network from the JSON input file."
    )

    # Positional arguments
    parser.add_argument(
        'filepath',
        type=str,
        help="Path to the svZeroDSolver input JSON file."
    )

    parser.add_argument(
        'output_dir',
        type=str,
        help="Directory to store the visualization results."
    )

    # Optional argument
    parser.add_argument(
        '--export_csv',
        action='store_true',
        help="If specified, export the results as CSV files."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the dirgraph function with the provided arguments
    return dirgraph(args.filepath, args.output_dir, args.export_csv)

if __name__ == '__main__':
    results, parameters, G = main()

# Mapping vessel names to IDs. The block parameters are obtained from the user's input json file. 
vessel_name_to_vessel_id_map = {}
vessel_id_to_vessel_name_map = {}
vessel_params = {}
for vessel in parameters["vessels"]:
    vessel_name = vessel["vessel_name"]
    vessel_id = "V" + str(vessel["vessel_id"])
    vessel_name_to_vessel_id_map[vessel_name] = vessel_id
    vessel_id_to_vessel_name_map[vessel_id] = vessel_name

    inlet = None
    outlet = None
    if 'boundary_conditions' in vessel:
        if 'inlet' in vessel['boundary_conditions']:
            inlet = vessel['boundary_conditions']['inlet']
        if 'outlet' in vessel['boundary_conditions']:
            outlet = vessel['boundary_conditions']['outlet']

    vessel_length = vessel["vessel_length"]
    zero_d_element_values = vessel["zero_d_element_values"]
    vessel_params[vessel_id] = {'vessel_length': vessel_length,
                "zero_d_element_values": zero_d_element_values, 'inlet': inlet, 'outlet': outlet}

junction_info = {}
if 'junctions' in parameters:
    for junction in parameters['junctions']:
        junction_info[junction['junction_name']] = junction['junction_type']

valve_info = {}
if 'valves' in parameters:
    for valve in parameters['valves']:
        valve_info[valve['name']] = valve['params']

chamber_info = {}
if 'chambers' in parameters:
    for chamber in parameters['chambers']:
        chamber_info[chamber['name']] = chamber['values']

bc_info = {}
if 'boundary_conditions' in parameters:
    for bc in parameters['boundary_conditions']:
        bc_info[bc['bc_name']] = bc['bc_type']


pos = nx.nx_pydot.pydot_layout(G, prog='dot')

# Prepare Plotly graph using temporary lists
edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])

# Now create the go.Scatter object
edge_trace = go.Scatter(
    x=edge_x,
    y=edge_y,
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

# Prepare node data
vessel_x = []
vessel_y = []
vessel_text = []
junction_x = []
junction_y = []
junction_text = []
bc_x = []
bc_y = []
bc_text = []
chamber_x = []
chamber_y = []
chamber_text = []
valve_x = []
valve_y = []
valve_text = []


bc_nodename_to_realname = {}
bc_realname_to_nodename = {}
for node in G.nodes():
    x, y = pos[node]
    if node.startswith('valve'):
        valve_x.append(x)
        valve_y.append(y)
        valve_text.append(node)
    elif node.startswith('V'):
        vessel_x.append(x)
        vessel_y.append(y)
        vessel_text.append(node)
    elif node.startswith('J'):
        junction_x.append(x)
        junction_y.append(y)
        junction_text.append(node)
    elif node.startswith('BC'):
        bc_x.append(x)
        bc_y.append(y)
        # bc_text.append(node) #This will add BC14_outlet etc that structure
        bc_id = node.split('BC')[1].split('_')[0]
        location = node.split('_')[1]
        if bc_id.isdigit():  # it's a boundary condition to a vessel
            # node_parent = vessel_id_to_vessel_name_map["V" + bc_id]
            node_parent = "V" + str(bc_id)
            real_node_name = vessel_params[node_parent][location]
        else:  # it's a boundary condition to a valve
            if location == 'inlet':
                location = 'upstream_block'
            else:
                location = 'downstream_block'
            real_node_name = valve_info[bc_id][location]
        bc_nodename_to_realname[node] = real_node_name
        bc_realname_to_nodename[real_node_name] = node
        bc_text.append(real_node_name)
    else:
        chamber_x.append(x)
        chamber_y.append(y)
        chamber_text.append(node)


# Function to normalize names
def normalize_names(name):
    if name in vessel_name_to_vessel_id_map:
        return vessel_name_to_vessel_id_map[name]
    elif any(bc['bc_name'] == name for bc in parameters['boundary_conditions']):
        return bc_realname_to_nodename[name]
    return name

if "pressure_in" in results.columns:
    # Apply normalization to the entire 'name' column efficiently
    results['name'] = results['name'].apply(normalize_names)
    grouped_data = results.groupby('name')

else:
    # Split the name column into data_type, structure, and position
    results[['data_type', 'structure', 'position']] = results['name'].str.split(':', expand=True)

    # Initialize a list to hold structured data
    structured_data_list = []

    # Group by structure (vessel)
    for (struct, time), group in results.groupby(['structure', 'time']):
        # Create a dictionary to store the values for each vessel
        struct_data = {
            'name': struct,
            'time': time,
            'pressure_in': None,
            'pressure_out': None,
            'flow_in': None,
            'flow_out': None,
        }

        # Assign values to the appropriate columns based on type and position
        for _, row in group.iterrows():
            if row['data_type'] == 'pressure':
                if row['position'] == 'INLET':
                    struct_data['pressure_in'] = row['y']
                else:
                    struct_data['pressure_out'] = row['y']
            elif row['data_type'] == 'flow':
                if row['position'] == 'INLET':
                    struct_data['flow_in'] = row['y']
                else:
                    struct_data['flow_out'] = row['y']
            else:
                struct_data[row['data_type']] = row['y']

        # Append the structured data to the list
        structured_data_list.append(struct_data)

    # Convert the list to a DataFrame
    df = pd.DataFrame(structured_data_list)
    df['name'] = df['name'].apply(normalize_names)
    grouped_data = df.groupby('name')

valve_trace = go.Scatter (
    name = 'Valve',
    x=valve_x,
    y=valve_y,
    text= valve_text,
    textposition='top center',
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        symbol='circle',
        size=10,
        color='purple',
        line_width=2
    )
)

chamber_trace = go.Scatter (
    name = 'Chamber',
    x=chamber_x,
    y= chamber_y,
    text= chamber_text,
    textposition='top center',
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        symbol='circle',
        size=10,
        color='orange',
        line_width=2
    )
)

bc_trace = go.Scatter(
    name = 'Boundary Condition',
    x=bc_x,
    y=bc_y,
    text=bc_text,
    textposition='top center',
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        symbol='circle',
        size=10,
        color='green',
        line_width=2
    )
)

vessel_trace = go.Scatter(
    name = 'Vessel',
    x=vessel_x,
    y=vessel_y,
    text=vessel_text,
    textposition='top center',
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        symbol='circle',
        size=10,
        color='blue',
        line_width=2
    )
)

junction_trace = go.Scatter(
    name = 'Junction',
    x=junction_x,
    y=junction_y,
    text=junction_text,
    textposition='top center',
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        symbol='circle',
        size=10,
        color='red',
        line_width=2
    )
)


# Initialize Dash app
app = dash.Dash(__name__)

# Define the layout
app.layout = html.Div([
    # Button to open the modal
    html.Button('Instructions', id='open-modal-button', style = { 'width': '90px',
    'height': '30px'}, n_clicks=0),

    # Modal div
    html.Div([  # modal div
        html.Div([  # content div
            html.Div([
                'This platform allows you to visualize and analyze SimVascular’s ZeroD Solver results. When you click a node, available pressure and flow data will be displayed alongside the simulation parameters for that given node.',
            ]),
            html.Hr(),
            html.Button('Close', id='modal-close-button', n_clicks=0)
        ],
            style={
                'textAlign': 'center',
                'margin': '90px',
                'padding': '30px',
                'backgroundColor': 'white'
            },
            className='modal-content',
        ),
    ],
        id='modal',
        className='modal',
        style={
            'display': 'none',  # Initially hidden
            'position': 'fixed',
            'zIndex': 1002,
            'left': '0',
            'top': '0',
            'width': '100%',
            'height': '100%',
            'backgroundColor': 'rgba(0, 0, 0, 0.6)'
        }
    ),

    dcc.Graph(
        id='network-graph',
        figure={
            'data': [edge_trace, vessel_trace, chamber_trace, valve_trace, junction_trace, bc_trace],
            'layout': go.Layout(
                title='Network Graph',
                titlefont_size=30,
                title_x=0.5,
                showlegend=True,
                hovermode='closest',
                margin=dict(b=20, l=0, r=0, t=50),
                autosize=True,
                # width=1600,  # Set the width of the plot
                # height=950,
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)

            )
        }
        ,
        style={'width': '100%', 'height': '100%'}  # 100% width and 100% viewport height
        ),
    html.H1(id='joint-title', style = {'textAlign': 'center'}, children="Click on a node to see data"),
    html.Div(id='simulation_info', style={
        'backgroundColor': '#e5ecf6',
        'whiteSpace': 'pre-line',
        'width': '90%',
        'margin': '0 auto',
        'fontSize': '20px',
        'textAlign': 'center'
    }),
    html.Div(id='charts-container', style={'display': 'flex', 'flexDirection': 'column'})
], style={'marginLeft': 'auto', 'marginRight': 'auto'})



# Callback to toggle modal visibility
@app.callback(
    Output('modal', 'style'),
    [Input('open-modal-button', 'n_clicks'),
     Input('modal-close-button', 'n_clicks')],
    [State('modal', 'style')]
)
def toggle_modal(open_clicks, close_clicks, current_style):
    if open_clicks > close_clicks:
        return {
            'display': 'block',
            'position': 'fixed',
            'zIndex': 1002,
            'left': '0',
            'top': '0',
            'width': '100%',
            'height': '100%',
            'backgroundColor': 'rgba(0, 0, 0, 0.6)'
        }
    return {
        'display': 'none',
        'position': 'fixed',
        'zIndex': 1002,
        'left': '0',
        'top': '0',
        'width': '100%',
        'height': '100%',
        'backgroundColor': 'rgba(0, 0, 0, 0.6)'
    }

@app.callback(
    [
        Output('joint-title', 'children'),
        Output('simulation_info', 'children'),
        Output('charts-container', 'children'),
        Output('network-graph', 'figure')
    ],
    [Input('network-graph', 'clickData')]
)

# This adds a box around the node that the user selected
def update_graphs(clickData):
    # Initialize empty figures to handle cases with no data or no node clicks
    empty_fig = go.Figure(layout=go.Layout(title="Click on a node to see data"))
    info_text = ''

    # if clickData:
    fig = go.Figure(
        data=[edge_trace, vessel_trace, chamber_trace, valve_trace, junction_trace, bc_trace],
        layout=go.Layout(
            title='Network Graph',
            titlefont_size=30,
            title_x=0.5,
            showlegend=True,
            hovermode='closest',
            margin=dict(b=20, l=0, r=0, t=50),
            autosize=True,
            width=1700,
            height=800,
            # width=1600,
            # height=950,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            shapes=[],
            paper_bgcolor='white',  # Set the paper background color to white
            plot_bgcolor='white'
        )
    )

    if clickData and 'points' in clickData and len(clickData['points']) > 0:
        point = clickData['points'][0]
        x_selected = point['x']
        y_selected = point['y']

        annotation = {
            'xref': 'x',
            'yref': 'y',
            'x': x_selected,
            'y': y_selected,
            'showarrow': False,
            'text': '',
            'font': {'color': 'RoyalBlue'},
            'align': 'center',
            'bordercolor': 'RoyalBlue',
            'borderwidth': 2,
            'borderpad': 4,
            'bgcolor': 'rgba(0,0,0,0)',
            'opacity': 0.8,
            'width': 10,
            'height': 10
        }
        fig.add_annotation(annotation)

        node_name = clickData['points'][0]['text']
        node_real_name = node_name

        if node_name.startswith('J'):
            if junction_info[node_name] == "NORMAL_JUNCTION":
                junction_type = "Normal Junction"
            else:
                junction_type = "Blood Vessel Junction"

            info_text = [html.B("Junction information", style={'fontSize': '22px', "text-decoration": "underline"}), html.Br(),
                         html.B("Junction ID: "), f"{node_name}", html.Br(),
                         html.B("Junction Type: "), f"{junction_type}", html.Br()]
            return f"Data for {node_name}", info_text, [dcc.Graph(figure=empty_fig, style={'flex': '1'})], fig

        if node_name.startswith('BC') and node_name not in grouped_data.groups:
            real_node_name = bc_nodename_to_realname[node_name]
            info_text = [
                html.B("Boundary Condition: "), f"{real_node_name.capitalize()}, ", html.Br(),
                html.B("Boundary Condition Type: "), f"{bc_info[real_node_name].capitalize()}", html.Br()
            ]
            return f"No simulation data for {node_name}", info_text, [dcc.Graph(figure=empty_fig, style={'flex': '1'})], fig

        if node_name in grouped_data.groups:
            if node_name.startswith('V'):
                node_real_name = vessel_id_to_vessel_name_map[node_name]
                vessel_info = vessel_params[node_name]
                zero_d_values = vessel_info['zero_d_element_values']

                info_text = [
                    html.B("Vessel information", style={'fontSize': '22px', "text-decoration": "underline"}), html.Br(),
                    html.B("Vessel ID: "), f"{node_name}", html.Br(),
                    html.B("Vessel Length: "), f"{vessel_info['vessel_length']}", html.Br(),
                    "\n",
                    html.B("Zero D Element Values", style = {'fontSize': '22px', "text-decoration": "underline"} ), html.Br(),
                ]
                if 'C' in zero_d_values:
                    info_text.extend([html.B("C: "), f"{zero_d_values['C']}, ", html.Br()])
                if 'L' in zero_d_values:
                    info_text.extend([html.B("L: "), f"{zero_d_values['L']}, ", html.Br()])
                if 'R_poiseuille' in zero_d_values:
                    info_text.extend([html.B("R poiseuille: "), f"{zero_d_values['R_poiseuille']}, ", html.Br()])
                if 'stenosis_coefficient' in zero_d_values:
                    info_text.extend([html.B("Stenosis Coefficient: "), f"{zero_d_values['stenosis_coefficient']}"])

            # Add valve information
            elif node_name.startswith('valve'):
                valve_params = valve_info[node_name]
                info_text = [
                    html.B("Valve information", style={'fontSize': '22px', "text-decoration": "underline"}), html.Br(),
                    html.B("Valve ID: "), f"{node_name}", html.Br(),
                    "\n",
                    html.B("Valve Paramaters", style={'fontSize': '22px', "text-decoration": "underline"}),
                    html.Br(),
                    html.B("Rmax: "), f"{valve_params['Rmax']}, ", html.Br(),
                    html.B("Rmin: "), f"{valve_params['Rmin']}, ", html.Br(),
                    html.B("Steepness: "), f"{valve_params['Steepness']}, ", html.Br(),
                    html.B("Upstream Block: "), f"{valve_params['upstream_block'].capitalize()}, ", html.Br(),
                    html.B("Downstream Block: "), f"{valve_params['downstream_block'].capitalize()}", html.Br()
                ]
            elif node_name.startswith('BC'):
                real_node_name = bc_nodename_to_realname[node_name]
                info_text = [
                    html.B("Boundary Condition: "), f"{real_node_name.capitalize()}, ", html.Br(),
                    html.B("Boundary Condition Type: "), f"{bc_info[real_node_name].capitalize()}", html.Br()
                ]

            # elif not node_name.startswith('BC'):
            else:
                chamber_params = chamber_info[node_name]
                info_text = [
                    html.B("Chamber information", style={'fontSize': '22px', "text-decoration": "underline"}), html.Br(),
                    html.B("Chamber Name: "), f"{node_name.capitalize()}", html.Br(),
                    "\n",
                    html.B("Chamber Paramaters", style={'fontSize': '22px', "text-decoration": "underline"}),
                    html.Br(),
                    html.B("Emax: "), f"{chamber_params['Emax']}, ", html.Br(),
                    html.B("Emin: "), f"{chamber_params['Emin']}, ", html.Br(),
                    html.B("Vrd: "), f"{chamber_params['Vrd']}, ", html.Br(),
                    html.B("Vrs: "), f"{chamber_params['Vrs']}, ", html.Br(),
                    html.B("t_active: "), f"{chamber_params['t_active']}, ", html.Br(),
                    html.B("t_twitch: "), f"{chamber_params['t_twitch']}, ", html.Br(),
                    html.B("Impedance: "), f"{chamber_params['Impedance']}", html.Br(),
                ]

            group = grouped_data.get_group(node_name)

            structured_data = {}
            charts = []

            # Iterate over columns to handle data directly based on column names
            for column in group.columns:
                if column not in ['name', 'time']:  # Exclude non-data columns
                    # Collect all values in the column for plotting
                    structured_data[column] = group[column].values

            # Generate plots for each data key found

            for data_key, data_values in structured_data.items():
                if pd.isnull(data_values).all():
                    continue
                data_fig = go.Figure()
                normalized_key = data_key.replace('_', ' ').title()
                title = normalized_key + " for " + node_name.capitalize()

                data_fig.add_trace(go.Scatter(
                    x=group['time'],  # Assuming 'time' is your x-axis data
                    y=data_values,
                    mode='lines+markers',
                    name=title # Beautify the legend name
                ))
                # Configure plot layout
                if data_key == 'pressure_in' or data_key == 'pressure_out':
                    units = 'Pressure [dynes/cm^2]'
                elif data_key == 'flow_in' or data_key == 'flow_out':
                    units = 'Flow [mL/s]'
                else:
                    units = ''

                data_fig.update_layout(
                    title={
                        'text': f"{title}",
                        'x': 0.5,
                        'xanchor': 'center',
                        'yanchor': 'top',
                        'font': {'size': 20, 'color': 'black', 'weight': 'bold'}
                    },
                    title_x=0.5,
                    xaxis=dict(title='Time'),
                    yaxis=dict(title=f'{units}')
                )
                # Add plot to charts list
                charts.append(dcc.Graph(figure=data_fig, style={'flex': '1'}))

            joint_title = f'Data for {node_real_name}'
            return joint_title, info_text, charts, fig
    return f"No data for this node. Please select another", info_text, [dcc.Graph(figure=empty_fig, style={'flex': '1'})], fig


if __name__ == '__main__':
    app.run_server(debug=True)
