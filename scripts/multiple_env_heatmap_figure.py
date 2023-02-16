# Figure showing routes that pedestrians took acrossing all simulation runs in each environment
import json
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import networkx as nx
from matplotlib import pyplot as plt

import batch_data_utils as bd_utils


#####################################
#
#
# Load config file
#
#
#####################################
with open(".//config.json") as f:
    config = json.load(f)


#################################
#
#
# Functions
#
#
#################################
def load_environment_data(env_config, env_timestamp, data_dir = "..\\output\\batch\\model_run_data\\"):
    gis_data_dir = os.path.join(env_config['gis_data_dir'], 'processed_gis_data')
    pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
    pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
    or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
    or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])
    itn_links_file = os.path.join(gis_data_dir, config['mastermap_itn_processed_direction_file'])
    itn_nodes_file = os.path.join(gis_data_dir, config['mastermap_node_processed_file'])
    crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])
    ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])

    gdfORLinks = gpd.read_file(or_links_file)
    gdfORNodes = gpd.read_file(or_nodes_file)
    gdfPaveNodes = gpd.read_file(pavement_nodes_file)
    gdfPedODs = gpd.read_file(ped_ods_file)

    ped_routes_file = os.path.join(data_dir, "pedestrian_routes.{}.csv".format(env_timestamp))

    weight_params = range(0, 100, 100)
    dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(None, None, None, None, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)

    # Get start and end nodes
    gdfStartNodes = gdfPaveNodes.loc[ gdfPaveNodes['fid'].isin(dfPedRoutes['start_node'])]
    gdfEndNodes = gdfPedODs.loc[gdfPedODs['fid'] =='od_0']

    return [gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes, gdfStartNodes, gdfEndNodes]

def paths_heatmap_data(gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes):
    or_graph = nx.Graph()
    gdfORLinks['length'] = gdfORLinks['geometry'].length.map(lambda x: int(x)).replace(0,1)
    gdfORLinks['edge_data'] = gdfORLinks.apply(lambda row: {'length':row['length'], 'fid':row['fid']}, axis=1)
    edges = gdfORLinks.loc[:, ['MNodeFID', 'PNodeFID', 'edge_data']].values
    or_graph.add_edges_from(edges)

    # Using the geographical coordinates of the nodes when plotting them
    points_pos = gdfORNodes.set_index('node_fid')
    points_pos['x'] = points_pos['geometry'].map(lambda g: g.coords[0][0])
    points_pos['y'] = points_pos['geometry'].map(lambda g: g.coords[0][1])
    node_posistions = list(zip(points_pos['x'], points_pos['y']))
    dict_node_pos = dict(zip(points_pos.index, node_posistions))

    # Get count of traversed edges across all runs
    dfPathORLinks = bd_utils.explode_data(dfPedRoutes.reindex(columns = ['run','ID','FullStrategicPathString']), explode_col='FullStrategicPathString')
    n_paths = dfPedRoutes['run'].unique().shape[0] * dfPedRoutes['ID'].unique().shape[0]
    edge_traverse_counts = dfPathORLinks['FullStrategicPathString'].value_counts()

    path_links = dfPathORLinks['FullStrategicPathString'].unique()
    edgedata = np.array([edge_traverse_counts[i] / n_paths for i in path_links])

    #tp_links = edge_traverse_counts.index()

    edgelist = []
    for edge_id in path_links:
        for e in list(or_graph.edges(data=True)):
            if edge_id == e[-1]['fid']:
                edgelist.append(e)

    return or_graph, dict_node_pos, edgelist, edgedata

def multi_paths_heatmap_plot(environment_data_dict, title = "Proportion of pedestrians traversing each road link", title_size = 24, cmap = 'Reds', labelsize=15, cbar_pad=0.03, label_pad = 20, figsize=(20,20)):
    fig , axs = plt.subplots(1,3,figsize=figsize)

    #  Varying density along a streamline
    for i, (k, v) in enumerate(environment_data_dict.items()):
        or_graph, dict_node_pos, edgelist, edgedata, gdfORLinks, gdfORNodes, gdfPaveNodes, gdfStartNodes, gdfEndNodes = v # unpack data from dictionary value

        vmin, vmax = pd.Series(edgedata).quantile([0.0, 0.99])

        print("N values greater that vmax:{}".format( (pd.Series(edgedata)>vmax).value_counts()[True]))

        ax = bd_utils.figure_rl_paths_heatmap(fig, axs[i], gdfORLinks, gdfStartNodes, gdfEndNodes, or_graph, dict_node_pos, edgelist, edgedata, plt.get_cmap(cmap), None, None, None, labelsize, fig_config, vlims = [vmin, vmax], cbar_pad = cbar_pad)

    outpath = os.path.join(img_dir, 'paths_heatmap.png')
    fig.suptitle(title, size=title_size)
    fig.savefig(outpath)

    return outpath
    
#################################
#
#
# Globals
#
#
#################################
plt.style.use('dark_background')
plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

with open("figure_config.json") as f:
    fig_config = json.load(f)



data_dir = "..\\output\\batch\\model_run_data\\"
img_dir = "..\\output\\img\\"

toygrid169_timestamp = "2022.Dec.11.12_56_58"
quadgrid100_timestamp = "2022.Dec.13.01_40_30"
cc_timestamp = "2022.Dec.19.16_42_32"

toygrid_config_path = "C:\\Users\\Obi Sargoni\\Documents\\CASA\\OSPavements\\configs\\config_toygrid_block_141nodes_buffered.json"
quadgrid_config_path = "C:\\Users\\Obi Sargoni\\Documents\\CASA\\OSPavements\\configs\\config_quadgrid_block_145nodes_buffer.json"
cc_config_path = "C:\\Users\\Obi Sargoni\\Documents\\CASA\\OSPavements\\configs\\config_claphamcommon.json"

environment_details =   {  'ToyGrid169':[toygrid_config_path, toygrid169_timestamp],
                            'QuadGrid100':[quadgrid_config_path, quadgrid100_timestamp],
                            'ClaphamCommon':[cc_config_path, cc_timestamp]}

environment_data_dict = {}


#####################################
#
#
# Load GIS and model output data
#
#
#####################################

for env_name, (config_path, env_timestamp) in environment_details.items():

    # Load data for this environment and store in dictionary
    with open(config_path) as f:
        config = json.load(f)
        gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes, gdfStartNodes, gdfEndNodes = load_environment_data(config, env_timestamp)
        or_graph, dict_node_pos, edgelist, edgedata = paths_heatmap_data(gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes)
        environment_data_dict[env_name] = [or_graph, dict_node_pos, edgelist, edgedata, gdfORLinks, gdfORNodes, gdfPaveNodes, gdfStartNodes, gdfEndNodes]
        dfPedRoutes=None

output_road_network_fig_path = os.path.join(img_dir, "paths_heatmap.png")


# Combine into a single dictionary

#################################
#
#
# Get section of data to plot
#
#
#################################
outpath = multi_paths_heatmap_plot(environment_data_dict, title = "Proportion of pedestrians agents traversing each road link", title_size=24, cmap = 'plasma', labelsize=15, cbar_pad=0.03, label_pad = 20, figsize=(30,12))
