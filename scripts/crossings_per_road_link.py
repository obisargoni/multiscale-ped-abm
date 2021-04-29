# Script to analyse which road links pedestrian agents cross on

import json
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from datetime import datetime as dt

import batch_data_utils as bd_utils

#####################################
#
# Globals
#
#####################################

with open(".//gis_data_processing//config.json") as f:
    config = json.load(f)

gis_data_dir = os.path.join(config['gis_data_dir'], "processed_gis_data")
data_dir = config['batch_data_dir'] #"..\\output\\batch\\model_run_data\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")
output_directory = "..\\output\\processed_data\\"

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])


# Output paths
img_dir = "..\\output\\img\\"
crossing_network_fig = os.path.join(img_dir, "crossing_network.png")

output_data_dir = "..\\output\\processed_data\\"

project_crs = {'init': 'epsg:27700'}

file_datetime_string = "2020.Feb.20.14_58_33"
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")

# Load batch result data and gis data
file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings")
ped_crossings_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

dfPedCrossings = pd.read_csv(ped_crossings_file)

gdfPaveNetwork = gpd.read_file(pavement_links_file)
gdfORLinks = gpd.read_file(or_links_file)
gdfORNodes = gpd.read_file(or_nodes_file)


# Filter to just include paveent links used for crossing
dfPedCrossings = dfPedCrossings.loc[ ~dfPedCrossings['TraversedPavementLinkID'].isnull()]

dfPedCrossings = dfPedCrossings.merge(gdfPaveNetwork, left_on = 'TraversedPavementLinkID', right_on = 'fid', how = 'left', indicator=True)
assert dfPedCrossings.loc[ dfPedCrossings['_merge']!='both'].shape[0]==0
dfPedCrossings.drop('_merge', axis=1, inplace=True)

# Aggregate
dfCrossingCounts = dfPedCrossings.groupby(['run', 'pedRLID']).apply(lambda g: g.shape[0]).reset_index()
dfCrossingCounts.rename(columns = {0:'cross_count'}, inplace=True)

# Join with pavement Road link data
gdfCrossingCounts = pd.merge(gdfORLinks, dfCrossingCounts, left_on = 'fid', right_on = 'pedRLID', how = 'left')
gdfCrossingCounts['cross_count'] = gdfCrossingCounts['cross_count'].fillna(0)


# Plot
import networkx as nx
from matplotlib import cm # for generating colour maps
from matplotlib import pyplot as plt

# Get networkx graph
G = nx.Graph()
edge_data = gdfCrossingCounts.loc[:, ['MNodeFID', 'PNodeFID', 'fid', 'cross_count']].values
G.add_edges_from(edge_data[:,:2], fid=edge_data[:,2], cross_count = edge_data[:,3])


# Using the geographical coordinates of the nodes when plotting them
points_pos = gdfORNodes.set_index('node_fid')
points_pos['x'] = points_pos['geometry'].map(lambda g: g.coords[0][0])
points_pos['y'] = points_pos['geometry'].map(lambda g: g.coords[0][1])
points_pos['pos'] = list(zip(points_pos['x'], points_pos['y']))

# Get dictionary of node id to (lat:lon)
geo_pos = points_pos['pos'].to_dict()

# Get edge colour map based on number of crossings

# unpack betweenness values for each node into a list
cmap = cm.get_cmap('viridis') # get a preset colourmap

# Normalise betweenness centrality measures so that they lie between 0-1.
# This is for generating colours and not for analysis
max_cc = float(max(edge_data[:,3]))
cc_norm = [i/max_cc for i in edge_data[:,3]]
edge_palette = [cmap(i) for i in cc_norm]
#widths = [i*6 for i in list_edge_betcen_norm]

plt.figure(figsize = (15,15))
nx.draw_networkx_nodes(G, geo_pos, nodelist=G.nodes(), node_color = 'grey', node_size = 1, alpha = 0.5)
nx.draw_networkx_edges(G, geo_pos, edgelist=G.edges(), width = 3, edge_color = edge_palette, alpha=1)
plt.axis('off')
plt.savefig(crossing_network_fig)