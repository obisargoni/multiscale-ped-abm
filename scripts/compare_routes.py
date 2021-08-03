# Script to analyse which road links pedestrian agents cross on

import json
import os
import itertools
import pandas as pd
import numpy as np
import geopandas as gpd
import re
import networkx as nx
from datetime import datetime as dt
from shapely.geometry import Point
import similaritymeasures as sim

import batch_data_utils as bd_utils

#####################################
#
# Globals
#
#####################################
project_crs = {'init': 'epsg:27700'}

with open(".//config.json") as f:
    config = json.load(f)

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])


# Model output data
file_datetime_string = "2021.Jul.22.17_31_50"
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")
file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime)
ped_crossings_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime, suffix = 'batch_param_map')
batch_file = bd_utils.most_recent_directory_file(data_dir, file_re)


#####################################
#
#
# Load Data
#
#
#####################################

# Data from model run
dfPedCrossings = pd.read_csv(ped_crossings_file)
dfRun = pd.read_csv(os.path.join(data_dir, batch_file))

# GIS Data
gdfPaveLinks = gpd.read_file(pavement_links_file)
gdfPaveNodes = gpd.read_file(pavement_nodes_file)
gdfORLinks = gpd.read_file(or_links_file)
gdfORNodes = gpd.read_file(or_nodes_file)

# Get networkx graph and node positions
pavement_graph = nx.Graph()
gdfPaveLinks['length'] = gdfPaveLinks['geometry'].length.map(lambda x: int(x))
gdfPaveLinks['edge_data'] = gdfPaveLinks.apply(lambda row: {'length':row['length'], 'fid':row['fid']}, axis=1)
edges = gdfPaveLinks.loc[:, ['MNodeFID', 'PNodeFID', 'edge_data']].values
pavement_graph.add_edges_from(edges)

# Using the geographical coordinates of the nodes when plotting them
points_pos = gdfPaveNodes.set_index('fid')
points_pos['x'] = points_pos['geometry'].map(lambda g: g.coords[0][0])
points_pos['y'] = points_pos['geometry'].map(lambda g: g.coords[0][1])
node_posistions = list(zip(points_pos['x'], points_pos['y']))
dict_node_pos = dict(zip(points_pos.index, node_posistions))


#################################
#
#
# Functions
#
#
#################################
def adjacent_edge_node(n, e):
    if n == e[0]:
        return e[1]
    elif n == e[1]:
        return e[0]
    else:
        return None

def node_path_from_edge_path(edge_id_path, start_node, end_node, pavement_graph):
    '''Get a node path from an edge path using the graph the path is on
    '''

    edge_path = []
    for e_id in edge_id_path:
        for e in pavement_graph.edges(data=True):
            if e[-1]['fid'] == e_id:
                edge_path.append(e[:2])

    # Work way backwards through edges to reconstruct path
    node_path = [end_node]
    prev = end_node
    edge_path_reverse = edge_path[::-1]
    for i, e in enumerate(edge_path_reverse):
        candidate = adjacent_edge_node(prev, e)

        # Decide whether this node was part of the path by looking at next edge
        if i >= len(edge_path_reverse)-1:
            # In this case can't check ahead for what the next candidate will be so just set to the current candidate
            next_candidate = candidate 
        else:
            # Now check that candidate connects to next edge
            next_candidate = adjacent_edge_node(candidate, edge_path_reverse[i+1])

        if (next_candidate is not None) & (candidate not in node_path):
             node_path.append(candidate)
             prev = candidate

        if prev == start_node:
            break

    # If ped agent got removed from the simulation it would not have reached the end node and the path will not have any other nodes added
    # In this case return null
    if len(node_path)==1:
        return None

    return node_path[::-1]

def get_strategic_path_pavement_edges(strategic_path, gdfORLinks, gdfPaveNodes):
    filter_pavement_edges = []
    for or_link in strategic_path:
        or_nodes = gdfORLinks.loc[ gdfORLinks['fid'] == or_link, ['MNodeFID', 'PNodeFID']].values[0]

        pavement_nodes = gdfPaveNodes.loc[ gdfPaveNodes['juncNodeID'].isin(or_nodes), 'fid'].values

        # Now get edges
        for u, v in itertools.product(pavement_nodes, pavement_nodes):
            try:
                e_id = pavement_graph[u][v]['fid']
                filter_pavement_edges.append((u,v))
            except KeyError:
                pass
    return filter_pavement_edges

def shortest_path_within_strategic_path(strategic_path, gdfORLinks, gdfPaveNodes, pavement_graph, start_node, end_node, weight = 'length'):

    filter_pavement_edges = get_strategic_path_pavement_edges(strategic_path, gdfORLinks, gdfPaveNodes)

    sub_pavement_graph = nx.edge_subgraph(pavement_graph, filter_pavement_edges)

    try:
        dijkstra_path = nx.dijkstra_path(sub_pavement_graph, start_node, end_node, weight = weight)
    except nx.exception.NodeNotFound as e:
        print(start_node, end_node, strategic_path)
        return None

    return dijkstra_path

def compare_node_paths(npa, npb, dict_node_pos, distance_function = sim.frechet_dist):
    pos_a = [dict_node_pos[i] for i in npa]
    pos_b = [dict_node_pos[i] for i in npb]

    d = distance_function(pos_a, pos_b)
    return d


######################################
#
#
# Extract pedestrian pavement node routes
#
#
######################################

# Function to get edge route from data
dfPedRoutes = dfPedCrossings.groupby(['run', 'ID'])['CurrentPavementLinkID'].apply(lambda s: list(s.unique())).reset_index()
dfPedRoutes.rename(columns = {'CurrentPavementLinkID':'edge_path'}, inplace=True)

# Get strategic path link list. First check that there is a single strategic path per run-ped
assert (dfPedCrossings.groupby(['run','ID'])['FullStrategicPathString'].apply(lambda s: s.drop_duplicates().shape[0]).value_counts().index == 1).all()
dfPedStratPaths = dfPedCrossings.groupby(['run','ID'])['FullStrategicPathString'].apply(lambda s: s.drop_duplicates().values[0].split(':')[1:]).reset_index()
dfPedRoutes = pd.merge(dfPedRoutes, dfPedStratPaths, on = ['run', 'ID'])

dfPedOD = dfPedCrossings.groupby(['run', 'ID'])['StartPavementJunctionID', 'DestPavementJunctionID'].apply(lambda df: df.drop_duplicates()).reset_index()
dfPedRoutes = pd.merge(dfPedRoutes, dfPedOD, left_on = ['run', 'ID'], right_on = ['run', 'ID'])


dfPedRoutes['node_path'] = dfPedRoutes.apply(lambda row: node_path_from_edge_path(row['edge_path'], row['StartPavementJunctionID'], row['DestPavementJunctionID'], pavement_graph), axis=1)

dfPedRoutes = dfPedRoutes.loc[ ~dfPedRoutes['node_path'].isnull()]

######################################
#
#
# Calculate shortest distance based paths
#
#
######################################
dfPedRoutes['dist_sp'] = dfPedRoutes.apply(lambda row: nx.dijkstra_path(pavement_graph, row['StartPavementJunctionID'], row['DestPavementJunctionID'], weight = 'length'), axis=1)

# Calculating shortest path from start ot end node can produce a path that travels along a different set of OR road links
# For a more constrained comparison need to limit the pavement network to just the edges along the startegic path
dfPedRoutes['dist_sp_filter'] = dfPedRoutes.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['StartPavementJunctionID'], row['DestPavementJunctionID'], weight = 'length'), axis=1)

######################################
#
#
# Compare paths using a similarity metric
#
#
######################################
dfPedRoutes['frechet_distance'] = dfPedRoutes.apply(lambda row: compare_node_paths(row['node_path'], row['dist_sp_filter'], dict_node_pos, distance_function = sim.frechet_dist), axis=1)

dfFrechet = dfPedRoutes.groupby('run')['frechet_distance'].describe()
dfFrechet = pd.merge(dfFrechet, dfRun, left_index = True, right_on = 'run')
dfFrechet.loc[:, ['tacticalPlanHorizon', 'mean','50%','max']].sort_values(by = 'tacticalPlanHorizon')


# In some cases the first pavement edge does not get included in the data. I think due to the ped origin being v close to pavement junction
# In this case need to alter what the start node is
# So if node path doesn't include start node, calculate alternative path using first node in path.
