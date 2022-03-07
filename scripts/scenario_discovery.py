# Script to run data analysis for scenario discovery analysis

import json
import os
import itertools
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from datetime import datetime as dt
import matplotlib.pyplot as plt

import ema_workbench
import sys
sys.path.append(".\\sample")
from SALibRepastParams import num_levels, params, random_seed, init_problem, calc_second_order
from sobol_plot import plot_sobol_indices, save_second_order_sobol_indices

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

#####################################
#
#
# Choose data to analyze and sensitivity analysis setting
#
#
#####################################
file_datetime_string = config['file_datetime_string']
vehicle_density_timestamp = config['vehicle_density_timestamp']
setting = config['setting']

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"

data_paths = bd_utils.get_data_paths(file_datetime_string, data_dir)
ped_crossings_file = data_paths["ped_crossings_file"]
ped_routes_file = data_paths["ped_routes_file"]
cross_events_file = data_paths["cross_events_file"]
ped_locations_file = data_paths["ped_locations_file"]
vehicle_rls_file = data_paths["vehicle_rls_file"]
batch_file = data_paths["batch_file"] 

output_paths = bd_utils.get_ouput_paths(file_datetime_string, data_dir)
output_ped_routes_file=             output_paths["output_ped_routes_file"]
output_route_length_file=           output_paths["output_route_length_file"]
output_sp_similarity_length_path=   output_paths["output_sp_similarity_length_path"]
output_cross_events_path=           output_paths["output_cross_events_path"]


#####################################
#
#
# Load Data
#
#
#####################################

# Data from model run
dfRun = pd.read_csv(os.path.join(data_dir, batch_file)).sort_values(by = 'run')
dfRun['minCrossing'] = dfRun['minCrossing'].astype(int)

# GIS Data
gdfPaveLinks = gpd.read_file(pavement_links_file)
gdfPaveNodes = gpd.read_file(pavement_nodes_file)
gdfORLinks = gpd.read_file(or_links_file)
gdfORNodes = gpd.read_file(or_nodes_file)
gdfITNLinks = gpd.read_file(itn_links_file)
gdfCAs = gpd.read_file(crossing_alternatives_file)

# Get networkx graph and node positions
pavement_graph = nx.Graph()
gdfPaveLinks['length'] = gdfPaveLinks['geometry'].length.map(lambda x: int(x)).replace(0,1) # Repalce 0 with 1 to prevent links with 0 weight. Only affects 4 links. (gdfPaveLinks['geometry'].length < 1).value_counts()
gdfPaveLinks['edge_data'] = gdfPaveLinks.apply(lambda row: {'length':row['length'], 'fid':row['fid']}, axis=1)
edges = gdfPaveLinks.loc[:, ['MNodeFID', 'PNodeFID', 'edge_data']].values
pavement_graph.add_edges_from(edges)

# Using the geographical coordinates of the nodes when plotting them
points_pos = gdfPaveNodes.set_index('fid')
points_pos['x'] = points_pos['geometry'].map(lambda g: g.coords[0][0])
points_pos['y'] = points_pos['geometry'].map(lambda g: g.coords[0][1])
node_posistions = list(zip(points_pos['x'], points_pos['y']))
dict_node_pos = dict(zip(points_pos.index, node_posistions))

######################################
#
#
# Process data to get ped routes and shortest path routes
#
#
######################################
# Not currently using vehicle density in the calculation of route lenght so just set to 0 everywhere
gdfPaveLinks['AvVehDen'] = 0.0
weight_params = range(0, 100, 100)

dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = ped_routes_file)
dfCrossEvents = bd_utils.load_and_clean_cross_events(gdfPaveLinks, cross_events_path = cross_events_file)

# Data aggregated to run level, used to calculate sensitivity indices
dfRouteCompletion = bd_utils.agg_route_completions(dfPedRoutes, dfRun, output_path = output_route_completion_path)

# Important for sensitivity analysis to compare the same set of trips between runs.
# To do this need to exclude peds ID that get removed from all other runs
dfPedRoutesConsistentPeds = dfPedRoutes.loc[ ~dfPedRoutes['ID'].isin(dfPedRoutes_removedpeds['ID'])]
dfCrossEventsConsistentPeds = dfCrossEvents.loc[ ~dfCrossEvents['ID'].isin(dfPedRoutes_removedpeds['ID'])]

dfLinkCrossCounts = bd_utils.get_road_link_pedestrian_crossing_counts(dfCrossEventsConsistentPeds, gdfPaveLinks)


print("\nCalculating/Loading Output Metrics")
dfRouteLength = bd_utils.get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, output_path = output_route_length_file)
dfSPSim = bd_utils.get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'dice_dist', output_path = output_sp_similarity_path)
dfSPSimLen = bd_utils.get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'path_length', output_path = output_sp_similarity_length_path)
dfConflicts = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsMarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unsignalised'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsUnmarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unmarked'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDirect = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='direct_cross'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonal = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='diag_cross'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonalUm = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ (dfCrossEventsConsistentPeds['linkType']=='diag_cross') & (dfCrossEventsConsistentPeds['CrossingType']=='unmarked')], dfLinkCrossCounts, ttc_col = 'TTC')

conflicts_data = {'all':dfConflicts, 'marked':dfConflictsMarked, 'unmarked':dfConflictsUnmarked, 'direct':dfConflictsDirect, 'diag':dfConflictsDiagonal, 'diag_um':dfConflictsDiagonalUm}

