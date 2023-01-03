# Script to run data analysis for scenario discovery analysis

import json
import os
import itertools
import pandas as pd
import numpy as np
import geopandas as gpd
import networkx as nx
import re
from datetime import datetime as dt

import sys
sys.path.append(".\\sample")
from SALibRepastParams import num_levels, params, random_seed, init_problem, calc_second_order, policies

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
setting = config['setting']

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
nbins = None
bin_dist = 2
ttc_threshold = config['ttc_threshold']

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])
itn_links_file = os.path.join(gis_data_dir, config['mastermap_itn_processed_direction_file'])
itn_nodes_file = os.path.join(gis_data_dir, config['mastermap_node_processed_file'])
crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])

data_paths = bd_utils.get_data_paths(file_datetime_string, data_dir)
ped_crossings_file = data_paths["pedestrian_pave_link_crossings"]
ped_routes_file = data_paths["pedestrian_routes"]
veh_routes_file = data_paths["vehicle_routes"]
cross_events_file = data_paths["cross_events"]
ped_locations_file = data_paths["pedestrian_locations"]
vehicle_rls_file = data_paths["vehicle_road_links"]
batch_file = data_paths["batch_file"] 

output_paths = bd_utils.get_ouput_paths(file_datetime_string, data_dir, nbins = bin_dist, ttc_threshold=ttc_threshold)
output_ped_routes_file = output_paths["output_ped_routes_file"]
output_route_length_file = output_paths["output_route_length_file"]
output_ped_distdurs_file = output_paths["output_ped_distdurs_file"]
output_veh_distdurs_file = output_paths["output_veh_distdurs_file"]
output_sp_similarity_length_path = output_paths["output_sp_similarity_length_path"]
output_cross_events_path = output_paths["output_cross_events_path"]
output_cross_entropy = output_paths["output_cross_entropy"]
output_link_cross_entropy = output_paths["output_link_cross_entropy"]
output_cross_conflicts = output_paths["output_cross_conflicts"]
output_ped_trip_length = output_paths["output_ped_trip_length"]

output_sd_data = output_paths["output_sd_data"]


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
gdfPedODs = gpd.read_file(ped_ods_file)

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

dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)
dfCrossEvents = bd_utils.load_and_clean_cross_events(gdfPaveLinks, cross_events_path = cross_events_file)
#dfRouteCompletion = bd_utils.agg_route_completions(dfPedRoutes, dfRun, output_path = output_route_completion_path)

# Important for sensitivity analysis to compare the same set of trips between runs.
# To do this need to exclude peds ID that get removed from all other runs
dfPedRoutesConsistentPeds = dfPedRoutes.loc[ ~dfPedRoutes['ID'].isin(dfPedRoutes_removedpeds['ID'])]
dfCrossEventsConsistentPeds = dfCrossEvents.loc[ ~dfCrossEvents['ID'].isin(dfPedRoutes_removedpeds['ID'])]

dfLinkCrossCounts = bd_utils.get_road_link_pedestrian_crossing_counts(dfCrossEventsConsistentPeds, gdfPaveLinks)


print("\nCalculating/Loading Output Metrics")
#dfRouteLength = bd_utils.get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, output_path = output_route_length_file)

dfPedTripDD = bd_utils.agg_trip_distance_and_duration(dfPedRoutes_removedpeds['ID'], dfRun, ped_routes_file, output_ped_distdurs_file)
dfVehTripDD = bd_utils.agg_trip_distance_and_duration(None, dfRun, veh_routes_file, output_veh_distdurs_file)
dfRouteLength = bd_utils.get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, output_path = output_route_length_file)
dfCrossLocEntropy = bd_utils.calculate_average_link_level_crossing_location_entropy(dfCrossEventsConsistentPeds, dfPedRoutesConsistentPeds.reindex(columns = ['run','ID','node_path']), gdfPaveLinks, gdfPaveNodes, gdfORLinks, dfRun, nbins = nbins, bin_dist = bin_dist, output_path = output_link_cross_entropy)
dfCrossCounts = dfCrossEventsConsistentPeds.merge(dfRun.reindex(columns = ['run','nPeds']), on='run').groupby('run').apply(lambda df: df.shape[0] / df['nPeds'].values[0]).reset_index().rename(columns = {0:'crossCountPP'})
dfCrossLocEntropy = pd.merge(dfCrossLocEntropy, dfCrossCounts, on='run', how = 'outer')

# Helpful to visualise the crossing coordiantes
'''
dfCrossEventsBins = bd_utils.get_crossing_locations_and_bins(dfCrossEvents, dfPedRoutesConsistentPeds.reindex(columns = ['run','ID','node_path']), gdfPaveLinks, gdfPaveNodes, gdfORLinks, nbins)
gdfCrossEventsBins = gpd.GeoDataFrame(dfCrossEventsBins, geometry = 'rl_cross_point', crs = {'init' :'epsg:27700'})
gdfCrossEventsBins = pd.merge(gdfCrossEventsBins, dfRun, on='run')
gdfCrossEventsBins.drop(['fid', 'node_path'], axis=1, inplace=True)
gdfCrossEventsBins.loc[ gdfCrossEventsBins['informalCrossing']==True].to_file(os.path.join(data_dir, 'cross_locs_informal.{}.{}.gpkg'.format(nbins, file_datetime_string)), driver='GPKG')
gdfCrossEventsBins.loc[ gdfCrossEventsBins['informalCrossing']==False].to_file(os.path.join(data_dir, 'cross_locs_no_informal.{}.{}.gpkg'.format(nbins, file_datetime_string)), driver='GPKG')
'''
dfConflicts = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds, dfRun, dfLinkCrossCounts, output_cross_conflicts, ttc_col = 'TTC', ttc_threshold=ttc_threshold)
#dfConflictsUnmarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unmarked'], dfLinkCrossCounts, ttc_col = 'TTC')
#dfConflictsDiagonalUm = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ (dfCrossEventsConsistentPeds['linkType']=='diag_cross') & (dfCrossEventsConsistentPeds['CrossingType']=='unmarked')], dfLinkCrossCounts, ttc_col = 'TTC')
#conflicts_data = {'all':dfConflicts, 'unmarked':dfConflictsUnmarked, 'diag_um':dfConflictsDiagonalUm}

dfPedTL = bd_utils.calculate_ped_trip_distance(dfPedRoutesConsistentPeds, dfCrossEventsConsistentPeds, gdfPaveLinks, gdfPaveNodes, dfRun, output_path = output_ped_trip_length)

print("\nAggregating Metrics for Policy Analysis")

# Merge pedestrian and vehicle distances and durations together
#dfDD = pd.merge(dfPedTripDD, dfVehTripDD.reindex(columns = ['run','DistPA','DurPA']), on='run', indicator=True, how = 'outer', suffixes = ("Ped", "Veh"))
#assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
#dfDD.drop('_merge', axis=1, inplace=True)

# Merge in crossing location entropy data and ped distance travelled
dfDD = pd.merge(dfRouteLength, dfCrossLocEntropy.reindex(columns = ['run','crossCountPP','mean_link_cross_entropy']), on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)
dfDD = pd.merge(dfDD, dfPedTripDD.reindex(columns = ['run','DistPA']), on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)
dfDD = pd.merge(dfDD, dfVehTripDD.reindex(columns = ['run','speed']).rename(columns={'speed':'speedVeh'}), on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)
dfDD = pd.merge(dfDD, dfConflicts.reindex(columns = ['run','conflict_count']), on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)
dfDD = pd.merge(dfDD, dfPedTL, on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)


dfDD.to_csv(output_sd_data, index=False)