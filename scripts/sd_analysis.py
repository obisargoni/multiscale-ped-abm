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
vehicle_density_timestamp = config['vehicle_density_timestamp']
setting = config['setting']

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
nbins = config['crossing_nbins']

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

output_paths = bd_utils.get_ouput_paths(file_datetime_string, vehicle_density_timestamp, data_dir, nbins = nbins)
output_ped_routes_file = output_paths["output_ped_routes_file"]
output_route_length_file = output_paths["output_route_length_file"]
output_ped_distdurs_file = output_paths["output_ped_distdurs_file"]
output_veh_distdurs_file = output_paths["output_veh_distdurs_file"]
output_sp_similarity_length_path = output_paths["output_sp_similarity_length_path"]
output_cross_events_path = output_paths["output_cross_events_path"]
output_cross_entropy = output_paths["output_cross_entropy"]

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

#dfLinkCrossCounts = bd_utils.get_road_link_pedestrian_crossing_counts(dfCrossEventsConsistentPeds, gdfPaveLinks)


print("\nCalculating/Loading Output Metrics")
#dfRouteLength = bd_utils.get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, output_path = output_route_length_file)

dfPedTripDD = bd_utils.agg_trip_distance_and_duration(dfPedRoutes_removedpeds['ID'], dfRun, ped_routes_file, output_ped_distdurs_file)
dfVehTripDD = bd_utils.agg_trip_distance_and_duration(None, dfRun, veh_routes_file, output_veh_distdurs_file)
dfCrossLocEntropy = bd_utils.calculate_crossing_location_entropy(dfCrossEventsConsistentPeds, dfPedRoutesConsistentPeds.reindex(columns = ['run','ID','node_path']), gdfPaveLinks, gdfPaveNodes, gdfORLinks, dfRun, nbins = nbins, output_path = output_cross_entropy)

# Helpful to visualise the crossing coordiantes
'''
dfCrossEventsBins = bd_utils.get_crossing_locations_and_bins(dfCrossEvents, dfPedRoutesConsistentPeds.reindex(columns = ['run','ID','node_path']), gdfPaveLinks, gdfPaveNodes, gdfORLinks, nbins)
gdfCrossEventsBins = gpd.GeoDataFrame(dfCrossEventsBins, geometry = 'rl_cross_point', crs = {'init' :'epsg:27700'})
gdfCrossEventsBins = pd.merge(gdfCrossEventsBins, dfRun, on='run')
gdfCrossEventsBins.drop(['fid', 'node_path'], axis=1, inplace=True)
gdfCrossEventsBins.loc[ gdfCrossEventsBins['informalCrossing']==True].to_file(os.path.join(data_dir, 'cross_locs_informal.{}.{}.gpkg'.format(nbins, file_datetime_string)), driver='GPKG')
gdfCrossEventsBins.loc[ gdfCrossEventsBins['informalCrossing']==False].to_file(os.path.join(data_dir, 'cross_locs_no_informal.{}.{}.gpkg'.format(nbins, file_datetime_string)), driver='GPKG')
'''
#dfConflicts = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds, dfLinkCrossCounts, ttc_col = 'TTC')
#dfConflictsUnmarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unmarked'], dfLinkCrossCounts, ttc_col = 'TTC')
#dfConflictsDiagonalUm = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ (dfCrossEventsConsistentPeds['linkType']=='diag_cross') & (dfCrossEventsConsistentPeds['CrossingType']=='unmarked')], dfLinkCrossCounts, ttc_col = 'TTC')
#conflicts_data = {'all':dfConflicts, 'unmarked':dfConflictsUnmarked, 'diag_um':dfConflictsDiagonalUm}


print("\nAggregating Metrics for Policy Analysis")

# Merge pedestrian and vehicle distances and durations together
dfDD = pd.merge(dfPedTripDD, dfVehTripDD.reindex(columns = ['run','DistPA','DurPA']), on='run', indicator=True, how = 'outer', suffixes = ("Ped", "Veh"))
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)

# Merge in crossing location entropy data
dfDD = pd.merge(dfDD, dfCrossLocEntropy.reindex(columns = ['run','cross_entropy']), on='run', indicator=True, how = 'outer')
assert dfDD.loc[ dfDD['_merge']!='both'].shape[0]==0
dfDD.drop('_merge', axis=1, inplace=True)

dfDD.to_csv(output_sd_data, index=False)




########################################
#
#
# Processing data to compare runs under different policies
#
# This involves matching up runs with the same parameter settings but under different policy conditions
#
#######################################
print("\nRunning Scenario Discovery Analysis")

dfDD = pd.read_csv(output_sd_data)

# get policy parameter and split the data into groups for different policies
policy_param = list(policies.keys())[0]
policy_values = policies[policy_param]
scenario_param_cols =  [i for i in params if i!=policy_param]

# Now group by scenario and aggregate to find difference in outputs between policy conditions
for c in scenario_param_cols:
	dfDD[c] = dfDD[c].astype(str) # Helps with grouping, makes matching doubles easier

dfPolicyDiff = dfDD.groupby(scenario_param_cols).agg( 	PedDistDiff = pd.NamedAgg(column = "DistPAPed", aggfunc=lambda s: s.values[0] - s.values[1]),
														VehDistDiff = pd.NamedAgg(column = "DistPAVeh", aggfunc=lambda s: s.values[0] - s.values[1]),
														PedDurDiff = pd.NamedAgg(column = "DurPAPed", aggfunc=lambda s: s.values[0] - s.values[1]),
														VehDurDiff = pd.NamedAgg(column = "DurPAVeh", aggfunc=lambda s: s.values[0] - s.values[1]),
														CrossEntDiff = pd.NamedAgg(column = "cross_entropy", aggfunc=lambda s: s.values[0] - s.values[1]),
														PedDistDiffFrac = pd.NamedAgg(column = "DistPAPed", aggfunc=lambda s: (s.values[0] - s.values[1]) / s.values[0]),
														CrossEntDiffFrac = pd.NamedAgg(column = "cross_entropy", aggfunc=lambda s: (s.values[0] - s.values[1]) / s.values[0]),
														CountRuns = pd.NamedAgg(column = "run", aggfunc=lambda s: s.shape[0]),
														RunsStr = pd.NamedAgg(column = "run", aggfunc=lambda s: ":".join(str(i) for i in s.tolist())),
													).reset_index()

for c in scenario_param_cols:
	dfPolicyDiff[c] = dfPolicyDiff[c].astype(float)

# Check that there are expected number of runs per scenario
assert (dfPolicyDiff['CountRuns']==2).all()

##############################
#
#
# Data mining techniques applied to results to distinguish scenarios
#
#
##############################

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns

assert 'latin' in config['setting'] # expect LH desig to be used when doing SD

with open("figure_config.json") as f:
    fig_config = json.load(f)


#
# Initial exploratory analysis of multiple outcomes
#

experiments = dfDD.loc[:, params]

outcome_vars = ['DistPAPed','cross_entropy']
outcomes = {k:dfDD[k].values for k in outcome_vars}

#
# Create pairs plot
#
data = dfDD.loc[:, outcome_vars].rename(columns = {"DistPAPed":"Average Pedestrian Trip Distance", "cross_entropy":"Crossing Location Entropy"})
data['Informal Crossing'] = experiments['informalCrossing']
sns.pairplot(data, hue='Informal Crossing', vars=["Average Pedestrian Trip Distance","Crossing Location Entropy"])
plt.savefig(os.path.join(img_dir, 'pair_plot.{}bins.{}.png'.format(nbins, file_datetime_string)))



#
# paths heatmap data
#

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
n_paths = dfPedRoutes['FullStrategicPathString'].unique().shape[0]
edge_traverse_counts = dfPathORLinks['FullStrategicPathString'].value_counts()

path_links = dfPathORLinks['FullStrategicPathString'].unique()
edgedata = np.array([edge_traverse_counts[i] / n_paths for i in path_links])

#tp_links = edge_traverse_counts.index()

edgelist = []
for edge_id in path_links:
    for e in list(or_graph.edges(data=True)):
        if edge_id == e[-1]['fid']:
            edgelist.append(e)


# Get start and end nodes
gdfStartNodes = gdfPaveNodes.loc[ gdfPaveNodes['fid'].isin(dfPedRoutes['start_node'])]
gdfEndNodes = gdfPedODs.loc[gdfPedODs['fid'] =='od_0']

#
# KDE + paths heatmap plot
#
#plt.style.use('dark_background')
fig = plt.figure(figsize=(20,20))
gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1, 2], width_ratios = [1,1])

#  Varying density along a streamline
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[1, :])

sns.kdeplot(ax=ax0, data=data, x="Average Pedestrian Trip Distance", hue="Informal Crossing")
sns.kdeplot(ax=ax1, data=data, x="Crossing Location Entropy", hue="Informal Crossing")


ax2 = bd_utils.figure_rl_paths_heatmap(fig, ax2, gdfORLinks, gdfStartNodes, gdfEndNodes, or_graph, dict_node_pos, edgelist, edgedata, plt.get_cmap('Reds'), 'Paths Heatmap', "Frequency of link traversal divided by total number of trips", {'fontsize':15}, 15, fig_config, cbar_pad = 0.03)

fig.savefig(os.path.join(img_dir, 'kde_path_heatmap.{}bins.{}.png'.format(nbins, file_datetime_string)))
plt.style.use('default')