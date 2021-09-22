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
from scipy import stats
import matplotlib.pyplot as plt

import batch_data_utils as bd_utils

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
                break

    # Work way backwards through edges to reconstruct path
    node_path = [end_node]
    prev = end_node
    edge_path_reverse = edge_path[::-1]
    for i, e in enumerate(edge_path_reverse):
        candidate = adjacent_edge_node(prev, e)

        # The candidate node wll be none in cases where the edge path did not get to the destination node due to the ped being removed from the simulation.
        if candidate is None:
            break

        # Decide whether this node was part of the path by looking at next edge
        if i == len(edge_path_reverse)-1:
            # If at the end of the edge path need different rule

            if candidate == start_node:
                # In this case the first pavement link is included in the edge path and so the path ends at the start node
                # Set next_candidate to be the candidate node so that the candidate node is added to the node path
                next_candidate = candidate
            else:
                # Check if the candidate is adjacent to the start node. if not then this edge is a default edge not traversed. Don't include in route.
                if start_node in pavement_graph.neighbors(candidate):
                    next_candidate = candidate
                else:
                    next_candidate = None

        else:
            # Now check that candidate connects to next edge
            next_candidate = adjacent_edge_node(candidate, edge_path_reverse[i+1])

        if (next_candidate is not None) & (candidate not in node_path):
             node_path.append(candidate)
             prev = candidate

        if prev == start_node:
            break
        

    # If ped agent got removed from the simulation it would not have reached the end node and the path will not have any other nodes added
    # In this case return empty list
    if len(node_path)==1:
        return []

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

    return tuple(dijkstra_path)

def node_paths_overlap(npa, npb, normalised=True):
    '''An alternative methods for comparing paths expressed as sequences of nodes. 
    Simply count number of nodes that are not in both paths and optionally normalise by
    total number of nodes in both paths. A simpler method that is not influecenced by node position.
    '''

    n_diff = len(set(npa).symmetric_difference(set(npb)))
    if normalised==False:
        return n_diff
    else:
        n_tot = len(npa) + len(npb)
        return float(n_diff) / n_tot

def compare_node_paths(graph, npa, npb, dict_node_pos, distance_function = sim.frechet_dist, account_for_path_length = False, weight = None):
    '''

    accound_for_path_length: It's possible that the two paths have the same path length but different nodes. This occurs when there are multiple shortest paths.
                            By checking the path length identify this and set distance between paths to 0 if path lengths are equal.

    weight: If account_for_path_length is True this gives the weight attribute to use to calculate path weight.
    '''

    if distance_function is None:
        d = node_paths_overlap(npa, npb)
    else:
        pos_a = [dict_node_pos[i] for i in npa]
        pos_b = [dict_node_pos[i] for i in npb]
        d = distance_function(pos_a, pos_b)

    if account_for_path_length:
        lengtha = nx.path_weight(graph, npa, weight=weight)
        lengthb = nx.path_weight(graph, npb, weight=weight)
        if abs(lengtha-lengthb)<0.00000001:
            d = 0.0
    return d

def get_pedestrian_run_durations(dfPedCrossings):
    # Get duration for each run - defined as total time pedestrians are in the simulation.
    dfPedStart = dfPedCrossings.groupby('run')['tick'].min().reset_index()
    dfPedEnd = dfPedCrossings.groupby('run')['tick'].max().reset_index()
    dfDurs = pd.merge(dfPedStart, dfPedEnd, on = 'run', suffixes = ('_start', '_end'))
    dfDurs['duration'] = dfDurs['tick_end'] - dfDurs['tick_start']
    return dfDurs

def get_road_link_vehicle_density(dfRunDurations, gdfITNLinks, data_file, output_path):

    if os.path.exists(output_path) == False:
        # Alternatative method for getting average vehicle counts
        dfVehRls = pd.read_csv(data_file)

        # Merge with start and end pedestrian times to and remove rows that lie outside these times
        dfVehRls = pd.merge(dfVehRls, dfRunDurations, on = 'run')
        dfVehRls = dfVehRls.loc[ (dfVehRls['tick']>=dfVehRls['tick_start']) & (dfVehRls['tick']<=dfVehRls['tick_end'])]

        # Merge with ITN links to get lookup to ped rl ID
        gdfITNLinks = gdfITNLinks.reindex(columns = ['fid','pedRLID'])
        dfVehRls = dfVehRls.reindex(columns = ['tick','run','ID','duration','CurrentRoadLinkID'])
        dfVehRls = pd.merge( dfVehRls, gdfITNLinks, left_on = 'CurrentRoadLinkID', right_on = 'fid', how = 'left')

        # get average count of vehicles on each ped road link, by first getting count per tick then summing and averaging this.
        VehCountTick = dfVehRls.groupby(['run', 'duration', 'pedRLID','tick'])['ID'].apply(lambda ids: ids.unique().shape[0]).reset_index().rename(columns = {'ID':'VehCount'})
        dfVehRls = None
        gdfITNLinks = None

        VehCountAv = VehCountTick.groupby(['run', 'pedRLID']).apply( lambda df: df['VehCount'].sum() / df['duration'].values[0]).reset_index().rename(columns = {0:'AvVehCount'})

        VehCountAv.to_csv(output_path, index=False)
    else:
        VehCountAv = pd.read_csv(output_path)
    
    return VehCountAv

def get_ped_routes(dfPedCrossings, gdfPaveLinks, weight_params):

    print("\nExtracting Pedestrian Agent Routes")
    
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
    dfPedStratPaths = dfPedCrossings.groupby(['run','ID'])['FullStrategicPathString'].apply(lambda s: tuple(s.drop_duplicates().values[0].split(':')[1:])).reset_index()
    dfPedRoutes = pd.merge(dfPedRoutes, dfPedStratPaths, on = ['run', 'ID'])

    dfPedOD = dfPedCrossings.groupby(['run', 'ID'])['StartPavementJunctionID', 'DestPavementJunctionID'].apply(lambda df: df.drop_duplicates()).reset_index()
    dfPedRoutes = pd.merge(dfPedRoutes, dfPedOD, left_on = ['run', 'ID'], right_on = ['run', 'ID'])

    ## Need to keep these in or keep a record of what rows got dropped in order to calculate SIs later
    # Or avoid dropping.
    dfPedRoutes['node_path'] = dfPedRoutes.apply(lambda row: node_path_from_edge_path(row['edge_path'], row['StartPavementJunctionID'], row['DestPavementJunctionID'], pavement_graph), axis=1)
    dfPedRoutes_removedpeds = dfPedRoutes.loc[ dfPedRoutes['node_path'].map(lambda x: len(x)==0)]
    removed_peds_index = dfPedRoutes_removedpeds.index

    ######################################
    #
    #
    # Identify unique OD trips
    #
    #
    ######################################

    # In some cases the first pavement edge does not get included in the data. I think due to the ped origin being v close to pavement junction
    # In this case need to alter what the start node is
    # So if node path doesn't include start node, calculate alternative path using first node in path.

    print("Checking for node_path missing the start pavement node (excluding peds removed from simulation)")
    dfPedRoutes['missing_start_node'] = np.nan
    dfPedRoutes.loc[~dfPedRoutes.index.isin(removed_peds_index), 'missing_start_node'] = dfPedRoutes.loc[~dfPedRoutes.index.isin(removed_peds_index)].apply(lambda row: row['StartPavementJunctionID'] not in row['node_path'], axis=1)
    print(dfPedRoutes['missing_start_node'].value_counts(dropna=True))

    print("Checking instances where the first node in node_path matches the start pavement node (excluding peds removed from simulation)")
    dfPedRoutes['sn_at_start'] = np.nan
    dfPedRoutes.loc[~dfPedRoutes.index.isin(removed_peds_index), 'sn_at_start'] = dfPedRoutes.loc[~dfPedRoutes.index.isin(removed_peds_index)].apply(lambda row: row['StartPavementJunctionID'] == row['node_path'][0], axis=1)
    print(dfPedRoutes['sn_at_start'].value_counts(dropna=True))

    # Use first not in path as start node, as opposed to the node the ped actually starts their journey from since these only differ in a handful of cases.
    dfPedRoutes['start_node'] = dfPedRoutes.apply(lambda row: row['node_path'][0] if len(row['node_path'])>0 else row['StartPavementJunctionID'], axis=1)
    dfPedRoutes['end_node'] = dfPedRoutes.apply(lambda row: row['node_path'][-1] if len(row['node_path'])>0 else row['DestPavementJunctionID'], axis=1)

    # Find the unique set of start and end node and calculate shortest paths between these. Then merge into the ped routes data.
    dfUniqueStartEnd = dfPedRoutes.loc[:, ['start_node', 'end_node', 'FullStrategicPathString']].drop_duplicates()

    ######################################
    #
    #
    # Calculate shortest paths using weight accounting for vehicle density.
    #
    # Compare to the distance weighted shortest path by comparing the means of the frechet distances between shortest path and the pedestrians actual path.
    #
    ######################################

    # Set dfPedRoutes to just have columns we are interested in
    dfPedRoutes = dfPedRoutes.reindex(columns = [   'run', 'ID', 'edge_path', 'FullStrategicPathString',
                                                    'StartPavementJunctionID', 'DestPavementJunctionID', 'node_path',
                                                    'start_node', 'end_node'])

    for k in weight_params:
        weight_name = "weight{}".format(k)
        gdfPaveLinks['cross_cost'] = gdfPaveLinks['AvVehDen'].fillna(0) * k
        gdfPaveLinks[weight_name] = gdfPaveLinks['length'] + gdfPaveLinks['cross_cost']

        # Weight pavement network crossing links by average vehicle flow
        weight01_attributes = gdfPaveLinks.set_index( gdfPaveLinks.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
        nx.set_edge_attributes(pavement_graph, weight01_attributes, name = weight_name)

        dfUniqueStartEnd['sp_{}'.format(k)] = dfUniqueStartEnd.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'], weight = weight_name), axis=1)


    dfPedRoutes = pd.merge(dfPedRoutes, dfUniqueStartEnd, on = ['start_node', 'end_node', 'FullStrategicPathString'])

    return dfPedRoutes, dfPedRoutes_removedpeds

def get_shortest_path_similarity(dfPedRoutes, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = None, output_path = "sp_similarity.csv"):

    ######################################
    #
    #
    # Compare these various shortest path routes to the ABM routes by calculating a metric of difference between their node paths
    #
    #
    ######################################

    if os.path.exists(output_path)==False:
        dfSPSim = pd.DataFrame()
        for k in weight_params:
            dfPedRoutes['comp_value_{}'.format(k)] = dfPedRoutes.apply(lambda row: compare_node_paths(pavement_graph, row['node_path'], row['sp_{}'.format(k)], dict_node_pos, distance_function = distance_function), axis=1)
            
            # This is only meaning full for the shortest path unweighted by vehicle traffic, where we can expect the path to match the ABM tactical path, and therefore compare path lengths to check for equivalence.
            dfPedRoutes['comp_path_weight_{}'.format(k)] = dfPedRoutes.apply(lambda row: compare_node_paths(pavement_graph, row['node_path'], row['sp_{}'.format(k)], dict_node_pos, distance_function = distance_function, account_for_path_length=True, weight='length'), axis=1)

            df = dfPedRoutes.groupby('run')['comp_value_{}'.format(k)].describe().reset_index()
            df['k'] = k

            # T test to compare means
            # Compare frechet distances for different edge weights. Used to test whether the additional weighting of crossing links better matches pedestrian routes.
            dfTTests = dfPedRoutes.groupby('run').apply(lambda df: stats.ttest_rel(df['comp_value_0'], df['comp_value_{}'.format(k)])[1]).reset_index().rename(columns = {0:'ttp'})
            df = pd.merge(df, dfTTests, on='run')

            dfSPSim = pd.concat([dfSPSim, df])

        # rebuild index
        dfSPSim.index = np.arange(dfSPSim.shape[0])

        # To make comparisons between the k=0 route comp mean and k>0 route comp means, merge the k=0 value in
        dfSPSim = pd.merge(dfSPSim, dfSPSim.loc[ dfSPSim['k']==0, ['run', 'mean']], on = 'run', suffixes = ('', '0'))

        # Identify cases where the weighted crossing link shortest path better matches the ABM path
        dfSPSim['is_sig_lower'] = (dfSPSim['mean0'] > dfSPSim['mean']) & (dfSPSim['ttp']<0.05)

        # Merge in parameter values
        dfSPSim = pd.merge(dfRun, dfSPSim, on = 'run')

        # Save date for future use
        dfSPSim.to_csv(output_path, index=False)
    else:
        dfSPSim = pd.read_csv(output_path)

    return dfSPSim

def agg_route_completions(dfPedRoutes, dfRun, output_path = 'route_completions.csv'):
    if os.path.exists(output_path)==False:
        dfPedRoutes['completed_journey'] = dfPedRoutes['node_path'].map(lambda x: int(len(x)>0))
        dfCompletions = dfPedRoutes.groupby('run').apply( lambda df: df['completed_journey'].sum() / float(df.shape[0])).reset_index().rename(columns = {0:'frac_completed_journeys'})
        dfCompletions = pd.merge(dfRun, dfCompletions, on = 'run')
        dfCompletions.to_csv(output_path, )
    else:
        dfCompletions = pd.read_csv(output_path)

    return dfCompletions

def factor_map(problem, X, Y, threshold):
    # Identify outputs that are above and below threshold
    b = np.where(Y > threshold)[0]
    b_ = np.where(~(Y>threshold))[0]

    # Compare parameter values between these two groups of scenarios using ks test
    ks_results = []
    for i, f in enumerate(problem['names']):
        Xi = X[b, i]
        Xi_ = X[b_, i]

        m = np.mean(Xi)
        s = np.std(Xi)
        m_ = np.mean(Xi_)
        s_ = np.std(Xi_)

        # KS test of similarity of distributions
        D, p_ks = stats.kstest(Xi, Xi_, alternative  = 'two_sided')

        # T test
        T, p_t = stats.ttest_ind(Xi, Xi_, alternative = 'two-sided')

        # Gather results into a dictionary
        ks_result = {'name':f, 'm':m, 's':s, 'm_':m_, 's_':s_, 'D':D, 'p_ks':p_ks, 'T':T, 'p_t':p_t}
        ks_results.append(ks_result)
    
    dfks = pd.DataFrame(ks_results).sort_values(by = 'p_ks')

    # Calculate correlations between parameters in behavioural group
    Xb = X[b,:]
    corrs = []
    for i, j in itertools.combinations(range(problem['num_vars']), 2):
        if i==j:
            continue

        Xbi = Xb[:,i]
        Xbj = Xb[:,j]
        lr = stats.linregress(Xbi, Xbj)

        corr = {'fi':problem['names'][i], 'fj':problem['names'][j], 'r':lr.rvalue, 'p':lr.pvalue, 'data':np.array([Xbi, Xbj])}
        corrs.append(corr)

    dfcorrs = pd.DataFrame(corrs)

    return dfks, dfcorrs

def morris_si_bar_figure(dfsi, fig_title):
    f, axs = plt.subplots(1,2, figsize = (20,10))
    axs[0].bar(range(dfsi.shape[0]), dfsi['mu_star'], width=0.8, yerr = dfsi['mu_star_conf'], align='center')
    axs[1].bar(range(dfsi.shape[0]), dfsi['sigma'], width=0.8, align='center')
    axs[0].set_xticks(range(dfsi.shape[0]))
    axs[1].set_xticks(range(dfsi.shape[0]))
    axs[0].set_xticklabels(dfsi['names'], rotation = 45)
    axs[1].set_xticklabels(dfsi['names'], rotation = 45)
    axs[0].set_title("mu star")
    axs[1].set_title("sigma")
    f.suptitle(fig_title)
    return f

#####################################
#
#
# Choose data to analyze and sensitivity analysis setting
#
#
#####################################
file_datetime_string = "2021.Sep.22.08_17_14"
vehicle_density_timestamp = "2021.Sep.20.15_18_00"
setting = 'monte_carlo_filtering'

#####################################
#
# File locations
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
itn_links_file = os.path.join(gis_data_dir, config['mastermap_itn_processed_direction_file'])
itn_nodes_file = os.path.join(gis_data_dir, config['mastermap_node_processed_file'])
crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])


# Model output data
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")
file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime)
ped_crossings_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_locations", file_datetime = file_datetime)
ped_locations_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("vehicle_road_links", file_datetime = file_datetime)
vehicle_rls_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime, suffix = 'batch_param_map')
batch_file = bd_utils.most_recent_directory_file(data_dir, file_re)

# output paths for processed data
output_vehicle_density_file = os.path.join(data_dir, "av_vehicle_density.{}.csv".format(vehicle_density_timestamp))
output_sp_similarity_path = os.path.join(data_dir, "sp_similarity.{}.csv".format(file_datetime_string))
output_route_completion_path = os.path.join(data_dir, "route_completions.{}.csv".format(file_datetime_string))

#####################################
#
#
# Load Data
#
#
#####################################

# Data from model run
dfPedCrossings = pd.read_csv(ped_crossings_file)
dfRun = pd.read_csv(os.path.join(data_dir, batch_file)).sort_values(by = 'run')

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

weight_params = range(0, 1000, 100)

######################################
#
#
# Process data to get ped routes and shortest path routes
#
#
######################################
# First process vehicle road link data, since this is largest file
dfRunDurations =  get_pedestrian_run_durations(dfPedCrossings)
VehCountAv = get_road_link_vehicle_density(dfRunDurations, gdfITNLinks, vehicle_rls_file, output_vehicle_density_file)
gdfPaveLinks = pd.merge(gdfPaveLinks, VehCountAv, left_on = 'pedRLID', right_on = 'pedRLID', how = 'left')
gdfPaveLinks.loc[ gdfPaveLinks['pedRLID'].isin(gdfCAs['roadLinkID'].unique()), 'AvVehDen'] = 0.0

dfPedRoutes, dfPedRoutes_removedpeds = get_ped_routes(dfPedCrossings, gdfPaveLinks, weight_params)
dfRouteCompletion = agg_route_completions(dfPedRoutes, dfRun, output_path = output_route_completion_path)
dfSPSim = get_shortest_path_similarity(dfPedRoutes, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = None, output_path = output_sp_similarity_path)


######################################
#
#
# Import problem definition used for sampling
#
#
######################################
from SALib.analyze import morris
import sys
sys.path.append(".\\sample")
from SALibRepastParams import num_levels, problem, random_seed

######################################
#
#
# Factor mapping
#
#
######################################

if setting == 'monte_carlo_filtering':
    print("\nPerforming Factor Mapping")
    X_rc = dfRouteCompletion.loc[:, problem['names']].values
    Y_rc = dfRouteCompletion.loc[:, 'frac_completed_journeys'].values

    # Identify factors that are significantly different between scenarios where peds complete journeys and those where they don't
    threshold = 0.0
    dfks, dfcorrs = factor_map(problem, X_rc, Y_rc, threshold)

######################################
#
#
# Calculate sensitivity indices
#
#
######################################

if setting == "morris_factor_fixing":

    print("\nCalculating sensitivity indices - Route completions")

    Sis = morris.analyze(problem, X_rc, Y_rc, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)

    # Gather into a dataframe
    dfcompsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_compsi = morris_si_bar_figure(dfcompsi, "Jouney Completion SIs")
    f_compsi.show()


    print("\nCalculating sensitivity indices - Comparison to shortest path")

    # Get array of parameter values and output values
    k = 0
    X = dfSPSim.loc[ dfSPSim['k']==k, problem['names']].values
    Y = dfSPSim.loc[ dfSPSim['k']==k, 'mean'].values

    Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)

    # Gather into a dataframe
    dfspsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_spsi = morris_si_bar_figure(dfspsi, "Shortest Path SIs")
    f_spsi.show()