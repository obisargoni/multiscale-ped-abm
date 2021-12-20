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

def node_path_dice_distance(npa, npb):
    '''Calculates the Sørensen–Dice index of similarity between two discrete sets. Returns 1-similarity score
    '''
    n_sim = len(set(npa).intersection(set(npb)))
    n_tot = len(npa) + len(npb)
    sim = (2*n_sim) / float(n_tot)
    return 1-sim

def compare_node_paths(graph, npa, npb, dict_node_pos, distance_function = 'frechet', account_for_path_length = False, weight = None):
    '''

    accound_for_path_length: It's possible that the two paths have the same path length but different nodes. This occurs when there are multiple shortest paths.
                            By checking the path length identify this and set distance between paths to 0 if path lengths are equal.

    weight: If account_for_path_length is True this gives the weight attribute to use to calculate path weight.
    '''

    if distance_function == 'dice_dist':
        d = node_path_dice_distance(npa, npb)
    elif distance_function == 'path_length':
        lengtha = nx.path_weight(graph, npa, weight=weight)
        lengthb = nx.path_weight(graph, npb, weight=weight)
        d = abs(lengtha-lengthb) / lengthb
    elif distance_function == 'frechet':
        pos_a = [dict_node_pos[i] for i in npa]
        pos_b = [dict_node_pos[i] for i in npb]
        d = sim.frechet_dist(pos_a, pos_b)
    else:
        return None


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
        gdfITNLinks = gdfITNLinks.reindex(columns = ['fid','pedRLID', 'length'])
        dfVehRls = dfVehRls.reindex(columns = ['tick','run','ID','duration','CurrentRoadLinkID'])
        dfVehRls = pd.merge( dfVehRls, gdfITNLinks, left_on = 'CurrentRoadLinkID', right_on = 'fid', how = 'left')

        # get average count of vehicles on each ped road link, by first getting count per tick then summing and averaging this.
        VehCountTick = dfVehRls.groupby(['run', 'duration', 'pedRLID','tick'])['ID'].apply(lambda ids: ids.unique().shape[0]).reset_index().rename(columns = {'ID':'VehCount'})

        # get total lenth of ITN links per each pedRLID
        AggITNLengths = gdfITNLinks.groupby(['pedRLID'])['length'].sum().reset_index()

        dfVehRls = None
        gdfITNLinks = None

        VehCountAv = VehCountTick.groupby(['run', 'pedRLID']).apply( lambda df: df['VehCount'].sum() / df['duration'].values[0]).reset_index().rename(columns = {0:'AvVehCount'})
        VehCountAv = pd.merge(VehCountAv, AggITNLengths, on = 'pedRLID')
        VehCountAv['AvVehDen'] = VehCountAv['AvVehCount'] / VehCountAv['length']
        VehCountAv.drop('length', axis=1, inplace=True)

        VehCountAv.to_csv(output_path, index=False)
    else:
        VehCountAv = pd.read_csv(output_path)

    return VehCountAv

def unstack_sp_string(df, sp_string_col = 'FullStrategicPathString'):
    links = df[sp_string_col].values[0].split(":")[1:]
    return pd.DataFrame({'pedRLID':links})

def get_road_link_pedestrian_counts(dfPedCrossings, gdfPaveLinks):

    dfUnstacked = dfPedCrossings.groupby(['run','ID']).apply(unstack_sp_string).reset_index()

    dfORPedCounts = dfUnstacked.groupby(['run','pedRLID'])['ID'].apply(lambda s: s.unique().shape[0]).reset_index()
    dfORPedCounts.rename(columns={'ID':'ped_count'}, inplace=True)

    # Now merge with lookup from or link to pave link
    dfORPedCounts = pd.merge(dfORPedCounts, gdfPaveLinks.reindex(columns = ['fid', 'pedRLID']).drop_duplicates(), on = 'pedRLID')

    return dfORPedCounts

def get_road_link_pedestrian_crossing_counts(dfCrossEvents, gdfPaveLinks):
    '''Count the number of pedestrians crossing on each pavement road link. Aggregate this to OR Road Link to get total number of crossings on that road link.
    Then provide lookup from pavement link to OR link so that crossings on a pavement link can be normalised by total number of crossings on the corresponding road link.
    '''

    dfCrossCounts = dfCrossEvents.groupby(['run','CurrentPavementLinkID'])['ID'].apply(lambda s: s.unique().shape[0]).reset_index()
    dfCrossCounts.rename(columns={'ID':'cross_count'}, inplace=True)

    # Now merge with lookup from or link to pave link
    dfCrossCounts = pd.merge(dfCrossCounts, gdfPaveLinks.reindex(columns = ['fid', 'pedRLID']).drop_duplicates(), left_on = 'CurrentPavementLinkID', right_on = 'fid')

    # Aggregate to get road link cross counts
    dfRLCrossCounts = dfCrossCounts.groupby(['run','pedRLID'])['cross_count'].apply(lambda s: s.sum()).reset_index()

    # Now merge this to lookup from OR link to pavement link
    dfRLCrossCounts = pd.merge(dfRLCrossCounts, gdfPaveLinks.reindex(columns = ['fid', 'pedRLID']).drop_duplicates(), on = 'pedRLID')

    return dfRLCrossCounts

def get_ped_routes(dfPedCrossings, gdfPaveLinks, weight_params, output_path = "ped_routes.csv"):

    split_path = os.path.splitext(output_path)
    routes_removed_path = split_path[0] + "_removed_peds" + split_path[1]

    if os.path.exists(output_path):
        # Load data from file
        dfPedRoutes = pd.read_csv(output_path)
        dfPedRoutes_removedpeds = pd.read_csv(routes_removed_path)

        # Convert string node paths to lists
        dfPedRoutes['node_path'] = dfPedRoutes['node_path'].map(lambda np: np.strip('"[').strip(']"').replace("'","").split(", "))
        dfPedRoutes_removedpeds['node_path'] = dfPedRoutes_removedpeds['node_path'].map(lambda np: np.strip('"[').strip(']"').replace("'","").split(", "))

        for k in weight_params:
            col = 'sp_{}'.format(0)
            dfPedRoutes[col] = dfPedRoutes[col].map(lambda np: np.strip('"(').strip(')"').replace("'","").split(", "))
    else:
        # Otherwise create data

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

        dfPedRoutes.to_csv(output_path, index=False)
        dfPedRoutes_removedpeds.to_csv(routes_removed_path, index=False)

    return dfPedRoutes, dfPedRoutes_removedpeds

def load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = "ped_routes.csv"):

    split_path = os.path.splitext(ped_routes_path)
    routes_removed_path = split_path[0] + "_removed_peds" + split_path[1]

    if os.path.exists(ped_routes_path)==False:
        print("\nPedestrian routes path not found.\n")

    else:
        # Otherwise create data

        print("\nLoading Pedestrian Agent Routes")


        # Load data from file
        dfPedRoutes = pd.read_csv(ped_routes_path)

        # Convert strategic path to list
        dfPedRoutes['FullStrategicPathString'] = dfPedRoutes['FullStrategicPathString'].map(lambda s: tuple(dict.fromkeys(s.strip(":").split(":"))))

        # Convert string tactical path to list
        dfPedRoutes['edge_path'] = dfPedRoutes['FullTacticalRouteString'].map(lambda s: tuple(dict.fromkeys(s.strip(":").split(":"))))

        # Convert edge path to node path
        dfPedRoutes['node_path'] = dfPedRoutes.apply(lambda row: node_path_from_edge_path(row['edge_path'], row['StartPavementJunctionID'], row['DestPavementJunctionID'], pavement_graph), axis=1)

        ## Need to keep these in or keep a record of what rows got dropped in order to calculate SIs later
        # Or avoid dropping.
        dfPedRoutes_removedpeds = dfPedRoutes.loc[ dfPedRoutes['node_path'].map(lambda x: x[-1])!=dfPedRoutes['DestPavementJunctionID']]
        removed_peds_index = dfPedRoutes_removedpeds.index

        # Find the unique set of start and end node and calculate shortest paths between these. Then merge into the ped routes data.
        dfPedRoutes['start_node'] = dfPedRoutes.apply(lambda row: row['node_path'][0] if len(row['node_path'])>0 else row['StartPavementJunctionID'], axis=1)
        dfPedRoutes['end_node'] = dfPedRoutes.apply(lambda row: row['node_path'][-1] if len(row['node_path'])>0 else row['DestPavementJunctionID'], axis=1)
        dfUniqueStartEnd = dfPedRoutes.loc[:, ['start_node', 'end_node', 'FullStrategicPathString']].drop_duplicates()

                # Set dfPedRoutes to just have columns we are interested in
        dfPedRoutes = dfPedRoutes.reindex(columns = [   'run', 'ID', 'FullStrategicPathString', 'edge_path',
                                                        'StartPavementJunctionID', 'DestPavementJunctionID', 'node_path',
                                                        'start_node', 'end_node'])

        for k in weight_params:
            weight_name = "weight{}".format(k)
            gdfPaveLinks['cross_cost'] = gdfPaveLinks['AvVehDen'].fillna(0) * k
            gdfPaveLinks[weight_name] = gdfPaveLinks['length'] + gdfPaveLinks['cross_cost']

            # Weight pavement network crossing links by average vehicle flow
            weight_attributes = gdfPaveLinks.set_index( gdfPaveLinks.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
            nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

            dfUniqueStartEnd['sp_{}'.format(k)] = dfUniqueStartEnd.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'], weight = weight_name), axis=1)


        dfPedRoutes = pd.merge(dfPedRoutes, dfUniqueStartEnd, on = ['start_node', 'end_node', 'FullStrategicPathString'])

        #dfPedRoutes.to_csv(output_path, index=False)
        #dfPedRoutes_removedpeds.to_csv(routes_removed_path, index=False)

    return dfPedRoutes, dfPedRoutes_removedpeds

def get_ped_cross_events(dfPedCrossings, gdfPaveLinks, output_path = "cross_events.csv"):
    '''Method to aggregate crossing events from the Ped Crossings dataset, to produce dataframe with single row per crossing event.
    Crossing event defined by crossing type, location and minimum TTC during crossing.
    '''

    if os.path.exists(output_path)==False:

        # Process raw ped crossings data to
        # - drop rows that don't correspond to a crossing event
        # - aggregate TTC data to choose the lowest TTC value per crossing event
        # - this produces dataset with 1 row per crossing event
        # - Drop duplicates, expect just one crossing per ped per pavement link
        # - Check that there is one crossing per link per ped

        # check there won't be any ttC data lost when excluding 'none' crossing events
        assert dfPedCrossings.loc[ (dfPedCrossings['ChosenCrossingTypeString']=='none') & (~dfPedCrossings['TTC'].isnull()) ].shape[0]==0
        assert dfPedCrossings['CurrentPavementLinkID'].isnull().any()==False

        dfCrossEvents = dfPedCrossings.loc[dfPedCrossings['ChosenCrossingTypeString']!='none'].reindex(columns = ['run', 'ID', 'FullStrategicPathString', 'ChosenCrossingTypeString', 'CurrentPavementLinkID', 'CrossingCoordsString', 'TTC'])
        dfCrossEvents = dfCrossEvents.drop_duplicates()

        # Group by run, ID and CurrentPavementLinkID to find crossing coordinates and lowest TTC
        dfCrossEvents['TTC'] = dfCrossEvents.groupby(['run', 'ID', 'ChosenCrossingTypeString', 'CurrentPavementLinkID'])['TTC'].transform(lambda s: s.min())
        dfCrossEvents['CrossingCoordsString'] = dfCrossEvents.groupby(['run', 'ID', 'ChosenCrossingTypeString', 'CurrentPavementLinkID'])['CrossingCoordsString'].transform(lambda s: max(s.dropna(), key=len) if ~s.isnull().all() else None)

        # Drop duplicates again now that TTC and crossing coord processed
        dfCrossEvents = dfCrossEvents.drop_duplicates()

        # Merge with pave links to get link type
        dfLinkTypes = gdfPaveLinks.reindex(columns = ['fid','linkType']).drop_duplicates()
        dfCrossEvents = pd.merge(dfCrossEvents, dfLinkTypes, left_on = 'CurrentPavementLinkID', right_on = 'fid', how = 'left')
        dfCrossEvents.drop('fid', axis=1,inplace=True)

        # Check that there is a single crossing event per ped per pavement link.
        cross_per_ped_link = dfCrossEvents.groupby(['run', 'ID', 'CurrentPavementLinkID']).apply(lambda df: df.shape[0])
        assert cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]==0

        dfCrossEvents.to_csv(output_path, index=False)
    else:
        dfCrossEvents = pd.read_csv(output_path)

    return dfCrossEvents

def get_shortest_path_similarity(dfPedRoutes, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'dice_dist', exclude_stuck_peds = True, output_path = "sp_similarity.csv"):

    ######################################
    #
    #
    # Compare these various shortest path routes to the ABM routes by calculating a metric of difference between their node paths
    #
    #
    ######################################

    if exclude_stuck_peds:
        stuck_peds_index = dfPedRoutes.loc[ dfPedRoutes['node_path'].map(lambda x: len(x)==0)].index
        dfPedRoutes.drop(stuck_peds_index, inplace=True)

    if os.path.exists(output_path)==False:
        dfSPSim = pd.DataFrame()
        for k in weight_params:
            dfPedRoutes['comp_value_{}'.format(k)] = dfPedRoutes.apply(lambda row: compare_node_paths(pavement_graph, row['node_path'], row['sp_{}'.format(k)], dict_node_pos, distance_function = distance_function, weight='length'), axis=1)

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

def get_run_total_route_length(dfPedRoutes, dfRun, pavement_graph, exclude_stuck_peds = True, output_path = "run_route_length.csv"):

    ######################################
    #
    #
    # Compare these various shortest path routes to the ABM routes by calculating a metric of difference between their node paths
    #
    #
    ######################################

    if exclude_stuck_peds:
        stuck_peds_index = dfPedRoutes.loc[ dfPedRoutes['node_path'].map(lambda x: len(x)==0)].index
        dfPedRoutes.drop(stuck_peds_index, inplace=True)

    if os.path.exists(output_path)==False:
        dfPedRoutes['route_length'] = dfPedRoutes.apply(lambda row: nx.path_weight(pavement_graph, row['node_path'], weight='length'), axis=1)

        dfRouteLength = dfPedRoutes.groupby('run')['route_length'].sum().reset_index()

        # Merge in parameter values
        dfRouteLength = pd.merge(dfRun, dfRouteLength, on = 'run')

        # Save date for future use
        dfRouteLength.to_csv(output_path, index=False)
    else:
        dfRouteLength = pd.read_csv(output_path)

    return dfRouteLength

def agg_route_completions(dfPedRoutes, dfRun, output_path = 'route_completions.csv'):
    if os.path.exists(output_path)==False:
        dfPedRoutes['completed_journey'] = dfPedRoutes['node_path'].map(lambda x: int(len(x)>0))
        dfCompletions = dfPedRoutes.groupby('run').apply( lambda df: df['completed_journey'].sum() / float(df.shape[0])).reset_index().rename(columns = {0:'frac_completed_journeys'})
        dfCompletions = pd.merge(dfRun, dfCompletions, on = 'run')
        dfCompletions.to_csv(output_path, )
    else:
        dfCompletions = pd.read_csv(output_path)

    return dfCompletions

def agg_cross_conflicts(dfCrossEvents, dfLinkCrossCounts, ttc_col = 'TTC', ttc_threshold = 3):
    '''Aggregate crossing events to create indicators of conflict for each run. This involves findings the total number of conflicts per run and the
    mean TTC per run.
    '''

    # Join with pavement links data and aggregate to OR road link ID to get number of peds per OR road link for normalising crossing counts


    calc_conflict_count = lambda s: s.dropna().loc[s<ttc_threshold].shape[0]
    calc_mean_ttc = lambda s: s.dropna().loc[s<ttc_threshold].mean()
    calc_var_ttc = lambda s: s.dropna().loc[s<ttc_threshold].var()

    dfRunConflicts = dfCrossEvents.groupby("run").agg(  conflict_count=pd.NamedAgg(column=ttc_col, aggfunc=calc_conflict_count),
                                                        meanTTC=pd.NamedAgg(column=ttc_col, aggfunc=calc_mean_ttc),
                                                        varTTC=pd.NamedAgg(column=ttc_col, aggfunc=calc_var_ttc),
                                                        ).reset_index()

    # get conflict counts normalsied by numbers of peds on road links
    dfPaveLinkConflictCounts = dfCrossEvents.groupby(['run','CurrentPavementLinkID']).agg( conflict_count=pd.NamedAgg(column=ttc_col, aggfunc=calc_conflict_count),).reset_index()
    dfPaveLinkConflictCounts = pd.merge(dfPaveLinkConflictCounts, dfLinkCrossCounts, left_on = ['run', 'CurrentPavementLinkID'], right_on = ['run', 'fid'], indicator=True)
    assert dfPaveLinkConflictCounts.loc[ dfPaveLinkConflictCounts['_merge']!='both'].shape[0]==0
    dfPaveLinkConflictCounts['norm_conflict_count'] = dfPaveLinkConflictCounts['conflict_count'] / dfPaveLinkConflictCounts['cross_count']

    # now aggregate trun
    dfRunConflictsNorm = dfPaveLinkConflictCounts.groupby("run").agg(   meanNormCC=pd.NamedAgg(column='norm_conflict_count', aggfunc=np.mean),
                                                                        varNormCC =pd.NamedAgg(column='norm_conflict_count', aggfunc=np.var),).reset_index()

    dfRunConflicts = pd.merge(dfRunConflicts, dfRunConflictsNorm, on='run')

    # Merge with run params
    dfRunConflicts = pd.merge(dfRun, dfRunConflicts, on='run')

    return dfRunConflicts

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

def morris_si_bar_figure_w_sigma(dfsi, fig_title):
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

def morris_si_bar_figure(dfsi, fig_title, ylabel, xticklabels):
    f, ax = plt.subplots(1,1, figsize = (10,10))
    ax.bar(range(dfsi.shape[0]), dfsi['mu_star'], width=0.8, yerr = dfsi['mu_star_conf'], align='center')
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels, rotation = 45)
    ax.set_ylabel(ylabel)
    f.suptitle(fig_title)
    return f

def sobol_si_bar_figure(dfsi, fig_title, xticklabels):
    f, ax = plt.subplots(1,1, figsize = (10,10))
    bar_width=0.4
    x_pos = np.arange(dfsi.shape[0])
    ax.bar(x_pos, dfsi['S1'], width=bar_width, yerr = dfsi['S1_conf'], align='center', label="S1")
    ax.bar(x_pos+bar_width, dfsi['ST'], width=bar_width, yerr = dfsi['ST_conf'], align='center', label="ST")

    ax.set_xticks(x_pos + bar_width / 2)
    ax.set_xticklabels(xticklabels, rotation=45)
    ax.legend()

    f.suptitle(fig_title)
    return f

def batch_run_scatter(df_data, groupby_columns, parameter_sweep_columns, value_col, rename_dict, cmap, title = None, cbarlabel = None, output_path = None):

    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of 'addVehicleTicks'
    p = len(df_data[groupby_columns[0]].unique())
    q = len(df_data[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    f,axs = plt.subplots(p, q, figsize=(20,10), sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            df = grouped.get_group(group_key)

            # Select the corresponding axis
            ax = axs[pi, qi]
            im = ax.scatter(df[parameter_sweep_columns[0]], df[parameter_sweep_columns[1]], c=df[value_col], cmap=cmap, norm = None, vmin=0.0, vmax=1.0)

            # optionally add line indiceting e-g region where threhold can be met
            e = np.linspace(min(df[parameter_sweep_columns[0]])+0.001, max(df[parameter_sweep_columns[0]]),60)
            g = 1 - 1/e
            inds = np.where(g>=0)[0]

            ax.plot(e[inds], g[inds], color='black')
            #im = ax.scatter(e, g, c='black')

            ax.set_ylabel(rename_dict[parameter_sweep_columns[1]])
            ax.set_xlabel(rename_dict[parameter_sweep_columns[0]])

    # Add colourbar
    smap = plt.cm.ScalarMappable(cmap='viridis', norm=None)
    cbar = f.colorbar(smap, ax=axs, fraction=0.1, shrink = 0.8)

    cbar_fontdict = {"size":14}
    cbar.ax.tick_params(labelsize=cbar_fontdict['size']-3)
    cbar.ax.set_ylabel(rename_dict[value_col], rotation=-90, labelpad = 15, fontdict = cbar_fontdict)

    # Now add text annotations to indicate the scenario
    for i in range(p):
        ki = key_indices[i, 0]
        group_key = keys[ki]
        ax = axs[i, 0]

        s = "{}".format(rename_dict[group_key[0]])
        plt.text(-0.25,0.5, s, fontsize = 11, transform = ax.transAxes)


    for j in range(q):
        ki = key_indices[-1, j]
        group_key = keys[ki]

        ax = axs[-1, j]

        s = "{}".format(rename_dict[group_key[1]])
        plt.text(0.45,-0.25, s, fontsize = 11, transform = ax.transAxes)

    if title is not None:
        f.suptitle(title, fontsize=16, y = 1)
    if cbarlabel is not None:
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    if output_path is not None:
        plt.savefig(output_path)

    return f, axs

def tile_rgba_data(df_group, row, col, value_col = 'unmarked_pcnt', alpha_col = None, cmap = plt.cm.viridis):

    # Values used for pixel colours
    colour_data = df_group.reindex(columns = [row, col, value_col]).set_index([row, col]).unstack().sort_index(ascending=False)

    # Get rgb data from values
    norm = plt.Normalize()
    rgba = cmap(norm(colour_data.values))

    # Replace alpha
    if alpha_col is not None:
        alpha_data = df_group.reindex(columns = [row, col, alpha_col]).set_index([row, col]).unstack()
        rgba[:,:,3] = alpha_data.values

    row_labels = colour_data.index
    col_labels = colour_data.columns.get_level_values(1)

    return rgba, row_labels, col_labels

def batch_run_tile_plot(df_data, groupby_columns, parameter_sweep_columns, value_col, rename_dict, cmap, title = None, cbarlabel = None, output_path = None, figsize=(20,10)):

    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of 'addVehicleTicks'
    p = len(df_data[groupby_columns[0]].unique())
    q = len(df_data[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    f,axs = plt.subplots(p, q, figsize=figsize, sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            df = grouped.get_group(group_key)

            # get extent of x and y values, helps with setting axis ticks
            x_min, x_max = df[parameter_sweep_columns[0]].min(), df[parameter_sweep_columns[0]].max()
            y_min, y_max = df[parameter_sweep_columns[1]].min(), df[parameter_sweep_columns[1]].max()

            extent = [x_min , x_max, 0 , 2]

            # Select the corresponding axis
            ax = axs[pi, qi]

            rgba_data, row_labels, col_labels = tile_rgba_data(df, parameter_sweep_columns[1], parameter_sweep_columns[0], value_col = value_col, alpha_col = None, cmap = cmap)

            # Plot the tile
            im = ax.imshow(rgba_data, extent=extent)

            # Create Major and Minor ticks

            from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_major_formatter('{x:.0f}')
            ax.xaxis.set_minor_locator(MultipleLocator(0.25))

            #ax.set_xticks(np.arange(rgba_data.shape[1]))
            #ax.yaxis.set_major_locator(MultipleLocator(1.5))
            #ax.yaxis.set_major_formatter('{x:.0f}')
            ax.set_yticks([0.5, 1.5])
            ax.set_yticklabels(row_labels[::-1])

            # ... and label them with the respective list entries.
            #ax.set_xticklabels(col_labels)
            #ax.set_yticklabels(row_labels)

            '''
            # optionally add line indiceting e-g region where threhold can be met
            e = np.linspace(min(df[parameter_sweep_columns[0]])+0.001, max(df[parameter_sweep_columns[0]]),60)
            g = 1 - 1/e
            inds = np.where(g>=0)[0]

            ax.plot(e[inds], g[inds], color='black')
            '''
            e08 = 1 / (1-0.8)
            e09 = 1 / (1-0.9)

            ax.plot([e08,e08], [0,1],color='black',linewidth=2)
            ax.plot([e09,e09], [1,2],color='black',linewidth=2)

            ax.set_ylabel(rename_dict[parameter_sweep_columns[1]])
            ax.set_xlabel(rename_dict[parameter_sweep_columns[0]])

    # Add colourbar
    smap = plt.cm.ScalarMappable(cmap='viridis', norm=None)
    cbar = f.colorbar(smap, ax=axs, fraction=0.1, shrink = 0.8)

    cbar_fontdict = {"size":14}
    cbar.ax.tick_params(labelsize=cbar_fontdict['size']-3)
    cbar.ax.set_ylabel(rename_dict[value_col], rotation=-90, labelpad = 15, fontdict = cbar_fontdict)

    # Now add text annotations to indicate the scenario
    for i in range(p):
        ki = key_indices[i, 0]
        group_key = keys[ki]
        ax = axs[i, 0]

        s = "{}".format(rename_dict[group_key[0]])
        plt.text(-0.25,0.5, s, fontsize = 11, transform = ax.transAxes)


    for j in range(q):
        ki = key_indices[-1, j]
        group_key = keys[ki]

        ax = axs[-1, j]

        s = "{}".format(rename_dict[group_key[1]])
        plt.text(0.45,-0.55, s, fontsize = 11, transform = ax.transAxes)

    if title is not None:
        f.suptitle(title, fontsize=16, y = 1)
    if cbarlabel is not None:
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    if output_path is not None:
        plt.savefig(output_path)

    return f, axs


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


#####################################
#
# File locations
#
#####################################
project_crs = {'init': 'epsg:27700'}

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

file_re = bd_utils.get_file_regex("pedestrian_routes", file_datetime = file_datetime)
ped_routes_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_locations", file_datetime = file_datetime)
ped_locations_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("vehicle_road_links", file_datetime = file_datetime)
vehicle_rls_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime, suffix = 'batch_param_map')
batch_file = bd_utils.most_recent_directory_file(data_dir, file_re)

# output paths for processed data
output_ped_routes_file = os.path.join(data_dir, "ped_routes.{}.csv".format(file_datetime_string))
output_vehicle_density_file = os.path.join(data_dir, "av_vehicle_density.{}.csv".format(vehicle_density_timestamp))
output_route_length_file = os.path.join(data_dir, "run_route_length.{}.csv".format(file_datetime_string))
output_sp_similarity_path = os.path.join(data_dir, "sp_similarity.{}.csv".format(file_datetime_string))
output_sp_similarity_length_path = os.path.join(data_dir, "path_length_sp_similarity.{}.csv".format(file_datetime_string))
output_route_completion_path = os.path.join(data_dir, "route_completions.{}.csv".format(file_datetime_string))
output_cross_events_path = os.path.join(data_dir, "cross_events.{}.csv".format(file_datetime_string))
output_ks_factormap = os.path.join(data_dir , "ks_factor_map.{}.csv".format (file_datetime_string))
output_corr_factormap = os.path.join(data_dir , "corr_factor_map.{}.csv".format (file_datetime_string))

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

weight_params = range(0, 100, 100)

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

dfPedRoutes, dfPedRoutes_removedpeds = get_ped_routes(dfPedCrossings, gdfPaveLinks, weight_params, output_ped_routes_file)
dfCrossEvents = get_ped_cross_events(dfPedCrossings, gdfPaveLinks, output_path = output_cross_events_path)

dfPedRoutes, dfPedRoutes_removedpeds = load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = ped_routes_file)
dfCrossEvents = load_and_clean_cross_events(gdfPaveLinks, cross_events_path = cross_events_file)

# Data aggregated to run level, used to calculate sensitivity indices
dfRouteCompletion = agg_route_completions(dfPedRoutes, dfRun, output_path = output_route_completion_path)

# Important for sensitivity analysis to compare the same set of trips between runs.
# To do this need to exclude peds ID that get removed from all other runs
dfPedRoutesConsistentPeds = dfPedRoutes.loc[ ~dfPedRoutes['ID'].isin(dfPedRoutes_removedpeds['ID'])]
dfCrossEventsConsistentPeds = dfCrossEvents.loc[ ~dfCrossEvents['ID'].isin(dfPedRoutes_removedpeds['ID'])]

dfLinkCrossCounts = get_road_link_pedestrian_crossing_counts(dfCrossEventsConsistentPeds, gdfPaveLinks)


dfRouteLength = get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, exclude_stuck_peds = True, output_path = output_route_length_file)
dfSPSim = get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'dice_dist', exclude_stuck_peds = True, output_path = output_sp_similarity_path)
dfSPSimLen = get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'path_length', exclude_stuck_peds = True, output_path = output_sp_similarity_length_path)
dfConflicts = agg_cross_conflicts(dfCrossEventsConsistentPeds, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsMarked = agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['ChosenCrossingTypeString']=='unsignalised'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsUnmarked = agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['ChosenCrossingTypeString']=='unmarked'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDirect = agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='direct_cross'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonal = agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='diag_cross'], dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonalUm = agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ (dfCrossEventsConsistentPeds['linkType']=='diag_cross') & (dfCrossEventsConsistentPeds['ChosenCrossingTypeString']=='unmarked')], dfLinkCrossCounts, ttc_col = 'TTC')

######################################
#
#
# Import problem definition used for sampling
#
#
######################################
from SALib.analyze import morris, sobol
import sys
sys.path.append(".\\sample")
from SALibRepastParams import num_levels, params, random_seed, init_problem, calc_second_order
from sobol_plot import plot_sobol_indices
problem = init_problem(params = params)

rename_dict = { 'alpha':r"$\mathrm{\alpha}$",
                'lambda':r"$\mathrm{\lambda}$",
                "epsilon":r"$\mathrm{\epsilon}$",
                "gamma":r"$\mathrm{\gamma}$",
                "minCrossing": r"$\mathrm{MC}$",
                "tacticalPlanHorizon": r"$\mathrm{PH}$",
                "addVehicleTicks": r"$\mathrm{T_{veh}}$",
                "addPedTicks": r"$\mathrm{T_{ped}}$",
                "pedSpeedSeed": r"$\mathrm{Seed_{pSpeed}}$",
                "pedMassSeed": r"$\mathrm{Seed_{pMass}}$",
                "caSampleSeed": r"$\mathrm{Seed_{CA}}$",
                "vehODSeed": r"$\mathrm{Seed_{veh}}$",
                "timeThreshold": r"$\mathrm{\tau}$"
                }

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

    dfks.to_csv(output_ks_factormap)
    dfcorrs.to_csv(output_corr_factormap)

######################################
#
#
# Calculate sensitivity indices
#
#
######################################

if setting == "morris_factor_fixing":

    print("\nCalculating sensitivity indices - Route completions")

    X_rc = dfRouteCompletion.loc[:, problem['names']].values
    Y_rc = dfRouteCompletion.loc[:, 'frac_completed_journeys'].values
    Sis = morris.analyze(problem, X_rc, Y_rc, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)

    # Gather into a dataframe
    dfcompsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_compsi = morris_si_bar_figure(dfcompsi, r"Jouney Completion $\mathrm{\mu^*}$", 'Fraction journeys completed', dfcompsi['names'].replace(rename_dict))
    #f_compsi.show()
    f_compsi.savefig(os.path.join(img_dir, "route_completion_sis.{}.png".format(file_datetime_string)))


    print("\nCalculating sensitivity indices - Conflicts")

    conflicts_data = {'all':dfConflicts, 'marked':dfConflictsMarked, 'unmarked':dfConflictsUnmarked, 'direct':dfConflictsDirect, 'diag':dfConflictsDiagonal, 'diag_um':dfConflictsDiagonalUm}
    metrics = ['conflict_count', 'meanNormCC', 'varNormCC', 'meanTTC', 'varTTC']
    title_dict = {  'conflict_count':"Conflict Count", 'meanTTC':"Conflict TTC (mean)", "varTTC":"Conflict TTC (variance)",
                    'meanNormCC':'Normalised Conflict Counts (mean)', 'varNormCC': 'Normalised Conflict Counts (variance)'}
    for cat, dfC in conflicts_data.items():
        for metric in metrics:
            X = dfC.loc[:, problem['names']].values
            Y = dfC.loc[:, metric].values.astype(float)

            try:
                Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)
            except ValueError as e:
                print(e)
                print(cat)
                continue

            # Gather into a dataframe
            df = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)

            # Create figures
            f_ccsi = morris_si_bar_figure(df, r"{} - {} Crossing Sensitivity".format(title_dict[metric], cat), r"$\mathrm{\mu^*}$", df['names'].replace(rename_dict))
            #f_ccsi.show()
            f_ccsi.savefig(os.path.join(img_dir, "{}_{}_sis.{}.png".format(metric, cat, file_datetime_string)))


    print("\nCalculating sensitivity indices - Comparison to shortest path")

    # Get array of parameter values and output values
    grouped = dfSPSim.groupby("k")
    group_keys = list(grouped.groups.keys())
    for i, k in enumerate(group_keys):
        dfSPSim_k = grouped.get_group(k)
        X = dfSPSim_k.loc[:, problem['names']].values
        Y = dfSPSim_k.loc[:, 'mean'].values

        try:
            Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)
        except ValueError as e:
            print(e)
            print(k)
            continue

        # Gather into a dataframe
        dfspsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
        f_spsi = morris_si_bar_figure(dfspsi, r"Shortest Path Dice Distance Sensitivity", r"$\mathrm{\mu^*}$", dfspsi['names'].replace(rename_dict))
        #f_spsi.show()
        f_spsi.savefig(os.path.join(img_dir, "sp_similarity_sis_{}.{}.png".format(k, file_datetime_string)))

    # Repeat for shortest path similarity based on path length
    grouped = dfSPSimLen.groupby("k")
    group_keys = list(grouped.groups.keys())
    for i, k in enumerate(group_keys):
        dfSPSimLen_k = grouped.get_group(k)
        X = dfSPSimLen_k.loc[:, problem['names']].values
        Y = dfSPSimLen_k.loc[:, 'mean'].values

        try:
            Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)
        except ValueError as e:
            print(e)
            print(k)
            continue

        # Gather into a dataframe
        dfspsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
        f_spsi = morris_si_bar_figure(dfspsi, r"Path Length Sensitivity", r"$\mathrm{\mu^*}$", dfspsi['names'].replace(rename_dict))
        #f_spsi.show()
        f_spsi.savefig(os.path.join(img_dir, "sp_len_similarity_sis_{}.{}.png".format(k, file_datetime_string)))

    print("\nCalculating sensitivity indices - Total route length")

    # Get array of parameter values and output values
    X = dfRouteLength.loc[:, problem['names']].values
    Y = dfRouteLength.loc[:, 'route_length'].astype(float).values
    try:
        Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)
    except ValueError as e:
        print(e)
        print(k)

    # Gather into a dataframe
    dfRLSis = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_rlsi = morris_si_bar_figure(dfRLSis, "Total Path Length Sensitivity", r"$\mathrm{\mu^*}$", dfRLSis['names'].replace(rename_dict))
    f_rlsi.savefig(os.path.join(img_dir, "route_length_sis.{}.png".format(file_datetime_string)))

if setting == 'sobol_si':

    print("\nCalculating sobol indices - Conflicts")

    conflicts_data = {'all':dfConflicts, 'marked':dfConflictsMarked, 'unmarked':dfConflictsUnmarked, 'direct':dfConflictsDirect, 'diag':dfConflictsDiagonal, 'diag_um':dfConflictsDiagonalUm}
    metrics = ['conflict_count', 'meanNormCC']
    title_dict = {  'conflict_count':"Conflict Count", 'meanTTC':"Conflict TTC (mean)", "varTTC":"Conflict TTC (variance)",
                    'meanNormCC':'Normalised Conflict Counts (mean)', 'varNormCC': 'Normalised Conflict Counts (variance)'}
    for cat, dfC in conflicts_data.items():
        for metric in metrics:
            X = dfC.loc[:, problem['names']].values
            Y = dfC.loc[:, metric].values.astype(float)
            Sis = sobol.analyze(problem, Y, calc_second_order=calc_second_order, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)

            # Gather into a dataframe
            Sis['names'] = problem['names']
            Sis_filtered = {k:Sis[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf', 'names']}
            df = pd.DataFrame(Sis_filtered).sort_values(by='S1', ascending=False)

            # Create figures
            f_si = sobol_si_bar_figure(df, "{} Sobol Indices - {} crossings".format(title_dict[metric], cat), df['names'].replace(rename_dict))
            f_si.savefig(os.path.join(img_dir, "{}_{}_sobol1T.{}.png".format(metric, cat, file_datetime_string)))
            f_si.clear()

            if calc_second_order==True:
                f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.001, rename_dict = rename_dict)
                f_si.savefig(os.path.join(img_dir, "{}_{}_sobol2T.{}.png".format(metric, cat, file_datetime_string)))
                f_si.clear()

    print("Calculating Sobol Indices - Shortest Path Comparison")

    # Get array of parameter values and output values
    grouped = dfSPSim.groupby("k")
    group_keys = list(grouped.groups.keys())
    for i, k in enumerate(group_keys):
        dfSPSim_k = grouped.get_group(k)
        X = dfSPSim_k.loc[:, problem['names']].values
        Y = dfSPSim_k.loc[:, 'mean'].values

        Sis = sobol.analyze(problem, Y, calc_second_order=calc_second_order, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)

        # Gather into a dataframe
        Sis['names'] = problem['names']
        Sis_filtered = {k:Sis[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf', 'names']}
        df = pd.DataFrame(Sis_filtered).sort_values(by='S1', ascending=False)

        # Plot
        f_si = sobol_si_bar_figure(df, "Shortest Path Similarity Sobol Indices", df['names'].replace(rename_dict))
        f_si.savefig(os.path.join(img_dir, "sp_similarity_sobol_{}.{}.png".format(k, file_datetime_string)))
        f_si.clear()

        if calc_second_order==True:
            f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.001, rename_dict = rename_dict)
            f_si.savefig(os.path.join(img_dir, "sp_similarity_sobol_2ndorder_{}.{}.png".format(k, file_datetime_string)))
            f_si.clear()

    print("Calculating Sobol indices - Total route length")

    X = dfRouteLength.loc[:, problem['names']].values
    Y = dfRouteLength.loc[:, 'route_length'].astype(float).values
    try:
        Sis = sobol.analyze(problem, Y, calc_second_order=calc_second_order, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)
    except ValueError as e:
        print(e)
        print(k)


    # Gather into a dataframe
    Sis['names'] = problem['names']
    Sis_filtered = {k:Sis[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf', 'names']}
    df = pd.DataFrame(Sis_filtered).sort_values(by='S1', ascending=False)

    # Plot
    f_si = sobol_si_bar_figure(df, "Total Path Length Sensitivity", df['names'].replace(rename_dict))
    f_si.savefig(os.path.join(img_dir, "route_length_sobol.{}.png".format(file_datetime_string)))
    f_si.clear()

    # If second_order then produce plot show interdependence of parameter sensitivity
    if calc_second_order==True:
        f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.01, rename_dict = rename_dict)
        f_si.savefig(os.path.join(img_dir, "route_length_sobol_2ndorder_.{}.png".format(file_datetime_string)))
        f_si.clear()

#########################################
#
#
# Scatter plot of two variables, coloured by output metric
#
#
#########################################

if setting == "epsilon_gamma_scatter":

    fixed_columns = ['random_seed', 'addPedTicks', 'alpha','tacticalPlanHorizon', 'minCrossing']
    variable_columns = ['epsilon', 'gamma', 'lambda', 'addVehicleTicks']

    metric = 'frac_completed_journeys'
    groupby_columns = ['addVehicleTicks', 'lambda']
    parameter_sweep_columns = ['epsilon', 'gamma']

    output_path = os.path.join(img_dir, "fract_competed_eg.{}.png".format(file_datetime_string))
    fig_title = "Route Completions\n{} and {} parameter sweep".format(r"$\mathrm{\epsilon}$", r"$\mathrm{\gamma}$")

    rename_dict = { 'addVehicleTicks':"Ticks\nBetween\nVehicle\nAddition",
                'alpha':r"$\mathrm{\alpha}$",
                'lambda':r"$\mathrm{\lambda}$",
                "epsilon":r"$\mathrm{\epsilon}$",
                "gamma":r"$\mathrm{\gamma}$",
                0.5: r"$\mathrm{\lambda}=0.5$",
                1.5:r"$\mathrm{\lambda}=1.5$",
                5:"High\nVehicle\nFlow",
                50:"Low\nVehicle\nFlow",
                'frac_target_cross': 'Postpone crossing proportion',
                'frac_completed_journeys': 'Complete journey proportion'
                }

    f, ax = batch_run_scatter(dfRouteCompletion, groupby_columns, parameter_sweep_columns, metric, rename_dict, 'viridis', title = fig_title, cbarlabel = None, output_path = output_path)

    # Measure numbers of agents crossing at a particular link
    target_links = ['pave_link_218_219', 'pave_link_217_219']

    # Aggregate cross events to get counts of peds crossing at a particular link

    # Only consider runs where some peds peds that didn't complete route
    runs_ped_complete = dfRouteCompletion.loc[ dfRouteCompletion['frac_completed_journeys']>0, 'run'].unique()

    # Need to select based on run and ID
    dfCrossEvents['run_ID'] = dfCrossEvents.apply(lambda row: (row['run'], row['ID']), axis=1)
    dfPedRoutes_removedpeds['run_ID'] = dfPedRoutes_removedpeds.apply(lambda row: (row['run'], row['ID']), axis=1)

    dfCrossEventsCompleteJourney = dfCrossEvents.loc[ ~dfCrossEvents['run_ID'].isin(dfPedRoutes_removedpeds['run_ID'])]

    # Get count of peds per run
    dfNPeds = dfCrossEventsCompleteJourney.groupby('run')['ID'].apply(lambda s: s.unique().shape[0]).reset_index().rename(columns = {'ID':'nPedsComplete'})

    dfCrossAtTarget = dfCrossEventsCompleteJourney.loc[dfCrossEventsCompleteJourney['CurrentPavementLinkID'].isin(target_links)]
    dfCrossAtTarget = dfCrossAtTarget.groupby('run')['ID'].apply(lambda s: s.unique().shape[0]).reset_index().rename(columns = {'ID':'n_target_cross'})

    dfCrossAtTarget = pd.merge(dfRun.loc[dfRun['run'].isin(runs_ped_complete)], dfCrossAtTarget, on='run', how = 'left')
    dfCrossAtTarget = pd.merge(dfCrossAtTarget, dfNPeds, on='run', how = 'inner')

    dfCrossAtTarget['n_target_cross'] = dfCrossAtTarget['n_target_cross'].fillna(0)
    dfCrossAtTarget['frac_target_cross'] = dfCrossAtTarget['n_target_cross'] / dfCrossAtTarget['nPedsComplete']

    # Now plot
    metric = 'frac_target_cross'
    output_path = os.path.join(img_dir, "postpone_crossing_eg.{}.png".format(file_datetime_string))
    fig_title = "Postpone Crossings\n{} and {} parameter sweep".format(r"$\mathrm{\epsilon}$", r"$\mathrm{\gamma}$")

    f, ax = batch_run_tile_plot(dfCrossAtTarget, groupby_columns, parameter_sweep_columns, metric, rename_dict, plt.cm.viridis, title = fig_title, cbarlabel = None, output_path = output_path, figsize=(20,5))
