from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import networkx as nx
from scipy import stats
import itertools
import re
from shapely.geometry import LineString
from shapely.geometry import Point
from matplotlib import pyplot as plt
from math import log10, floor

######################################
#
# Functions
#
######################################
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def get_file_regex(prefix, file_datetime = None, suffix = None, ext = "csv"):
    if file_datetime is None:
        year_re = r"\d{4}"
        month_re = r"[a-zA-Z]{3}"
        day_re = r"\d{2}"
        hr_re = r"\d{2}"
        min_re = r"\d{2}"
        sec_re = r"\d{2}"
    else:
        date_string = file_datetime.strftime("%Y.%b.%d.%H.%M.%S")
        year_re, month_re, day_re, hr_re, min_re, sec_re = date_string.split(".")

    if suffix is None:
        suffix_re = r""
    else:
        suffix_re = r"\." + suffix

    file_re = re.compile(prefix+   r"\.("+ year_re + 
                                        r")\.(" + month_re + 
                                        r")\.(" + day_re    + 
                                        r")\.(" + hr_re +
                                        r")_(" + min_re + 
                                        r")_(" + sec_re + 
                                        r")" +suffix_re + 
                                        r"\." + ext + 
                                        r"")
    return file_re


def dt_from_file_name(file_name, regex):
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt

def most_recent_directory_file(directory, file_regex):
    files = os.listdir(directory)
    filtered_files = [f for f in files if file_regex.search(f) is not None]
    filtered_files.sort(key = lambda x: dt_from_file_name(x, file_regex), reverse=True)
    chosen_path = ""
    if len(filtered_files)>0:
        chosen_path = filtered_files[0]
    else:
        print("No paths found for regex {}".format(file_regex))
    return chosen_path


######################################
#
# Aggregate crossing choices by run
#
######################################

def explode_data(df, explode_col = 'edge_path'):
    data = df.values
    cols = list(df.columns)
    explode_col_index = cols.index(explode_col)
    new_data = []
    for row in data:
        for i in row[explode_col_index]:
            new_row = list(row)
            new_row[explode_col_index] = i
            new_data.append(new_row)
    dfOut = pd.DataFrame(new_data, columns = cols)
    return dfOut

def get_peds_crossing_choice(series_choices):
    series_choices = series_choices.drop_duplicates()

    if series_choices.shape[0] != 1:
        series_choices = series_choices.replace({'none':np.nan}).dropna()

    
    try:
        assert series_choices.shape[0] == 1
    except Exception as e:
        print(series_choices)
        raise AssertionError
    
    return series_choices.values[0]

def get_data_paths(file_datetime_string, data_dir):
    file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")
    file_re = get_file_regex("pedestrian_pave_link_crossings", file_datetime = file_datetime)
    ped_crossings_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex(r"^pedestrian_routes", file_datetime = file_datetime)
    ped_routes_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("vehicle_routes", file_datetime = file_datetime)
    veh_routes_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("average_vehicle_count", file_datetime = file_datetime)
    av_vehicle_counts_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("cross_events", file_datetime = file_datetime)
    cross_events_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("pedestrian_locations", file_datetime = file_datetime)
    ped_locations_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("vehicle_road_links", file_datetime = file_datetime)
    vehicle_rls_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("pedestrian_routes", file_datetime = file_datetime, suffix = 'batch_param_map')
    batch_file = most_recent_directory_file(data_dir, file_re)

    output = {}
    output["pedestrian_pave_link_crossings"]=ped_crossings_file
    output["pedestrian_routes"]=ped_routes_file
    output["vehicle_routes"]=veh_routes_file
    output['av_vehicle_counts'] = av_vehicle_counts_file
    output["cross_events"]=cross_events_file
    output["pedestrian_locations"]=ped_locations_file
    output["vehicle_road_links"]=vehicle_rls_file
    output["batch_file"]=batch_file

    return output

def get_ouput_paths(file_datetime_string, data_dir, nbins = '', ttc_threshold = 3):
    # output paths for processed data
    paths = {}
    paths["output_ped_routes_file"] = os.path.join(data_dir, "ped_routes.{}.csv".format(file_datetime_string))
    paths["output_single_ped_links_file"] = os.path.join(data_dir, "single_ped_routes.{}.csv".format(file_datetime_string))
    paths["output_vehicle_density_file"] = os.path.join(data_dir, "av_vehicle_density.{}.csv".format(file_datetime_string))
    paths["output_route_length_file"] = os.path.join(data_dir, "run_route_length.{}.csv".format(file_datetime_string))
    paths["output_sp_similarity_path"] = os.path.join(data_dir, "sp_similarity.{}.csv".format(file_datetime_string))
    paths["output_sp_similarity_length_path"] = os.path.join(data_dir, "path_length_sp_similarity.{}.csv".format(file_datetime_string))
    paths["output_route_completion_path"] = os.path.join(data_dir, "route_completions.{}.csv".format(file_datetime_string))
    paths["output_cross_events_path"] = os.path.join(data_dir, "cross_events.{}.csv".format(file_datetime_string))
    paths["output_ks_factormap"] = os.path.join(data_dir , "ks_factor_map.{}.csv".format (file_datetime_string))
    paths["output_corr_factormap"] = os.path.join(data_dir , "corr_factor_map.{}.csv".format (file_datetime_string))
    paths["output_ped_distdurs_file"] = os.path.join(data_dir, "ped_durdists.{}.csv".format(file_datetime_string))
    paths["output_veh_distdurs_file"] = os.path.join(data_dir, "veh_durdists.{}.csv".format(file_datetime_string))
    paths['output_alt_routes_file'] = os.path.join(data_dir, "alt_model_paths.{}.csv".format(file_datetime_string))
    paths["output_sd_data"] = os.path.join(data_dir, "metrics_for_sd_analysis.{}.csv".format(file_datetime_string))
    paths["output_route_data"] = os.path.join(data_dir, "metrics_for_route_analysis.{}.csv".format(file_datetime_string))
    paths["output_cross_entropy"] = os.path.join(data_dir, "cross_loc_entropy.{}bins.{}.csv".format(nbins, file_datetime_string))
    paths["output_cross_conflicts"] = os.path.join(data_dir, "conflicts.{}.{}.csv".format(ttc_threshold, file_datetime_string))

    return paths

def crossing_percentages(row, c1 = 'unmarked', c2 = 'unsignalised', scale = 100):
    '''Calculates percentages that cross at either crossing as proportion of those that do crossing the road, not those that are undecided.
    '''

    crossing_total = row[[c1,c2]].sum()

    row[c1+'_pcnt'] = row[c1] / crossing_total * scale
    row[c2+'_pcnt'] = row[c2] / crossing_total * scale

    return row

def load_batch_data(data_dir, file_prefix, file_datetime = None, run_selection_dict = {}):

    file_re = get_file_regex(file_prefix, file_datetime = file_datetime)
    data_file = most_recent_directory_file(data_dir, file_re)

    df_data = pd.read_csv(os.path.join(data_dir, data_file))

    file_re = get_file_regex(file_prefix, file_datetime = file_datetime, suffix = 'batch_param_map')
    
    batch_file = most_recent_directory_file(data_dir, file_re)
    df_run = pd.read_csv(os.path.join(data_dir, batch_file))

    for col in run_selection_dict.keys():
        df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]

    return df_data, df_run

def adjacent_edge_node(n, e):
    if n == e[0]:
        return e[1]
    elif n == e[1]:
        return e[0]
    else:
        return None

def node_path_from_edge_path_old(edge_id_path, start_node, end_node, pavement_graph):
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

        if (next_candidate is not None):
             node_path.append(candidate)
             prev = candidate

    # If ped agent got removed from the simulation it would not have reached the end node and the path will not have any other nodes added
    # In this case return empty list
    if len(node_path)==1:
        return []

    return node_path[::-1]

def get_strategic_path_pavement_edges(strategic_path, gdfORLinks, gdfPaveNodes, pavement_graph):
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

    filter_pavement_edges = get_strategic_path_pavement_edges(strategic_path, gdfORLinks, gdfPaveNodes, pavement_graph)

    sub_pavement_graph = nx.edge_subgraph(pavement_graph, filter_pavement_edges)

    try:
        dijkstra_path = nx.dijkstra_path(sub_pavement_graph, start_node, end_node, weight = weight)
    except nx.exception.NodeNotFound as e:
        print(start_node, end_node, strategic_path)
        return None

    return tuple(dijkstra_path)

def simple_paths_within_strategic_path(strategic_path, gdfORLinks, gdfPaveNodes, pavement_graph, start_node, end_node):

    filter_pavement_edges = get_strategic_path_pavement_edges(strategic_path, gdfORLinks, gdfPaveNodes, pavement_graph)

    sub_pavement_graph = nx.edge_subgraph(pavement_graph, filter_pavement_edges)

    try:
        simple_paths = nx.all_simple_paths(sub_pavement_graph, start_node, end_node)
    except nx.exception.NodeNotFound as e:
        print(start_node, end_node, strategic_path)
        return None

    return tuple(simple_paths)

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
    print("\nCalculating Road Link Aggregated Crossing Counts")

    dfCrossCounts = dfCrossEvents.groupby(['run','TacticalEdgeID'])['ID'].apply(lambda s: s.unique().shape[0]).reset_index()
    dfCrossCounts.rename(columns={'ID':'cross_count'}, inplace=True)

    # Now merge with lookup from or link to pave link
    dfCrossCounts = pd.merge(dfCrossCounts, gdfPaveLinks.reindex(columns = ['fid', 'pedRLID']).drop_duplicates(), left_on = 'TacticalEdgeID', right_on = 'fid')

    # Aggregate to get road link cross counts
    dfRLCrossCounts = dfCrossCounts.groupby(['run','pedRLID'])['cross_count'].apply(lambda s: s.sum()).reset_index()

    # Now merge this to lookup from OR link to pavement link
    dfRLCrossCounts = pd.merge(dfRLCrossCounts, gdfPaveLinks.reindex(columns = ['fid', 'pedRLID']).drop_duplicates(), on = 'pedRLID')

    return dfRLCrossCounts

def get_shortest_path_similarity(dfPedRoutes, dfRun, pavement_graph, dict_node_pos, weight_params, distance_function = 'dice_dist', output_path = "sp_similarity.csv"):

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
            dfPedRoutes['comp_value_{}'.format(k)] = dfPedRoutes.apply(lambda row: compare_node_paths(pavement_graph, row['node_path'], row['sp_{}'.format(k)], dict_node_pos, distance_function = distance_function, weight='length'), axis=1)

            # This is only meaning full for the shortest path unweighted by vehicle traffic, where we can expect the path to match the ABM tactical path, and therefore compare path lengths to check for equivalence.
            dfPedRoutes['comp_path_weight_{}'.format(k)] = dfPedRoutes.apply(lambda row: compare_node_paths(pavement_graph, row['node_path'], row['sp_{}'.format(k)], dict_node_pos, distance_function = distance_function, account_for_path_length=True, weight='length'), axis=1)

            df = dfPedRoutes.groupby('run')['comp_value_{}'.format(k)].describe().reset_index()
            df['k'] = k

            # Also get count of 0 difference to sp
            dfQ1Counts = dfPedRoutes.groupby('run')['comp_value_{}'.format(k)].apply(lambda s: (s<0.00000001).value_counts()[True]).reset_index().rename(columns = {'comp_value_{}'.format(k):'cv{}zeroCount'.format(k)})
            df = pd.merge(df, dfQ1Counts, on = 'run', how = 'left')

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

def get_run_total_route_length(dfPedRoutes, dfRun, pavement_graph, output_path = "run_route_length.csv"):

    ######################################
    #
    #
    # Compare these various shortest path routes to the ABM routes by calculating a metric of difference between their node paths
    #
    #
    ######################################

    if os.path.exists(output_path)==False:
        dfPedRoutes['route_length'] = dfPedRoutes.apply(lambda row: nx.path_weight(pavement_graph, row['node_path'], weight='length'), axis=1)

        dfRouteLength = dfPedRoutes.groupby('run')['route_length'].sum().reset_index()

        # Also calculate route length per ped
        dfRouteLengthPerPed = dfPedRoutes.groupby('run')['route_length'].mean().reset_index().rename(columns = {'route_length':'route_length_pp'})

        # Merge in parameter values
        dfRouteLength = pd.merge(dfRun, dfRouteLength, on = 'run')
        dfRouteLength = pd.merge(dfRouteLength, dfRouteLengthPerPed, on = 'run')

        # Save date for future use
        dfRouteLength.to_csv(output_path, index=False)
    else:
        dfRouteLength = pd.read_csv(output_path)

    return dfRouteLength

def agg_route_completions(dfPedRoutes, dfRun, output_path = 'route_completions.csv'):
    if os.path.exists(output_path)==False:
        print("\nProcessing Pedestrian Route Completions")
        dfPedRoutes['completed_journey'] = dfPedRoutes['node_path'].map(lambda x: int(len(x)>0))
        dfCompletions = dfPedRoutes.groupby('run').apply( lambda df: df['completed_journey'].sum() / float(df.shape[0])).reset_index().rename(columns = {0:'frac_completed_journeys'})
        dfCompletions = pd.merge(dfRun, dfCompletions, on = 'run')
        dfCompletions.to_csv(output_path, )
    else:
        print("\nLoading Pedestrian Route Completions")
        dfCompletions = pd.read_csv(output_path)

    return dfCompletions

def agg_cross_conflicts(dfCrossEvents, dfRun, dfLinkCrossCounts, output_path, ttc_col = 'TTC', ttc_threshold = 3):
    '''Aggregate crossing events to create indicators of conflict for each run. This involves findings the total number of conflicts per run and the
    mean TTC per run.
    '''
    if os.path.exist(output_path):
        print("\nLoading crossing conflicts, {}, {}".format(ttc_col, ttc_threshold))
        dfRunConflicts = pd.read_csv(output_path)

    else:
        print("\nProcessing crossing conflicts, {}, {}".format(ttc_col, ttc_threshold))
        
        # Join with pavement links data and aggregate to OR road link ID to get number of peds per OR road link for normalising crossing counts
        calc_conflict_count = lambda s: s.dropna().loc[s<ttc_threshold].shape[0]
        calc_mean_ttc = lambda s: s.dropna().loc[s<ttc_threshold].mean()
        calc_var_ttc = lambda s: s.dropna().loc[s<ttc_threshold].var()

        dfRunConflicts = dfCrossEvents.groupby("run").agg(  conflict_count=pd.NamedAgg(column=ttc_col, aggfunc=calc_conflict_count),
                                                            meanTTC=pd.NamedAgg(column=ttc_col, aggfunc=calc_mean_ttc),
                                                            varTTC=pd.NamedAgg(column=ttc_col, aggfunc=calc_var_ttc),
                                                            ).reset_index()

        # get conflict counts normalsied by numbers of peds on road links
        dfPaveLinkConflictCounts = dfCrossEvents.groupby(['run','TacticalEdgeID']).agg( conflict_count=pd.NamedAgg(column=ttc_col, aggfunc=calc_conflict_count),).reset_index()
        dfPaveLinkConflictCounts = pd.merge(dfPaveLinkConflictCounts, dfLinkCrossCounts, left_on = ['run', 'TacticalEdgeID'], right_on = ['run', 'fid'], indicator=True)
        assert dfPaveLinkConflictCounts.loc[ dfPaveLinkConflictCounts['_merge']!='both'].shape[0]==0
        dfPaveLinkConflictCounts['norm_conflict_count'] = dfPaveLinkConflictCounts['conflict_count'] / dfPaveLinkConflictCounts['cross_count']

        # now aggregate trun
        dfRunConflictsNorm = dfPaveLinkConflictCounts.groupby("run").agg(   meanNormCC=pd.NamedAgg(column='norm_conflict_count', aggfunc=np.mean),
                                                                            varNormCC =pd.NamedAgg(column='norm_conflict_count', aggfunc=np.var),).reset_index()

        dfRunConflicts = pd.merge(dfRunConflicts, dfRunConflictsNorm, on='run')

        # Merge with run params
        dfRunConflicts = pd.merge(dfRun, dfRunConflicts, on='run')

    return dfRunConflicts

def get_processed_crossing_locations_data(data_dir, file_prefix, file_datetime = None):

    df_data, df_run = load_batch_data(data_dir, file_prefix, file_datetime = file_datetime)

    # Select just the crossing choice column
    desired_columns = ['ID', 'CrossingChoice','tick','run']
    df_data = df_data.reindex(columns = desired_columns)


    df_ped_cc = df_data.groupby(['run','ID'], group_keys=True)['CrossingChoice'].apply(get_peds_crossing_choice).reset_index()

    # Get count of peds each run
    ser_run_ped_counts = df_ped_cc.groupby('run')['ID'].apply(lambda s: s.unique().shape[0])
    df_ped_counts = pd.DataFrame({'run_npeds':ser_run_ped_counts})

    # Now group by run to get count of each crossing type
    df_cc_count = df_ped_cc.groupby(['run','CrossingChoice']).count().unstack()
    df_cc_count.columns = [c[1] for c in df_cc_count.columns]

    if 'none' in df_cc_count.columns:
        df_cc_count.rename(columns = {'none':'undecided'}, inplace=True)
    else:
        df_cc_count['undecided'] = np.nan
    df_cc_count.fillna(0, inplace = True)

    # Join to df of npeds per run and calculate percentages
    df_cc_count = pd.merge(df_cc_count, df_ped_counts, left_index = True, right_index = True)

    df_cc_count = df_cc_count.apply(crossing_percentages, axis = 1)

    df_cc_count['undecided_frac'] = df_cc_count['undecided'] / df_cc_count['run_npeds']

    df_cc_count = pd.merge(df_cc_count, df_run, left_index = True, right_on = 'run', how = 'inner')

    return df_cc_count

def get_pedestrian_run_durations(dfCrossEvents):
    # Get duration for each run - defined as total time pedestrians are in the simulation.
    dfPedStart = dfCrossEvents.groupby('run')['tick'].min().reset_index()
    dfPedEnd = dfCrossEvents.groupby('run')['tick'].max().reset_index()
    dfDurs = pd.merge(dfPedStart, dfPedEnd, on = 'run', suffixes = ('_start', '_end'))
    dfDurs['duration'] = dfDurs['tick_end'] - dfDurs['tick_start']
    return dfDurs

def get_road_link_vehicle_density_from_vehicle_routes(gdfITNLinks, data_file, output_path):

    if os.path.exists(output_path) == False:
        
        dfVehRoutes = pd.read_csv(data_file)

        # Convert strategic path to list
        dfVehRoutes['FullStrategicPathString'] = dfVehRoutes['FullStrategicPathString'].map(lambda s: tuple(dict.fromkeys(s.strip(":").split(":"))))

        # Then explode list of road links in order to get total number of vehicle per road link
        dfRunRLs = explode_data(dfVehRoutes.reindex(columns = ['run','FullStrategicPathString']), explode_col='FullStrategicPathString')
        dfVehCounts = dfRunRLs.groupby(['run','FullStrategicPathString']).apply(lambda df: df.shape[0]).reset_index()
        dfVehCounts.rename(columns = {0:'VehCount', 'FullStrategicPathString':'fid'}, inplace=True)

        # Merge in tick, which is measure of duration of run
        dfVehCounts = pd.merge(dfVehCounts, dfVehRoutes.reindex(columns=['run','tick']).drop_duplicates(), on='run')

        dfVehCounts['AvVehCount'] = dfVehCounts['VehCount'] / dfVehCounts['tick']

        # Merge with ITN links to get lookup to ped rl ID
        gdfITNLinks = gdfITNLinks.reindex(columns = ['fid','length', 'pedRLID'])

        dfVehRoutes = None

        dfVehCounts = pd.merge(dfVehCounts, gdfITNLinks, left_on = 'fid', right_on = 'fid', how = 'left')
        dfVehCounts['AvVehDen'] = dfVehCounts['AvVehCount'] / dfVehCounts['length']
        dfVehCounts.drop('length', axis=1, inplace=True)

        dfVehCounts.to_csv(output_path, index=False)
    else:
        dfVehCounts = pd.read_csv(output_path)

    return dfVehCounts

def get_road_link_vehicle_density_from_vehicle_counts(gdfITNLinks, data_file, output_path):

    if os.path.exists(output_path) == False:
        
        dfVehCounts = pd.read_csv(data_file)
        dfVehCounts.rename(columns = {'FID':'fid'}, inplace=True)

        # Merge with ITN links to get lookup to ped rl ID
        if 'length' not in gdfITNLinks.columns:
            gdfITNLinks['length'] = gdfITNLinks.geometry.length
        gdfITNLinks = gdfITNLinks.reindex(columns = ['fid','length', 'pedRLID'])

        dfVehCounts = pd.merge(dfVehCounts, gdfITNLinks, left_on = 'fid', right_on = 'fid', how = 'left')
        dfVehCounts['AvVehDen'] = dfVehCounts['AvVehCount'] / dfVehCounts['length']
        dfVehCounts.drop('length', axis=1, inplace=True)

        dfVehCounts.to_csv(output_path, index=False)
    else:
        dfVehCounts = pd.read_csv(output_path)

    return dfVehCounts

def edge_path_from_tactical_route_string_old(trs):
    ep = tuple(dict.fromkeys(trs.strip(":").split(":")))
    return ep

def edge_path_from_tactical_route_string(trs):

    # First get unique set of links in order of when they first appear in the path
    links = tuple(dict.fromkeys(trs.strip(":").split(":")))

    # Need to loop through the path anre remove duplicated links in a row.
    # The whole path can have duplicated links but not next to each other
    ep = []
    orig_ep = trs.strip(":").split(":")
    i=0
    ep.append(orig_ep[i])
    n=0
    while (n<len(orig_ep)):
        # if next link is different, add it to path and set counter to zero
        if (orig_ep[i]==orig_ep[n]):
            pass
        else:
            i=n
            ep.append(orig_ep[i])
        n+=1

    # Finally, identify whether the agent doubles back at any point in the path. Duplicate the double back link where this occurs
    ep_final = []
    prev_direction = True
    for i, e in enumerate(ep[:-1]):
        ep_final.append(e)

        e_next = ep[i+1]
        direction = links.index(e)<links.index(e_next) # true means ped moving to unseen link 

        if (prev_direction & ~direction): 
            # means that a double back has occurred, ped was moving to unseen links then moved to already traversed link
            ep_final.append(e)

        prev_direction = direction

    ep_final.append(ep[-1])

    return tuple(ep_final)

def load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, dfVehCounts = None, ped_routes_path = "ped_routes.csv", strategic_path_filter = True):

    d,f = os.path.split(ped_routes_path)
    output_path = os.path.join(d, 'processed_'+f)

    split_path = os.path.splitext(ped_routes_path)
    routes_removed_path = split_path[0] + "_removed_peds" + split_path[1]

    if os.path.exists(output_path):
        print("\nLoading Pedestrian Agent Routes")
        dfPedRoutes = pd.read_csv(output_path)
        dfPedRoutes_removedpeds = pd.read_csv(routes_removed_path)

        # Convert columns to tuples
        dfPedRoutes['FullStrategicPathString'] = dfPedRoutes['FullStrategicPathString'].map(lambda s: tuple(s.strip("('").strip("')").strip("',").split("', '")))
        dfPedRoutes['edge_path'] = dfPedRoutes['edge_path'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))
        dfPedRoutes['node_path'] = dfPedRoutes['node_path'].map(lambda s: s.strip("['").strip("']").split("', '"))

        for k in weight_params:
            dfPedRoutes['sp_{}'.format(k)] = dfPedRoutes['sp_{}'.format(k)].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))
    else:
        # Otherwise create data

        print("\nProcessing Pedestrian Agent Routes")


        # Load data from file
        dfPedRoutes = pd.read_csv(ped_routes_path)

        # Convert strategic path to list
        dfPedRoutes['FullStrategicPathString'] = dfPedRoutes['FullStrategicPathString'].map(lambda s: tuple(dict.fromkeys(s.strip(":").split(":"))))

        # Convert string tactical path to list
        dfPedRoutes['edge_path'] = dfPedRoutes['FullTacticalRouteString'].map(lambda s: edge_path_from_tactical_route_string(s))

        # Convert edge path to node path
        dfPedRoutes['node_path'] = dfPedRoutes.apply(lambda row: node_path_from_edge_path(row['edge_path'], row['StartPavementJunctionID'], row['DestPavementJunctionID'], pavement_graph), axis=1)

        ## Need to keep these in or keep a record of what rows got dropped in order to calculate SIs later
        # Or avoid dropping.
        dfPedRoutes_removedpeds = dfPedRoutes.loc[ dfPedRoutes['node_path'].map(lambda x: len(x))==0]
        removed_peds_index = dfPedRoutes_removedpeds.index

        # Find the unique set of start and end node and calculate shortest paths between these. Then merge into the ped routes data.
        dfPedRoutes['start_node'] = dfPedRoutes.apply(lambda row: row['node_path'][0] if len(row['node_path'])>0 else row['StartPavementJunctionID'], axis=1)
        dfPedRoutes['end_node'] = dfPedRoutes.apply(lambda row: row['node_path'][-1] if len(row['node_path'])>0 else row['DestPavementJunctionID'], axis=1)
        dfUniqueStartEnd = dfPedRoutes.loc[:, ['start_node', 'end_node', 'FullStrategicPathString']].drop_duplicates()

                # Set dfPedRoutes to just have columns we are interested in
        dfPedRoutes = dfPedRoutes.reindex(columns = [   'run', 'ID', 'FullStrategicPathString', 'edge_path',
                                                        'StartPavementJunctionID', 'DestPavementJunctionID', 'node_path',
                                                        'start_node', 'end_node', 'PostponeCounts'])

        if (dfVehCounts is not None) & (weight_params is not None):
            print("Calculating SPs with cross costs")
            dfAltPathsOrig = alt_paths_cross_weight(dfPedRoutes, dfUniqueStartEnd, gdfPaveLinks, gdfORLinks, gdfPaveNodes, dfVehCounts, pavement_graph, weight_params, strategic_path_filter)
        else:
            print("Calculating SPs without cross costs")
            dfAltPaths = alt_paths(dfPedRoutes, dfUniqueStartEnd, gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, strategic_path_filter)
            

        dfPedRoutes = pd.merge(dfPedRoutes, dfAltPaths, on = ['run', 'start_node', 'end_node', 'FullStrategicPathString'])

        dfPedRoutes.to_csv(output_path, index=False)
        dfPedRoutes_removedpeds.to_csv(routes_removed_path, index=False)

    return dfPedRoutes, dfPedRoutes_removedpeds

def alt_paths_cross_weight(dfPedRoutes, dfUniqueStartEnd, gdfPaveLinks, gdfORLinks, gdfPaveNodes, dfVehCounts, pavement_graph, weight_params, strategic_path_filter):
    dfAltPaths = pd.DataFrame()
    for run in dfPedRoutes['run'].unique():
        dfAltPathsRun = dfUniqueStartEnd.copy()
        dfAltPathsRun['run'] = run
        for k in weight_params:
            weight_name = "weight{}".format(k)

            # Create df of costs to assign to pavement network links based on length and vehicle flow
            dfLinkWeights = gdfPaveLinks.reindex(columns = ['fid', 'MNodeFID', 'PNodeFID', 'pedRLID', 'length'])
               
            dfRunVehCounts = dfVehCounts.loc[dfVehCounts['run']==run]
            dfLinkWeights = pd.merge(dfLinkWeights, dfRunVehCounts, on='pedRLID', how = 'left')
            dfLinkWeights['cross_cost'] = dfLinkWeights['AvVehDen'].fillna(0) * k

            dfLinkWeights[weight_name] = dfLinkWeights['length'] + dfLinkWeights['cross_cost']

            # Weight pavement network crossing links by average vehicle flow
            weight_attributes = dfLinkWeights.set_index( dfLinkWeights.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
            nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

            if strategic_path_filter:
                dfAltPathsRun['sp_{}'.format(k)] = dfAltPathsRun.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'], weight = weight_name), axis=1)
            else:
                dfAltPathsRun['sp_{}'.format(k)] = dfAltPathsRun.apply(lambda row: nx.dijkstra_path(pavement_graph, row['start_node'], row['end_node'], weight =  weight_name), axis=1)

        dfAltPaths = pd.concat([dfAltPaths, dfAltPathsRun])
    return dfAltPaths

def alt_paths(dfPedRoutes, dfUniqueStartEnd, gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, strategic_path_filter):
    
    dfAltPathsBase = dfUniqueStartEnd.copy()
    weight_name = "weight{}".format(0)

    # Create df of costs to assign to pavement network links based on length and vehicle flow
    dfLinkWeights = gdfPaveLinks.reindex(columns = ['fid', 'MNodeFID', 'PNodeFID', 'pedRLID', 'length'])
    dfLinkWeights['cross_cost'] = 0
    dfLinkWeights[weight_name] = dfLinkWeights['length'] + dfLinkWeights['cross_cost']

    # Weight pavement network crossing links by average vehicle flow
    weight_attributes = dfLinkWeights.set_index( dfLinkWeights.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
    nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

    if strategic_path_filter:
        dfAltPathsBase['sp_{}'.format(0)] = dfAltPathsBase.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'], weight = weight_name), axis=1)
    else:
        dfAltPathsBase['sp_{}'.format(0)] = dfAltPathsBase.apply(lambda row: nx.dijkstra_path(pavement_graph, row['start_node'], row['end_node'], weight =  weight_name), axis=1)

    dfAltPaths = pd.DataFrame()
    for run in dfPedRoutes['run'].unique():
        dfAltPathsRun = dfAltPathsBase.copy()
        dfAltPathsRun['run'] = run          
        dfAltPaths = pd.concat([dfAltPaths, dfAltPathsRun])
    return dfAltPaths


def median_ped_pavement_link_counts(dfPedRoutes, start_node = None, output_path = 'single_ped_links.csv'):
    '''Selects all of the tactical paths traverse by one pedestrian agetn across all simulation runs. From this calculates number of times each link is 
    traversed. Used to visualise path heterogeneity.
    '''

    # Calculate number of links in each strategic path and choose a pedestrian with a median length path
    if start_node is None:
        dfPedRoutes['sp_len'] = dfPedRoutes['FullStrategicPathString'].map(lambda x: len(x))
        pedID = dfPedRoutes.loc[ dfPedRoutes['sp_len'] == dfPedRoutes['sp_len'].median(), 'ID'].values[0]
    else:
        pedID = dfPedRoutes.loc[ dfPedRoutes['start_node'] == start_node, 'ID'].values[0]
    
    # Now get all ped routes for this ped so they can be visualised
    dfSinglePedRoutes = dfPedRoutes.loc[ dfPedRoutes['ID'] == pedID]
    assert dfSinglePedRoutes['run'].unique().shape[0] == dfPedRoutes['run'].unique().shape[0]
    
    # Select only columns of interest and explode out the pat links
    dfSinglePedRoutes = dfSinglePedRoutes.reindex(columns = ['ID','run','FullStrategicPathString','start_node','end_node', 'edge_path'])
    dfSinglePedRoutes = explode_data(dfSinglePedRoutes, explode_col='edge_path')

    # Save this subset and visualise
    dfSinglePedRoutes.to_csv(output_path, index=False)

    return dfSinglePedRoutes, pedID

def load_sp_model_shortest_paths(dfPedRoutes, dfRun, gdfORLinks, gdfPaveLinks, gdfPaveNodes, gdfCAs, pavement_graph, weight_params, dfVehCounts, alt_routes_path, strategic_path_filter = True):
    '''
    Method for calculating routes accoring to an alternative shortest path model. Useful for comparing my hierarchical model to
    '''
    p, ext = os.path.splitext(alt_routes_path)
    wmax = max(list(weight_params))
    wstep = list(weight_params)[1] - list(weight_params)[0]
    output_path = p.split(".")[0]+"_{}_{}_{}".format(wmax, wstep, strategic_path_filter)+".".join(p.split(".")[1:])+ext

    if os.path.exists(output_path)==True:
        print("\nLoading existing alternative model paths")

        dfAltPaths = pd.read_csv(output_path)

        # Reformal paths
        dfAltPaths['sp'] = dfAltPaths['sp'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))
        dfAltPaths['alt_path'] = dfAltPaths['alt_path'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")) if pd.isnull(s)==False else s)
    else:
        print("\nCalculating alternative model paths")

        dfUniqueStartEnd = dfPedRoutes.loc[:, ['start_node', 'end_node', 'FullStrategicPathString']].drop_duplicates()

        # To do that ideantify direct crossings that correspond to crossing infrastructure
        gdfPaveLinksDirectCrossings = gdfPaveLinks.loc[ gdfPaveLinks['linkType'] == 'direct_cross']
        gdfPaveLinksCAs = gpd.sjoin(gdfPaveLinks.loc[ gdfPaveLinks['linkType'] == 'direct_cross'], gdfCAs, op = 'within')
        direct_crossing_with_marked_cas = gdfPaveLinksCAs['fid'].unique()

        # Create alternative set of paths for each vehicle flow setting
        data = {'run':[], 'start_node':[], 'end_node':[], 'sp':[], 'k':[], 'j':[], 'ratio_min':[], 'ratio_max':[], 'ratio_median':[], 'alt_path':[], 'alt_path_length':[]}
        for run in dfRun.drop_duplicates(subset = 'addVehicleTicks')['run'].unique():
            if (int(run) % 10) == 0:
                print("Run:{}".format(run))
            for k in weight_params:
                for j in weight_params:
                    weight_name = "weight{}_{}".format(k, j)

                    # Calculate pavement link weights based on vehicle density
                    dfLinkWeights = gdfPaveLinks.reindex(columns = ['fid', 'MNodeFID', 'PNodeFID', 'pedRLID', 'length'])

                    if dfVehCounts is None:
                        dfLinkWeights['cross_cost'] = 0
                    else:
                        dfRunVehCounts = dfVehCounts.loc[dfVehCounts['run']==run]
                        dfLinkWeights = pd.merge(dfLinkWeights, dfRunVehCounts, on='pedRLID', how = 'left')

                        # Initialise cross cost as zero. Calculate mean veh den to account for one:many lookup from or to itn links
                        dfLinkWeights['AvVehDen'] = dfLinkWeights['AvVehDen'].fillna(0)
                        dfLinkWeights = dfLinkWeights.groupby(['fid', 'MNodeFID', 'PNodeFID', 'pedRLID', 'length'], dropna=False)['AvVehDen'].mean().reset_index() # Very important to set dropna=False otherwise non-crossing pavement links get dropped from the dataset.

                        # Then for road crossing link set crossing cost as a multiple of the average vehicle desnity on the link the crossing is on.
                        dfLinkWeights['cross_cost']=0
                        dfLinkWeights.loc[ dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'cross_cost'] = dfLinkWeights.loc[ dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'AvVehDen'] * k
                        dfLinkWeights.loc[ ~dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'cross_cost'] = dfLinkWeights.loc[ ~dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'AvVehDen'] * j
                        dfLinkWeights[weight_name] = dfLinkWeights['length'] + dfLinkWeights['cross_cost']

                    # sense check the weights by recording the median, min and max vaues of cross cost to length ratio for crossing links only
                    ratio = (dfLinkWeights.loc[~dfLinkWeights['pedRLID'].isnull(), 'cross_cost'] / dfLinkWeights.loc[~dfLinkWeights['pedRLID'].isnull(), 'length'])
                    min_ratio = ratio.min()
                    max_ratio = ratio.max()
                    median_ratio = ratio.median()

                    # Weight pavement network crossing links by average vehicle flow
                    weight_attributes = dfLinkWeights.set_index( dfLinkWeights.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
                    nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

                    # With network weights set can now calculate paths
                    for start_node, end_node, sp in dfUniqueStartEnd.values:
                        try:
                            if strategic_path_filter:
                               alt_model_path = shortest_path_within_strategic_path(sp, gdfORLinks, gdfPaveNodes, pavement_graph, start_node, end_node, weight = weight_name)
                               alt_path_length = nx.path_weight(pavement_graph, alt_model_path, 'length')
                            else:
                                alt_model_path = nx.dijkstra_path(pavement_graph, start_node, end_node, weight =  weight_name)
                                alt_path_length = nx.path_weight(pavement_graph, alt_model_path, 'length')
                        except Exception as err:
                            print(start_node, end_node, sp, run, k, j)
                            alt_model_path = None
                            alt_path_length = None

                        data['run'].append(run)
                        data['start_node'].append(start_node)
                        data['end_node'].append(end_node)
                        data['sp'].append(sp)
                        data['j'].append(j)
                        data['k'].append(k)
                        data['ratio_min'].append(min_ratio)
                        data['ratio_max'].append(max_ratio)
                        data['ratio_median'].append(median_ratio)
                        data['alt_path'].append(alt_model_path)
                        data['alt_path_length'].append(alt_path_length)

        dfAltPaths = pd.DataFrame(data)
        dfAltPaths.to_csv(output_path, index=False)

    return dfAltPaths

def load_and_clean_cross_events(gdfPaveLinks, cross_events_path = "cross_events.csv"):
    '''Method to aggregate crossing events from the Ped Crossings dataset, to produce dataframe with single row per crossing event.
    Crossing event defined by crossing type, location and minimum TTC during crossing.
    '''
    d,f = os.path.split(cross_events_path)
    output_path = os.path.join(d, 'processed_'+f)

    if os.path.exists(output_path)==False:

        print("\nProcessing Pedestrian Cross Events")

        dfCrossEvents = pd.read_csv(cross_events_path)

        dfCrossEvents = dfCrossEvents.reindex(columns = ['run', 'ID', 'FullStrategicPathString', 'CrossingType', 'TacticalEdgeID', 'CrossingCoordinatesString', 'tick', 'TTC'])

        # Group by run, ID and TacticalEdgeID to find min TTC per cross event
        dfCrossEvents['TTC'] = dfCrossEvents.groupby(['run', 'ID', 'CrossingType', 'TacticalEdgeID'])['TTC'].transform(lambda s: s.min())
        dfCrossEvents['tick_start'] = dfCrossEvents.groupby(['run', 'ID', 'CrossingType', 'TacticalEdgeID'])['tick'].transform(lambda s: s.min())
        dfCrossEvents['tick_end'] = dfCrossEvents.groupby(['run', 'ID', 'CrossingType', 'TacticalEdgeID'])['tick'].transform(lambda s: s.max())
        dfCrossEvents.drop('tick',axis=1, inplace=True)

        # Drop duplicates again now that TTC and crossing coord processed
        dfCrossEvents = dfCrossEvents.drop_duplicates(subset=['run', 'ID', 'CrossingType', 'TacticalEdgeID'])

        # Merge with pave links to get link type
        dfLinkTypes = gdfPaveLinks.reindex(columns = ['fid','linkType']).drop_duplicates()
        dfCrossEvents = pd.merge(dfCrossEvents, dfLinkTypes, left_on = 'TacticalEdgeID', right_on = 'fid', how = 'left')
        dfCrossEvents.drop('fid', axis=1,inplace=True)

        # Check that there is a single crossing event per ped per pavement link.
        cross_per_ped_link = dfCrossEvents.groupby(['run', 'ID', 'TacticalEdgeID']).apply(lambda df: df.shape[0])
        nMultiCross = cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]
        if nMultiCross!=0:
            print("WARNING: {} road links with multiple crossings per ped".format(nMultiCross))

        dfCrossEvents.to_csv(output_path, index=False)
    else:
        print("\nLoading Existing Processed Pedestrian Cross Events")
        dfCrossEvents = pd.read_csv(output_path)

    return dfCrossEvents

def linestring_from_crossing_coord_string(ccs, coord_regex = re.compile(r"(\d+.\d+)")):
    xys = coord_regex.findall(ccs)
    p1 = Point(*map(float, xys[:2]))
    p2 = Point(*map(float, xys[2:]))
    l = LineString([p1,p2])
    return l

def road_link_crossing_point(cross_line_geom, rl_geom):

    intersection = cross_line_geom.intersection(rl_geom)

    c_point=None
    if intersection.is_empty:
        midpoint = cross_line_geom.centroid
        d1 = midpoint.distance(Point(rl_geom.coords[0]))
        d2 = midpoint.distance(Point(rl_geom.coords[-1]))

        if d1<d2:
            c_point = rl_geom.coords[0]
        else:
            c_point = rl_geom.coords[-1]
    else:
        t = intersection.type
        if t=='Point':
            c_point = intersection.coords[0]
        else:
            c_point = intersection.centroid.coords[0]

    return Point(c_point)
def crossing_location(cross_point, rl_geom, rl_fid, or_first_node, gdfORLinks):
    road_link_start_coord_index = None
    if or_first_node==gdfORLinks.loc[gdfORLinks['fid']==rl_fid, 'MNodeFID'].values[0]:
        road_link_start_coord_index = 0
    else:
        road_link_start_coord_index = -1

    # Calculate distance of crossing point from start of road link
    d = cross_point.distance(Point(rl_geom.coords[road_link_start_coord_index]))

    return d


def crossing_location_quantile_bin(cross_point, rl_geom, rl_fid, or_first_node, gdfORLinks, nbins):
    d = crossing_location(cross_point, rl_geom, rl_fid, or_first_node, gdfORLinks)

    # Now bin quantile
    cross_bin = int(np.floor( (d/rl_geom.length) * nbins))
    if cross_bin>=nbins:
        cross_bin = nbins-1

    return cross_bin

def crossing_location_distance_bin(cross_point, rl_geom, rl_fid, or_first_node, gdfORLinks, bin_dist=2):
    d = crossing_location(cross_point, rl_geom, rl_fid, or_first_node, gdfORLinks)

    # Now bin quantile
    cross_bin = d//bin_dist

    return cross_bin

def get_crossing_locations_and_bins(dfCrossEvents, dfPedPaths, gdfPaveLinks, gdfPaveNodes, gdfORLinks, nbins, bin_dist):
    '''
    '''
    # Merge in ped paths so that the direction they walk along the road link being crossed can be determined
    dfCrossEvents = pd.merge(dfCrossEvents, dfPedPaths, on = ['run','ID'], how='left', indicator = True)
    assert dfCrossEvents.loc[ dfCrossEvents['_merge']!='both'].shape[0]==0
    dfCrossEvents.drop('_merge', axis=1, inplace=True)

    # Merge to Pave Links and to OR links to get length of road being crossed.
    dfCrossEvents = pd.merge(dfCrossEvents,  gdfPaveLinks.reindex(columns = ['pedRLID', 'fid', 'MNodeFID','PNodeFID']), left_on = 'TacticalEdgeID', right_on = 'fid', how = 'left', indicator=True)
    assert dfCrossEvents.loc[ dfCrossEvents['_merge']!='both'].shape[0]==0
    dfCrossEvents.drop('_merge', axis=1, inplace=True)

    # Merge to OR Links to get road length
    dfCrossEvents = pd.merge(dfCrossEvents,  gdfORLinks.reindex(columns = ['fid', 'geometry']), left_on = 'pedRLID', right_on = 'fid', how = 'left', indicator=True, suffixes = ('_pave', '_or'))
    assert dfCrossEvents.loc[ dfCrossEvents['_merge']!='both'].shape[0]==0
    dfCrossEvents.drop('_merge', axis=1, inplace=True)

    # Now identify which node of the crossing link ped started at
    dfCrossEvents['start_pave_node'] = dfCrossEvents.apply(lambda row: row['MNodeFID'] if row['node_path'].index(row['MNodeFID']) < row['node_path'].index(row['PNodeFID']) else row['PNodeFID'], axis=1)

    # Join OR junction node to this to find starting road link node
    dfCrossEvents = pd.merge(dfCrossEvents,  gdfPaveNodes.reindex(columns = ['fid', 'juncNodeID']), left_on = 'start_pave_node', right_on = 'fid', how = 'left', indicator=True)
    assert dfCrossEvents.loc[ dfCrossEvents['_merge']!='both'].shape[0]==0
    dfCrossEvents.drop('_merge', axis=1, inplace=True)

    # Get intersection/nearest point between crossing coord string and road link
    dfCrossEvents['cross_linestring'] = dfCrossEvents['CrossingCoordinatesString'].map(lambda x: linestring_from_crossing_coord_string(x))
    dfCrossEvents['rl_cross_point'] = dfCrossEvents.apply(lambda row: road_link_crossing_point(row['cross_linestring'], row['geometry']), axis=1)

    # Get bin of this position
    if nbins is None:
        dfCrossEvents['bin'] = dfCrossEvents.apply(lambda row: crossing_location_distance_bin(row['rl_cross_point'], row['geometry'], row['fid_or'], row['juncNodeID'], gdfORLinks, bin_dist), axis=1)
    else:
        dfCrossEvents['bin'] = dfCrossEvents.apply(lambda row: crossing_location_quantile_bin(row['rl_cross_point'], row['geometry'], row['fid_or'], row['juncNodeID'], gdfORLinks, nbins), axis=1)

    dfCrossEvents.drop(['cross_linestring', 'geometry'], axis=1, inplace=True)

    return dfCrossEvents

def calculate_crossing_location_entropy(dfCrossEvents, dfPedPaths, gdfPaveLinks, gdfPaveNodes, gdfORLinks, dfRun, nbins = None, bin_dist = 2, output_path = "crossing_location_entropy.csv"):
    '''Calculates an entropy measure of the heterogeneity of crossing locations along road links. Does this by binning locations into nbins for each road link 
    and calculating the probability of a crossing occuring in each bin. These probabilities are used to calculate the crossing entropy.
    '''

    if os.path.exists(output_path)==False:

        dfCrossEvents = get_crossing_locations_and_bins(dfCrossEvents, dfPedPaths, gdfPaveLinks, gdfPaveNodes, gdfORLinks, nbins, bin_dist)

        run_bin_counts = dfCrossEvents.groupby(['run'])['bin'].value_counts()
        run_bin_counts.name = 'bin_count'
        dfRunBinCounts = run_bin_counts.reset_index()

        dfRunBinCounts['total'] = dfRunBinCounts.groupby('run')['bin_count'].transform(lambda s: s.sum())
        dfRunBinCounts['pi'] = dfRunBinCounts['bin_count'] / dfRunBinCounts['total']
        dfRunBinCounts['pi_log_pi'] = dfRunBinCounts['pi']* np.log(dfRunBinCounts['pi'])

        # Also calculate entropy across all crossing locations, not aggregated by run
        sim_bin_counts = dfCrossEvents['bin'].value_counts()
        sim_bin_counts.name = 'sim_bin_count'
        dfSimBinCounts = sim_bin_counts.reset_index()
        dfSimBinCounts['sim_total'] = dfSimBinCounts['sim_bin_count'].sum()
        dfSimBinCounts['sim_pi'] = dfSimBinCounts['sim_bin_count'] / dfSimBinCounts['sim_total']
        dfSimBinCounts['sim_pi_log_pi'] = dfSimBinCounts['sim_pi']*np.log(dfSimBinCounts['sim_pi'])
        sim_cross_entropy = dfSimBinCounts['sim_pi_log_pi'].sum() * -1

        dfCrossRunEntropy = dfRunBinCounts.groupby('run')['pi_log_pi'].apply(lambda s: -sum(s)).reset_index().rename(columns = {'pi_log_pi':'cross_entropy'})
        dfCrossRunEntropy['all_run_cross_entropy'] = sim_cross_entropy

        # Finally merge with run parameters
        dfCrossRunEntropy = pd.merge(dfCrossRunEntropy, dfRun, on='run')

        dfCrossRunEntropy.to_csv(output_path, index=False)
    else:
        dfCrossRunEntropy = pd.read_csv(output_path)

    return dfCrossRunEntropy

def all_strategic_path_simple_paths(dfPedRoutes, gdfORLinks, gdfPaveNodes, pavement_graph, simple_paths_path = "simple_paths.csv", weight='length'):
    '''Function returns a dataframe with the IDs of pedestrian agents, their corresponding strategic paths, and all the simple pavement netowkr paths corresponding to each strategic path + start and end pavement node
    '''

    if os.path.exists(simple_paths_path):
        # Load existing data
        print("\nLoading Simple Paths")
        dfPedRoutes = pd.read_csv(output_path)
        dfPedRoutes_removedpeds = pd.read_csv(routes_removed_path)

        # Convert columns to tuples
        dfPedRoutes['FullStrategicPathString'] = dfPedRoutes['FullStrategicPathString'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))
        dfPedRoutes['simple_path'] = dfPedRoutes['edge_path'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))

    else:
        print("\nCalculating Simple Paths")

        # Get unique set of Strategis paths and start and end nodes
        dfUniqueSP = dfPedRoutes.reindex(columns = ['FullStrategicPathString', 'start_node', 'end_node']).drop_duplicates()

        dfSPSimplePaths = pd.DataFrame()
        for ir, row in dfUniqueSP.iterrows():
            sp_simple_paths = simple_paths_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'])
            n_paths = len(sp_simple_paths)
            df  = pd.DataFrame({"FullStrategicPathString": [row['FullStrategicPathString']]*n_paths,
                                "start_node":[row['start_node']]*n_paths,
                                "end_node":[row['end_node']]*n_paths,
                                "simple_path":sp_simple_paths})
            dfSPSimplePaths = pd.concat([dfSPSimplePaths, df])

        # Calculate path lengths
        dfSPSimplePaths['path_length'] = dfSPSimplePaths['simple_path'].map(lambda p: nx.path_weight(pavement_graph, p, weight))

        # Now merge back into original data to get ID
        dfSPSimplePaths = pd.merge(dfPedRoutes, dfSPSimplePaths, on = ['FullStrategicPathString', 'start_node', 'end_node'], how = 'outer', indicator=True)
        #assert dfSPSimplePaths.loc[ dfSPSimplePaths['_merge']!='both'].shape[0]==0
        print(dfSPSimplePaths.loc[ dfSPSimplePaths['_merge']!='both'].shape[0])

        # Save the data
        dfSPSimplePaths.to_csv(simple_paths_path)

    return dfSPSimplePaths

def agg_trip_distance_and_duration(agent_ids_to_exclude, dfRun, routes_path, output_path):
    '''Loads the raw routes data, cleans the data and aggregates trip distances and durations
    '''

    if os.path.exists(output_path)==False:
        # Load the data
        dfRoutes = pd.read_csv(routes_path)

        # Keep only the columns of interest
        dfRoutes = dfRoutes.reindex(columns = ["run", "ID", "JourneyDistance", "JourneyDuration"])
        dfRoutes['JourneySpeed'] = dfRoutes['JourneyDistance'] / dfRoutes['JourneyDuration']

        # Remove ids, typically these are agents that get stuck along their journey
        if agent_ids_to_exclude is not None:
            drop_indices = dfRoutes.loc[ dfRoutes['ID'].isin(agent_ids_to_exclude)].index
            dfRoutes.drop(drop_indices, inplace=True)

        assert dfRoutes.duplicated().any()==False

        # Aggregate trip distance and duration to the run level
        dfDursDists = dfRoutes.groupby("run").agg(  distance=pd.NamedAgg(column="JourneyDistance", aggfunc=lambda s: s.dropna().sum()),
                                                    duration=pd.NamedAgg(column="JourneyDuration", aggfunc=lambda s: s.dropna().sum()),
                                                    speed=pd.NamedAgg(column="JourneySpeed", aggfunc=lambda s: s.dropna().mean()),
                                                    nagents=pd.NamedAgg(column="ID", aggfunc=lambda s: s.unique().shape[0])
                                                ).reset_index()

        # Calculate per agent values
        dfDursDists['DistPA'] = dfDursDists['distance'] / dfDursDists['nagents']
        dfDursDists['DurPA'] = dfDursDists['duration'] / dfDursDists['nagents']

        # Merge in parameter values
        dfDursDists = pd.merge(dfRun, dfDursDists, on = 'run')

        # Save date for future use
        dfDursDists.to_csv(output_path, index=False)
    else:
        dfDursDists = pd.read_csv(output_path)

    return dfDursDists

def figure_rl_paths_heatmap(fig, ax, gdfORLink, gdfStartNodes, gdfEndNodes, graph, dict_node_pos, edgelist, edgedata, edge_cmap, title, cbar_title, title_font, labelsize, fig_config, vlims = None, cbar_pad=0.05, label_pad = 20):
    '''Function for creating figures illustrating tactical path finding
    '''
    if vlims is None:
        vmin = 0
        vmax = 1
    else:
        vmin, vmax = vlims
    xmin, ymin, xmax, ymax = gdfORLink.total_bounds

    #gdfORLink.plot(ax=ax, edgecolor = fig_config['road_link']['color'], linewidth=fig_config['road_link']['linewidth'], linestyle = '-')
    nx.draw_networkx_edges(graph, dict_node_pos, ax = ax, edgelist=edgelist, width = 3, edge_color = edgedata, edge_cmap=edge_cmap, alpha=0.8, edge_vmin = vmin, edge_vmax=vmax)

    gdfStartNodes.plot(ax=ax, edgecolor = fig_config['pavement_node']['path_color'], facecolor = fig_config['pavement_node']['path_color'], linewidth=fig_config['pavement_node']['linewidth'], zorder=8)
    gdfEndNodes.plot(ax=ax, edgecolor = 'red', facecolor = 'red', linewidth=fig_config['pavement_node']['linewidth'], zorder=8)

    #ax.set_title(title, fontdict = title_font, y = 1.03)

    ax.set_xlim(xmin-3, xmax+3)
    ax.set_ylim(ymin-7.5, ymax+7.5)

    ax.set_axis_off()

    # Add colour bar
    smap = plt.cm.ScalarMappable(cmap=edge_cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(smap, ax=ax, fraction=0.1, shrink = 0.3, location='bottom', pad = cbar_pad)
    cbar.ax.tick_params(labelsize=labelsize)
    cbar.ax.set_xticks([vmin, vmax])
    cbar.ax.set_xticklabels(["{:e}".format(round_sig(vmin, sig=3)).replace("0000",""), "{:e}".format(round_sig(vmax, sig=3)).replace("0000","")])
    cbar.ax.set_xlabel(cbar_title) #, rotation=0, labelpad = label_pad)
    return ax

def sobol_si_bar_subplot(ax, dfsi, fig_title, xticklabels):
    bar_width=0.4
    x_pos = np.arange(dfsi.shape[0])
    ax.bar(x_pos, dfsi['S1'], width=bar_width, yerr = dfsi['S1_conf'], align='center', label="S1")
    ax.bar(x_pos+bar_width, dfsi['ST'], width=bar_width, yerr = dfsi['ST_conf'], align='center', label="ST")

    ax.set_xticks(x_pos + bar_width / 2)
    ax.set_xticklabels(xticklabels, rotation=45)
    ax.legend()

    ax.set_title(fig_title)
    return ax

def output_metrics_descriptive_statistics(dfDD, policy_col, metrics):
    '''Group the output metrics by the policy col and perform aggregations to give:
    mean, std, normalised std, standard error in the mean.
    '''

    standard_error = lambda s: s.std() / np.sqrt(s.shape[0])
    standardised_std = lambda s: s.std() / s.mean()

    dfDesc = pd.DataFrame()
    for metric in metrics:
        dfDescM = dfDD.groupby(policy_col).agg( mean=pd.NamedAgg(column=metric, aggfunc=np.mean),
                                                err=pd.NamedAgg(column=metric, aggfunc=standard_error),
                                                std=pd.NamedAgg(column=metric, aggfunc=np.std),
                                                std_st=pd.NamedAgg(column=metric, aggfunc=standardised_std),)
        iterables = [[metric], dfDescM.columns]
        dfDescM.columns = pd.MultiIndex.from_product(iterables, names=["m", "desc"])
        dfDesc = pd.concat([dfDesc, dfDescM], axis=1)

    return dfDesc.loc[['always','sometimes','never']]