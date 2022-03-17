from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import networkx as nx
import itertools
import re
from shapely.geometry import LineString
from shapely.geometry import Point

######################################
#
# Functions
#
######################################

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

    file_re = get_file_regex("pedestrian_routes", file_datetime = file_datetime)
    ped_routes_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

    file_re = get_file_regex("vehicle_routes", file_datetime = file_datetime)
    veh_routes_file = os.path.join(data_dir, most_recent_directory_file(data_dir, file_re))

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
    output["cross_events"]=cross_events_file
    output["pedestrian_locations"]=ped_locations_file
    output["vehicle_road_links"]=vehicle_rls_file
    output["batch_file"]=batch_file

    return output

def get_ouput_paths(file_datetime_string, vehicle_density_timestamp,  data_dir):
    # output paths for processed data
    paths = {}
    paths["output_ped_routes_file"] = os.path.join(data_dir, "ped_routes.{}.csv".format(file_datetime_string))
    paths["output_vehicle_density_file"] = os.path.join(data_dir, "av_vehicle_density.{}.csv".format(vehicle_density_timestamp))
    paths["output_route_length_file"] = os.path.join(data_dir, "run_route_length.{}.csv".format(file_datetime_string))
    paths["output_sp_similarity_path"] = os.path.join(data_dir, "sp_similarity.{}.csv".format(file_datetime_string))
    paths["output_sp_similarity_length_path"] = os.path.join(data_dir, "path_length_sp_similarity.{}.csv".format(file_datetime_string))
    paths["output_route_completion_path"] = os.path.join(data_dir, "route_completions.{}.csv".format(file_datetime_string))
    paths["output_cross_events_path"] = os.path.join(data_dir, "cross_events.{}.csv".format(file_datetime_string))
    paths["output_ks_factormap"] = os.path.join(data_dir , "ks_factor_map.{}.csv".format (file_datetime_string))
    paths["output_corr_factormap"] = os.path.join(data_dir , "corr_factor_map.{}.csv".format (file_datetime_string))
    paths["output_ped_distdurs_file"] = os.path.join(data_dir, "ped_durdists.{}.csv".format(file_datetime_string))
    paths["output_veh_distdurs_file"] = os.path.join(data_dir, "veh_durdists.{}.csv".format(file_datetime_string))

    paths["output_sd_data"] = os.path.join(data_dir, "metrics_for_sd_analysis.{}.csv".format(file_datetime_string))

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
        dfRouteLength = pd.merge(dfRun, dfRouteLengthPerPed, on = 'run')

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

def load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = "ped_routes.csv", strategic_path_filter = True):

    d,f = os.path.split(ped_routes_path)
    output_path = os.path.join(d, 'processed_'+f)

    split_path = os.path.splitext(ped_routes_path)
    routes_removed_path = split_path[0] + "_removed_peds" + split_path[1]

    if os.path.exists(output_path):
        print("\nLoading Pedestrian Agent Routes")
        dfPedRoutes = pd.read_csv(output_path)
        dfPedRoutes_removedpeds = pd.read_csv(routes_removed_path)

        # Convert columns to tuples
        dfPedRoutes['FullStrategicPathString'] = dfPedRoutes['FullStrategicPathString'].map(lambda s: tuple(s.strip("('").strip("')").split("', '")))
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
        dfPedRoutes['edge_path'] = dfPedRoutes['FullTacticalRouteString'].map(lambda s: tuple(dict.fromkeys(s.strip(":").split(":"))))

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
                                                        'start_node', 'end_node'])

        for k in weight_params:
            weight_name = "weight{}".format(k)
            gdfPaveLinks['cross_cost'] = gdfPaveLinks['AvVehDen'].fillna(0) * k
            gdfPaveLinks[weight_name] = gdfPaveLinks['length'] + gdfPaveLinks['cross_cost']

            # Weight pavement network crossing links by average vehicle flow
            weight_attributes = gdfPaveLinks.set_index( gdfPaveLinks.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
            nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

            if strategic_path_filter:
                dfUniqueStartEnd['sp_{}'.format(k)] = dfUniqueStartEnd.apply(lambda row: shortest_path_within_strategic_path(row['FullStrategicPathString'], gdfORLinks, gdfPaveNodes, pavement_graph, row['start_node'], row['end_node'], weight = weight_name), axis=1)
            else:
                dfUniqueStartEnd['sp_{}'.format(k)] = dfUniqueStartEnd.apply(lambda row: nx.dijkstra_path(pavement_graph, row['start_node'], row['end_node'], weight =  weight_name), axis=1)

        dfPedRoutes = pd.merge(dfPedRoutes, dfUniqueStartEnd, on = ['start_node', 'end_node', 'FullStrategicPathString'])

        dfPedRoutes.to_csv(output_path, index=False)
        dfPedRoutes_removedpeds.to_csv(routes_removed_path, index=False)

    return dfPedRoutes, dfPedRoutes_removedpeds

def load_and_clean_cross_events(gdfPaveLinks, cross_events_path = "cross_events.csv"):
    '''Method to aggregate crossing events from the Ped Crossings dataset, to produce dataframe with single row per crossing event.
    Crossing event defined by crossing type, location and minimum TTC during crossing.
    '''
    d,f = os.path.split(cross_events_path)
    output_path = os.path.join(d, 'processed_'+f)

    if os.path.exists(output_path)==False:

        print("\nProcessing Pedestrian Cross Events")

        dfCrossEvents = pd.read_csv(cross_events_path)

        dfCrossEvents = dfCrossEvents.reindex(columns = ['run', 'ID', 'FullStrategicPathString', 'CrossingType', 'TacticalEdgeID', 'CrossingCoordinatesString', 'TTC'])

        # Group by run, ID and TacticalEdgeID to find min TTC per cross event
        dfCrossEvents['TTC'] = dfCrossEvents.groupby(['run', 'ID', 'CrossingType', 'TacticalEdgeID'])['TTC'].transform(lambda s: s.min())

        # Drop duplicates again now that TTC and crossing coord processed
        dfCrossEvents = dfCrossEvents.drop_duplicates()

        # Merge with pave links to get link type
        dfLinkTypes = gdfPaveLinks.reindex(columns = ['fid','linkType']).drop_duplicates()
        dfCrossEvents = pd.merge(dfCrossEvents, dfLinkTypes, left_on = 'TacticalEdgeID', right_on = 'fid', how = 'left')
        dfCrossEvents.drop('fid', axis=1,inplace=True)

        # Check that there is a single crossing event per ped per pavement link.
        cross_per_ped_link = dfCrossEvents.groupby(['run', 'ID', 'TacticalEdgeID']).apply(lambda df: df.shape[0])
        assert cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]==0

        dfCrossEvents.to_csv(output_path, index=False)
    else:
        print("\nLoading Existing Processed Pedestrian Cross Events")
        dfCrossEvents = pd.read_csv(output_path)

    return dfCrossEvents

def agg_trip_distance_and_duration(agent_ids_to_exclude, dfRun, routes_path, output_path):
    '''Loads the raw routes data, cleans the data and aggregates trip distances and durations
    '''

    if os.path.exists(output_path)==False:
        # Load the data
        dfRoutes = pd.read_csv(routes_path)

        # Keep only the columns of interest
        dfRoutes = dfRoutes.reindex(columns = ["run", "ID", "JourneyDistance", "JourneyDuration"])

        # Remove ids, typically these are agents that get stuck along their journey
        if agent_ids_to_exclude is not None:
            drop_indices = dfRoutes.loc[ dfRoutes['ID'].isin(agent_ids_to_exclude)].index
            dfRoutes.drop(drop_indices, inplace=True)

        assert dfRoutes.duplicated().any()==False

        # Aggregate trip distance and duration to the run level
        dfDursDists = dfRoutes.groupby("run").agg(  distance=pd.NamedAgg(column="JourneyDistance", aggfunc=lambda s: s.dropna().sum()),
                                                    duration=pd.NamedAgg(column="JourneyDuration", aggfunc=lambda s: s.dropna().sum()),
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