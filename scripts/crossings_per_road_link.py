# Script to analyse which road links pedestrian agents cross on

import json
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
import networkx as nx
from datetime import datetime as dt
from shapely.geometry import Point

import batch_data_utils as bd_utils

#####################################
#
# Globals
#
#####################################
project_crs = {'init': 'epsg:27700'}

with open("config.json") as f:
    config = json.load(f)

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])


# Model output data
file_datetime_string = "2021.Sep.07.12_38_34"
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
dfPedCrossingsRaw = pd.read_csv(ped_crossings_file)
dfRun = pd.read_csv(os.path.join(data_dir, batch_file))

# GIS Data
gdfPaveNetwork = gpd.read_file(pavement_links_file)
gdfORLinks = gpd.read_file(or_links_file)
gdfORNodes = gpd.read_file(or_nodes_file)

# Get networkx graph and node positions
G = nx.Graph()
edges = gdfORLinks.loc[:, ['MNodeFID', 'PNodeFID', 'fid']].values
G.add_weighted_edges_from(edges, weight='fid')

# Using the geographical coordinates of the nodes when plotting them
points_pos = gdfORNodes.set_index('node_fid')
points_pos['x'] = points_pos['geometry'].map(lambda g: g.coords[0][0])
points_pos['y'] = points_pos['geometry'].map(lambda g: g.coords[0][1])
node_posistions = list(zip(points_pos['x'], points_pos['y']))
dict_node_pos = dict(zip(points_pos.index, node_posistions))

#######################################
#
#
# Process data
#
# Calculating number of crossings per road link. Normalising by number of peds on road link.
#
#
########################################

# Process raw model output to 
# - drop rows that don't correspond to a crossing event
# - aggregate TTC data to choose the lowest TTC value per crossing event
# - this produces dataset with 1 row per crossing event
# - Drop duplicates, expect just one crossing per ped per pavement link - can't include crossing coords string bc sometimes it includes both crossing coords
# - Check that there is one crossing per link per ped

# check there won't be any ttC data lost when excluding 'none' crossing events
assert dfPedCrossingsRaw.loc[ (dfPedCrossingsRaw['ChosenCrossingTypeString']=='none') & (~dfPedCrossingsRaw['TTC'].isnull()) ].shape[0]==0
assert dfPedCrossingsRaw['CurrentPavementLinkID'].isnull().any()==False

dfCrossEvents = dfPedCrossingsRaw.loc[dfPedCrossingsRaw['ChosenCrossingTypeString']!='none'].reindex(columns = ['run', 'ID', 'FullStrategicPathString', 'ChosenCrossingTypeString', 'CurrentPavementLinkID', 'CrossingCoordsString', 'TTC'])
dfCrossEvents = dfCrossEvents.drop_duplicates()

# Group by run, ID and CurrentPavementLinkID to find crossing coordinates and lowest TTC
dfCrossEvents['TTC'] = dfCrossEvents.groupby(['run', 'ID', 'ChosenCrossingTypeString', 'CurrentPavementLinkID'])['TTC'].transform(lambda s: s.min())
dfCrossEvents['CrossingCoordsString'] = dfCrossEvents.groupby(['run', 'ID', 'ChosenCrossingTypeString', 'CurrentPavementLinkID'])['CrossingCoordsString'].transform(lambda s: max(s.dropna(), key=len) if ~s.isnull().all() else None)

# Drop duplicates again now that TTC and crossing coord processed
dfCrossEvents = dfCrossEvents.drop_duplicates()

# Check that there is a single crossing event per ped per pavement link.
cross_per_ped_link = dfCrossEvents.groupby(['run', 'ID', 'CurrentPavementLinkID']).apply(lambda df: df.shape[0])
assert cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]==0


#
# Now start aggregating to OR road link
#

# combine with pavement network data to aggregate crossings per OR link
dfORCrossEvents = dfCrossEvents.merge(gdfPaveNetwork, left_on = 'CurrentPavementLinkID', right_on = 'fid', how = 'left', indicator=True)
assert dfORCrossEvents.loc[ dfORCrossEvents['_merge']!='both'].shape[0]==0
dfORCrossEvents.drop('_merge', axis=1, inplace=True)
dfORCrossEvents['FullStrategicPathString'] = dfORCrossEvents['FullStrategicPathString'].map(lambda s: s.strip(':').split(':'))

# Aggregate crossing counts
dfCrossingCounts = dfORCrossEvents.groupby(['run', 'pedRLID']).apply(lambda g: g.shape[0]).reset_index()
dfCrossingCounts.rename(columns = {0:'cross_count'}, inplace=True)

# Aggregate unmarked crossing counts
dfUmCC = dfORCrossEvents.loc[dfORCrossEvents['ChosenCrossingTypeString']=='unmarked'].groupby(['run', 'pedRLID']).apply(lambda g: g.shape[0]).reset_index()
dfUmCC.rename(columns = {0:'um_cross_count'}, inplace=True)
dfCrossingCounts = pd.merge(dfCrossingCounts, dfUmCC, on = ['run', 'pedRLID'], how = 'left')

# Get number of pedestrians per road link, use this to normalise crossing counts - crossings per ped on link
# Given fixed ODs and flows should get same number of peds per road link but different numbers crossing.
road_links = pd.Series(np.concatenate(dfORCrossEvents['FullStrategicPathString'].values))
peds_per_link = road_links.value_counts()
peds_per_link.name = 'peds_per_link'

# Join with pavement Road link data and run data
gdfCrossingCounts = pd.merge(gdfORLinks, dfCrossingCounts, left_on = 'fid', right_on = 'pedRLID', how = 'left')
gdfCrossingCounts['cross_count'] = gdfCrossingCounts['cross_count'].fillna(0)
gdfCrossingCounts['um_cross_count'] = gdfCrossingCounts['um_cross_count'].fillna(0)
gdfCrossingCounts = pd.merge(gdfCrossingCounts, dfRun, on = 'run')

# Join with number of peds per link and calculate crossings per ped
gdfCrossingCounts = pd.merge(gdfCrossingCounts, peds_per_link, left_on = 'fid', right_index = True)
gdfCrossingCounts['cross_count_pp'] = gdfCrossingCounts['cross_count'] / gdfCrossingCounts['peds_per_link']
gdfCrossingCounts['um_cross_count_pp'] = gdfCrossingCounts['um_cross_count'] / gdfCrossingCounts['peds_per_link']


# Need to average over random seeds before creating colour map
cols = list(dfRun.columns) + ['pedRLID_x', 'peds_per_link', 'cross_count_pp', 'um_cross_count_pp']
gdfCrossingCounts = gdfCrossingCounts.reindex( columns = cols)

group_cols = ['epsilon', 'lambda', 'addPedTicks', 'alpha', 'addVehicleTicks', 'tacticalPlanHorizon', 'gamma', 'minCrossingProp', 'pedRLID_x']
gdfCrossingCountsAv = gdfCrossingCounts.groupby(group_cols)['cross_count_pp'].apply(lambda s: s.mean()).reset_index()
gdfCrossingCountsAv['um_cross_count_pp'] = gdfCrossingCounts.groupby(group_cols)['um_cross_count_pp'].transform(lambda s: s.mean()).fillna(0)
gdfCrossingCountsAv = pd.merge(gdfORLinks.reindex(columns = ['fid','geometry']), gdfCrossingCountsAv, left_on = 'fid', right_on = 'pedRLID_x', how = 'left')

# Map crossing counts to range for colormap
max_cc_pp = gdfCrossingCountsAv['cross_count_pp'].max()
gdfCrossingCountsAv['cmap_value'] = gdfCrossingCountsAv['cross_count_pp'].map(lambda c: 255*(c/max_cc_pp))

um_max_cc_pp = gdfCrossingCountsAv['um_cross_count_pp'].max()
gdfCrossingCountsAv['um_cmap_value'] = gdfCrossingCountsAv['um_cross_count_pp'].map(lambda c: 255*(c/max_cc_pp))



#
# Aggregate conflicts to pavement links
#
calc_conflict_count = lambda s: s.dropna().shape[0]
calc_mean_ttc = lambda s: s.dropna().mean()
calc_var_ttc = lambda s: s.dropna().var()

dfConflictCounts = dfORCrossEvents.groupby(['run', 'CurrentPavementLinkID']).agg(   conflict_count=pd.NamedAgg(column="TTC", aggfunc=calc_conflict_count),
                                                                                    meanTTC=pd.NamedAgg(column="TTC", aggfunc=calc_mean_ttc),
                                                                                    varTTC=pd.NamedAgg(column="TTC", aggfunc=calc_var_ttc),
                                                                                    )

# Now average over random seeds
dfConflictCounts = pd.merge(dfRun, dfConflictCounts, on = 'run')
param_cols = [c for c in dfRun.columns if c not in ['run', 'randomSeed']]


dfAvConflictCounts = dfConflictCounts.groupby(['run', 'CurrentPavementLinkID']).agg(    mean_conflict_count=pd.NamedAgg(column="conflict_count", aggfunc=np.mean),
                                                                                        totMeanTTC=pd.NamedAgg(column='meanTTC', aggfunc=np.mean),
                                                                                        totVarTTC=pd.NamedAgg(column="varTTC", aggfunc=np.mean),
                                                                                        )

# Merge with pavement links to get the geometries
gdfAvConflictCounts = pd.merge(dfAvConflictCounts, gdfPaveNetwork.reindex(columns = ['fid','geometry']), left_on = 'CurrentPavementLinkID', right_on = 'fid')
gdfAvConflictCounts['linkMidPoint'] = gdfAvConflictCounts['geometry'].map(lambda g: g.centroid)

#
# Also process data to get coordinates of crossing points so these can be mapped
#
coord_regex = re.compile(r"(\d{6}.\d+)")
dfCrossEvents['first_coord'] = dfCrossEvents['CrossingCoordsString'].map(lambda s: coord_regex.findall(s.split("),(")[0]))
dfCrossEvents['geometry'] = dfCrossEvents['first_coord'].map(lambda p:  Point(*map(float, p)))

dfPedCrossingsRun = pd.merge(dfCrossEvents, dfRun, on = 'run', how = 'left')

gdfPedCrossingRun = gpd.GeoDataFrame(dfPedCrossingsRun, geometry = 'geometry')
for run in gdfPedCrossingRun['run'].unique():
    gdfSub = gdfPedCrossingRun.loc[ gdfPedCrossingRun['run']==run, ['run','ChosenCrossingTypeString', 'geometry']]
    gdfSub.to_file("..\\output\\crossing_points_run{}.shp".format(run))

##########################################
#
#
# Plot network
#
#
##########################################

# Plot
from matplotlib import cm # for generating colour maps
from matplotlib import pyplot as plt

def road_network_figure(G, dict_node_pos, dict_edge_values, title, cmap_name = 'viridis', edge_width = 3, edge_alpha = 1):

    plt.style.use('dark_background')
    f, ax = plt.subplots(1,1,figsize = (15,15))

    ax = road_network_subfigure(ax, G, dict_node_pos, dict_edge_values, title, cmap_name = cmap_name, edge_width=edge_width, edge_alpha=edge_alpha)
    return f

def road_network_subfigure(ax, G, dict_node_pos, dict_edge_values, title, cmap_name = 'viridis', edge_width = 3, edge_alpha = 1, title_font = {'size': 12}):
    # Get edge colour map based on number of crossings
    cmap = cm.get_cmap(cmap_name)
    edge_palette = []
    for e in G.edges(data=True):
        fid = e[-1]['fid']
        if fid in dict_edge_values.keys():
            value = dict_edge_values[fid]
        else:
            value = 0

        edge_palette.append(cmap(value))

    nx.draw_networkx_nodes(G, dict_node_pos, ax = ax, nodelist=G.nodes(), node_color = 'grey', node_size = 1, alpha = 0.5)
    nx.draw_networkx_edges(G, dict_node_pos, ax = ax, edgelist=G.edges(), width = 3, edge_color = edge_palette, alpha=1)
    ax.set_title(title, fontdict = title_font)
    ax.axis('off')
    return ax

def batch_runs_road_network_figure(G, dict_node_pos, dfBatch, groupby_columns, fig_title, value_col = 'cmap_value', cmap_name = 'viridis', edge_width = 3, edge_alpha = 1, sub_title_font = {'size': 12}, title_font = {'size': 16}):
    '''Loop through batch run groups and get edge pallet data for each group. Use this to make road crossings
    figure for each group.
    '''

    grouped = dfBatch.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of
    p = len(dfBatch[groupby_columns[0]].unique())
    q = len(dfBatch[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    plt.style.use('dark_background')
    f, axs = plt.subplots(p, q, figsize = (20,20), sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            dfRun = grouped.get_group(group_key)
            dict_edge_values = dfRun.set_index('fid')[value_col].to_dict()

            rename_dict = {'0':'minimise distance', '1':'minimise crossings'}
            title = "Tactical planning horizon: {}".format(group_key[0])+r"$\degree$"+"\nTactical planning heuristic: {}".format(rename_dict[str(int(group_key[1]))])

            # Select the corresponding axis
            ax = axs[pi, qi]

            road_network_subfigure(ax, G, dict_node_pos, dict_edge_values, title, cmap_name = cmap_name, edge_width = 3, edge_alpha = 1, title_font = sub_title_font)

    f.suptitle(fig_title, fontdict = title_font)
    return f

def batch_runs_crossing_points_figure(G, dict_node_pos, dfBatch, groupby_columns, fig_title, point_size = 3, point_col_dict = {'unmarked':'red','unsignalised':'green'}, point_alpha = 0.5, edge_width = 3, edge_alpha = 1, sub_title_font = {'size': 12}, title_font = {'size': 16}):
    '''Loop through batch run groups and get edge pallet data for each group. Use this to make road crossings
    figure for each group.
    '''

    grouped = dfBatch.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of
    p = len(dfBatch[groupby_columns[0]].unique())
    q = len(dfBatch[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    plt.style.use('dark_background')
    f, axs = plt.subplots(p, q, figsize = (20,20), sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            dfRun = grouped.get_group(group_key)

            rename_dict = {'0':'minimise distance', '1':'minimise crossings'}
            title = "Tactical planning horizon: {}".format(group_key[0])+r"$\degree$"+"\nTactical planning heuristic: {}".format(rename_dict[str(int(group_key[1]))])

            # Select the corresponding axis
            ax = axs[pi, qi]

            # Draw the road network
            nx.draw_networkx_nodes(G, dict_node_pos, ax = ax, nodelist=G.nodes(), node_color = 'grey', node_size = 1, alpha = 0.5)
            nx.draw_networkx_edges(G, dict_node_pos, ax = ax, edgelist=G.edges(), width = 3, edge_color = 'grey', alpha=1)
            ax.set_title(title, fontdict = title_font)
            ax.axis('off')

            # Draw the crossing points
            x = dfRun.loc[ dfRun['ChosenCrossingTypeString']=='unmarked', 'geometry'].map(lambda p:p.x)
            y = dfRun.loc[ dfRun['ChosenCrossingTypeString']=='unmarked', 'geometry'].map(lambda p:p.y)
            ax.scatter(x, y, s = point_size, marker = 'o', alpha = point_alpha, c = point_col_dict['unmarked'])

            x = dfRun.loc[ dfRun['ChosenCrossingTypeString']=='unsignalised', 'geometry'].map(lambda p:p.x)
            y = dfRun.loc[ dfRun['ChosenCrossingTypeString']=='unsignalised', 'geometry'].map(lambda p:p.y)
            ax.scatter(x, y, s = point_size, marker = 'o', alpha = point_alpha, c = point_col_dict['unsignalised'])

    f.suptitle(fig_title, fontdict = title_font)
    return f

def batch_runs_conflicts_figure(G, dict_node_pos, dfBatch, groupby_columns, fig_title, geometry_col = 'geometry', size_col = 'size', size_param = 5, colour_col = 'colour', cmap = 'viridis', norm = None, point_alpha = 0.5, edge_width = 3, edge_alpha = 1, sub_title_font = {'size': 12}, title_font = {'size': 16}):
    '''Loop through batch run groups and get edge pallet data for each group. Use this to make road crossings
    figure for each group.
    '''

    grouped = dfBatch.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of
    p = len(dfBatch[groupby_columns[0]].unique())
    q = len(dfBatch[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    plt.style.use('dark_background')
    f, axs = plt.subplots(p, q, figsize = (20,20), sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            dfRun = grouped.get_group(group_key)

            rename_dict = {'0':'minimise distance', '1':'minimise crossings'}
            title = "Tactical planning horizon: {}".format(group_key[0])+r"$\degree$"+"\nTactical planning heuristic: {}".format(rename_dict[str(int(group_key[1]))])

            # Select the corresponding axis
            ax = axs[pi, qi]

            # Draw the road network
            nx.draw_networkx_nodes(G, dict_node_pos, ax = ax, nodelist=G.nodes(), node_color = 'midnightblue', node_size = 1, alpha = 0.5)
            nx.draw_networkx_edges(G, dict_node_pos, ax = ax, edgelist=G.edges(), width = 3, edge_color = 'midnightblue', alpha=1)
            ax.set_title(title, fontdict = sub_title_font)
            ax.axis('off')

            # Draw the crossing points
            x = dfRun[geometry_col].map(lambda p:p.x)
            y = dfRun[geometry_col].map(lambda p:p.y)
            sizes = dfRun[size_col]*size_param
            colours = dfRun[size_col]
            ax.scatter(x, y, s = sizes, cmap = cmap, norm = norm, marker = 'o', alpha = point_alpha, c = colours, edgecolors='none')

    # Add colourbar
    cbar_fontdict = {"size":14}
    smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = f.colorbar(smap, ax=axs, fraction=0.1, shrink = 0.8)
    cbar.ax.tick_params(labelsize=cbar_fontdict['size']-3)
    cbar.ax.set_ylabel('Mean TTC', rotation=-90, labelpad = 15, fontdict = cbar_fontdict)

    # Select some example point to create size legend from later
    size_data = dfBatch[size_col].values
    example_values = np.sort(size_data)[::len(size_data)//10][-3:]
    indices = [np.where(size_data==v)[0][0] for v in example_values]

    # Then create legend using an inset axis
    xmin, ymin, dx, dy = cbar.ax.get_position().bounds
    cax = f.add_axes([xmin+0.05, ymin+dy-dy/5, dx, dy/5])
    x = [0]*len(example_values)
    y = range(len(example_values))
    cax.scatter(x, y, s = example_values*size_param, c = 'white', edgecolors = 'none', marker = 'o')
    cax.yaxis.set_label_position("right")
    cax.yaxis.tick_right()
    cax.set_yticks(y)
    cax.set_yticklabels(np.round_(example_values, decimals=1), fontdict = {'size':cbar_fontdict['size']-3})
    cax.set_ylabel('Mean Conflict Count', rotation=-90, labelpad = 15, fontdict = cbar_fontdict)
    cax.set_xticks([])

    for pos in ['right', 'top', 'bottom', 'left']:
        cax.spines[pos].set_visible(False)


    f.suptitle(fig_title, fontsize = title_font['size'])
    return f

groupby_columns = ['tacticalPlanHorizon', 'minCrossingProp']
batch_fig = batch_runs_road_network_figure(G, dict_node_pos, gdfCrossingCounts, groupby_columns, "Crossings per pedestrian on road link", value_col = 'cmap_value', cmap_name = 'viridis', edge_width = 3, edge_alpha = 1)

batch_fig_path = os.path.join(img_dir, "network_figure_"+ os.path.split(ped_crossings_file)[1].replace(".csv", ".png"))
batch_fig.savefig(batch_fig_path)
batch_fig.show()


um_batch_fig = batch_runs_road_network_figure(G, dict_node_pos, gdfCrossingCounts, groupby_columns, "Informal crossings per pedestrian on road link", value_col = 'um_cmap_value', cmap_name = 'plasma', edge_width = 3, edge_alpha = 1)

um_batch_fig_path = os.path.join(img_dir, "network_figure_informal_"+ os.path.split(ped_crossings_file)[1].replace(".csv", ".png"))
um_batch_fig.savefig(um_batch_fig_path)
um_batch_fig.show()

points_batch_fig_path = os.path.join(img_dir, "crossing_points_figure"+ os.path.split(ped_crossings_file)[1].replace(".csv", ".png"))
points_batch_fig = batch_runs_crossing_points_figure(G, dict_node_pos, dfPedCrossingsRun, groupby_columns, "Crossing Locations", point_size = 3, point_col_dict = {'unmarked':'red', 'unsignalised':'green'}, point_alpha = 0.5, edge_width = 3, edge_alpha = 1, sub_title_font = {'size': 12}, title_font = {'size': 16})
points_batch_fig.savefig(points_batch_fig_path)
points_batch_fig.show()


gdfConflicts = gdfAvConflictCounts.loc[ gdfAvConflictCounts['alpha'] == 0.8]

conflicts_batch_fig_path = os.path.join(img_dir, "conflicts_figure"+ os.path.split(ped_crossings_file)[1].replace(".csv", ".png"))
conflicts_batch_fig = batch_runs_conflicts_figure(G, dict_node_pos, gdfConflicts, groupby_columns, "Conflict Locations", geometry_col = 'linkMidPoint', size_col = 'mean_conflict_count', colour_col = 'totMeanTTC', cmap = 'OrRd', norm = None, point_alpha = 0.7, edge_width = 3, edge_alpha = 1, sub_title_font = {'size': 12}, title_font = {'size': 20})
conflicts_batch_fig.savefig(conflicts_batch_fig_path)
conflicts_batch_fig.show()