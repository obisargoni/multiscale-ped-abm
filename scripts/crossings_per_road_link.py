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

with open(".//gis_data_processing//config.json") as f:
    config = json.load(f)

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])


# Model output data
file_datetime_string = "2021.Jun.04.12_48_52"
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
#
########################################

# Process raw model output to 
# - just include pavement links used for crossing and with a crossing coordinate and a crossing type
# - Drop duplicates, expect just one crossing per ped per pavement link - can't include crossing coords string bc sometimes it includes both crossing coords
# - Check that there is one crossing per link per ped
dfPedCrossings = dfPedCrossings.loc[(~dfPedCrossings['CurrentPavementLinkID'].isnull())]

dfPedCrossings = dfPedCrossings.drop_duplicates(subset = ['run', 'ID', 'ChosenCrossingTypeString', 'CurrentPavementLinkID'])

# Need to account for cases where two crossing types recorded for a pavement link.

cross_per_ped_link = dfPedCrossings.groupby(['run', 'ID', 'CurrentPavementLinkID']).apply(lambda df: df.shape[0])
#assert cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]==0

# Now combine with pavement network data to aggregate crossings per OR link

# Merge with pavement network to find the OR road link being crossed
dfPedCrossingsLinks = dfPedCrossings.merge(gdfPaveNetwork, left_on = 'CurrentPavementLinkID', right_on = 'fid', how = 'left', indicator=True)
assert dfPedCrossingsLinks.loc[ dfPedCrossingsLinks['_merge']!='both'].shape[0]==0
dfPedCrossingsLinks.drop('_merge', axis=1, inplace=True)
dfPedCrossingsLinks['FullStrategicPathString'] = dfPedCrossingsLinks['FullStrategicPathString'].map(lambda s: s.strip('-').split('-'))

# Aggregate crossing counts
dfCrossingCounts = dfPedCrossingsLinks.groupby(['run', 'pedRLID']).apply(lambda g: g.shape[0]).reset_index()
dfCrossingCounts.rename(columns = {0:'cross_count'}, inplace=True)

# Aggregate unmarked crossing counts
dfUmCC = dfPedCrossingsLinks.loc[dfPedCrossingsLinks['ChosenCrossingTypeString']=='unmarked'].groupby(['run', 'pedRLID']).apply(lambda g: g.shape[0]).reset_index()
dfUmCC.rename(columns = {0:'um_cross_count'}, inplace=True)
dfCrossingCounts = pd.merge(dfCrossingCounts, dfUmCC, on = ['run', 'pedRLID'], how = 'left')

# Get number of pedestrians per road link, use this to normalise crossing counts - crossings per ped on link
# Given fixed ODs and flows should get same number of peds per road link but different numbers crossing.
road_links = pd.Series(np.concatenate(dfPedCrossingsLinks['FullStrategicPathString'].values))
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

# Map crossing counts to range for colormap
max_cc_pp = gdfCrossingCounts['cross_count_pp'].max()
gdfCrossingCounts['cmap_value'] = gdfCrossingCounts['cross_count_pp'].map(lambda c: 255*(c/max_cc_pp))

um_max_cc_pp = gdfCrossingCounts['um_cross_count_pp'].max()
gdfCrossingCounts['um_cmap_value'] = gdfCrossingCounts['um_cross_count_pp'].map(lambda c: 255*(c/max_cc_pp))


# Also process data to get coordinates of crossing points so these can be mapped

# Filter out records without primary crossings
dfPedCrossings = dfPedCrossings.loc[    (~dfPedCrossings['CrossingCoordsString'].isnull())  &
                                        (dfPedCrossings['ChosenCrossingTypeString']!='none')]

cross_per_ped_link = dfPedCrossings.groupby(['run', 'ID', 'CurrentPavementLinkID']).apply(lambda df: df.shape[0])
assert cross_per_ped_link.loc[ cross_per_ped_link!=1].shape[0]==0

coord_regex = re.compile(r"(\d{6}.\d+)")
dfPedCrossings['first_coord'] = dfPedCrossings['CrossingCoordsString'].map(lambda s: coord_regex.findall(s.split("),(")[0]))
dfPedCrossings['geometry'] = dfPedCrossings['first_coord'].map(lambda p:  Point(*map(float, p)))

dfPedCrossingsRun = pd.merge(dfPedCrossings, dfRun, on = 'run', how = 'left')

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