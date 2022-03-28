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
import matplotlib.transforms as transforms

import batch_data_utils as bd_utils

#################################
#
#
# Functions
#
#
#################################
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

def sobol_second_order_si_bar_figure(dfsi, fig_title, rename_dict):

    xticklabels = dfsi['j'].replace(rename_dict) + " - " + dfsi['i'].replace(rename_dict)

    f, ax = plt.subplots(1,1, figsize = (10,10))
    #bar_width=0.4
    x_pos = np.arange(dfsi.shape[0])
    ax.bar(x_pos, dfsi['S2'], width=bar_width, yerr = dfsi['S2_conf'], align='center', label="S1")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(xticklabels, rotation=45)

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

def plot_layers(ax, config, pavement = None, carriageway = None, road_link = None, road_node = None, rays = None, pavement_link = None, pavement_node = None):
    '''Keyword aarguments are geodataframes containing the shapes to be plotted
    '''
    for i, (k, v) in enumerate(locals().items()):

        # Skip no keywork arguments
        if k in ['ax', 'config']:
            continue
        if v is not None:

            if k in ['pavement','carriageway']:
                v.plot(ax=ax, color = config[k]['color'], zorder=i)
            elif k in ['road_link', 'pavement_link']:
                v.plot(ax=ax, facecolor=config[k]['color'], edgecolor = config[k]['color'], linewidth=config[k]['linewidth'], zorder=i)
            elif k in ['road_node', 'pavement_node']:
                v.plot(ax=ax, facecolor=config[k]['color'], edgecolor = config[k]['color'], linewidth=config[k]['linewidth'], zorder=i)
            elif k in ['rays']:
                v.plot(ax=ax, color=config[k]['color'], linewidth=config[k]['linewidth'], zorder=i)
            else:
                v.plot(ax=ax)

    return ax

def figure_single_ped_tactical_paths(study_area_rls, origin_id, dest_id, sp, gdfTopoVeh, gdfTopoPed, gdfORNode, gdfORLink, gdfPedODs, pavement_graph, dict_node_pos, edgelist, edgedata, edge_cmap, ped_links_exclude, fig_config, rotation=0):
    '''Function for creating figures illustrating tactical path finding
    '''

    # Initialise figure
    f, ax = plt.subplots(1,1, figsize = (10,10))

    # Get study area gdfs
    gdfORLinkSA = gdfORLink.loc[ gdfORLink['fid'].isin(study_area_rls)]
    study_area_nodes = np.concatenate( [gdfORLinkSA['MNodeFID'].values, gdfORLinkSA['PNodeFID'].values] )
    gdfORNodeSA = gdfORNode.loc[ gdfORNode['node_fid'].isin(study_area_nodes)]
    gdfTopoPedSA = gdfTopoPed.loc[gdfTopoPed['roadLinkID'].isin(study_area_rls)]
    gdfTopoVehSA = gdfTopoVeh.loc[gdfTopoVeh['roadLinkID'].isin(study_area_rls)]

    # Select route layers
    gdfsp = gdfORLink.loc[ gdfORLink['fid'].isin( sp )]

    gdfods = gdfPedODs.loc[ gdfPedODs['fid'].isin( [origin_id, dest_id])]

    vmin = min(edgedata)
    vmax = max(edgedata)

    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(rotation)

    # plot these additional layers
    gdfsp.plot(ax=ax, edgecolor = 'black', linewidth=fig_config['road_link']['linewidth'], linestyle = '-', zorder=7, transform = rot+base)

    #nx.draw_networkx_nodes(G, dict_node_pos, ax = ax, nodelist=G.nodes(), node_color = 'grey', node_size = 1, alpha = 0.2)
    nx.draw_networkx_edges(pavement_graph, dict_node_pos, ax = ax, edgelist=edgelist, width = 3, edge_color = edgedata, edge_cmap=edge_cmap, alpha=0.8, edge_vmin = vmin, edge_vmax=vmax, transform = rot+base)

    gdfods.plot(ax=ax, edgecolor = fig_config['od']['color'], facecolor = fig_config['od']['color'], linewidth=fig_config['od']['linewidth'], zorder=9, transform = rot+base)

    # Add colour bar
    smap = plt.cm.ScalarMappable(cmap=edge_cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(smap, ax=axs, fraction=0.1, shrink = 0.8)
    cbar.ax.tick_params(labelsize=25)

    # Set limits
    xmin, ymin, xmax, ymax = gdfsp.total_bounds
    ax.set_xlim(xmin-3, xmax+3)
    ax.set_ylim(ymin-7.5, ymax+7.5)

    ax.set_axis_off()

    return f

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
vehicle_topographic_file = os.path.join(gis_data_dir, config['topo_vehicle_processed_file'])
pedestrian_topographic_file = os.path.join(gis_data_dir, config['topo_pedestrian_processed_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])

# Model output data
data_paths = bd_utils.get_data_paths(file_datetime_string, data_dir)
ped_routes_file = data_paths["pedestrian_routes"]
veh_routes_file = data_paths["vehicle_routes"]
av_vehicle_counts_file = data_paths['av_vehicle_counts']
cross_events_file = data_paths["cross_events"]
batch_file = data_paths["batch_file"] 

# output paths for processed data
output_paths = bd_utils.get_ouput_paths(file_datetime_string, vehicle_density_timestamp, data_dir)
output_ped_routes_file=             output_paths["output_ped_routes_file"]
output_single_ped_links_file=       output_paths["output_single_ped_links_file"]
output_vehicle_density_file=        output_paths["output_vehicle_density_file"]
output_route_length_file=           output_paths["output_route_length_file"]
output_sp_similarity_path=          output_paths["output_sp_similarity_path"]
output_sp_similarity_length_path=   output_paths["output_sp_similarity_length_path"]
output_route_completion_path=       output_paths["output_route_completion_path"]
output_cross_events_path=           output_paths["output_cross_events_path"]
output_ks_factormap=                output_paths["output_ks_factormap"]
output_corr_factormap=              output_paths["output_corr_factormap"]


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

weight_params = range(0, 100, 100)

######################################
#
#
# Process data to get ped routes and shortest path routes
#
#
######################################
dfVehCounts = bd_utils.get_road_link_vehicle_density_from_vehicle_counts(gdfITNLinks, av_vehicle_counts_file, output_vehicle_density_file)
dfVehCounts.rename(columns = {'fid':'itn_fid'}, inplace=True)

dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, range(0,100,100), dfVehCounts = dfVehCounts, ped_routes_path = ped_routes_file)
dfSinglePedPaths, ped_id_simple_paths = bd_utils.median_ped_pavement_link_counts(dfPedRoutes, output_path = output_single_ped_links_file)

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
dfConflicts = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds, dfRun, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsMarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unsignalised'], dfRun, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsUnmarked = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['CrossingType']=='unmarked'], dfRun, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDirect = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='direct_cross'], dfRun, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonal = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ dfCrossEventsConsistentPeds['linkType']=='diag_cross'], dfRun, dfLinkCrossCounts, ttc_col = 'TTC')
dfConflictsDiagonalUm = bd_utils.agg_cross_conflicts(dfCrossEventsConsistentPeds.loc[ (dfCrossEventsConsistentPeds['linkType']=='diag_cross') & (dfCrossEventsConsistentPeds['CrossingType']=='unmarked')], dfRun, dfLinkCrossCounts, ttc_col = 'TTC')

conflicts_data = {'all':dfConflicts, 'marked':dfConflictsMarked, 'unmarked':dfConflictsUnmarked, 'direct':dfConflictsDirect, 'diag':dfConflictsDiagonal, 'diag_um':dfConflictsDiagonalUm}

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
from sobol_plot import plot_sobol_indices, save_second_order_sobol_indices
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

# Create an excel writer to record the senstitivity indices
xlWriter = pd.ExcelWriter(os.path.join(data_dir, "{}_results.{}.xlsx".format(setting, file_datetime_string)), mode='w', engine_kwargs=None)


######################################
#
#
# Factor mapping
#
#
######################################

if 'monte_carlo_filtering' in setting:
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

if "morris_factor_fixing" in setting:

    print("\nCalculating sensitivity indices - Route completions")

    X_rc = dfRouteCompletion.loc[:, problem['names']].values
    Y_rc = dfRouteCompletion.loc[:, 'frac_completed_journeys'].values
    Sis = morris.analyze(problem, X_rc, Y_rc, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)

    # Gather into a dataframe
    dfcompsi = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_compsi = morris_si_bar_figure(dfcompsi, r"Jouney Completion $\mathrm{\mu^*}$", 'Fraction journeys completed', dfcompsi['names'].replace(rename_dict))
    f_compsi.savefig(os.path.join(img_dir, "route_completion_sis.{}.png".format(file_datetime_string)))
    f_compsi.clear()
    dfcompsi.to_excel(xlWriter, sheet_name = 'frac_completed_journeys')


    print("\nCalculating sensitivity indices - Conflicts")
    title_dict = config["title_dict"]
    for cat in config["conflict_categories"]:
        dfC = conflicts_data[cat]
        for metric in config["conflict_metrics"]:
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
            f_ccsi.savefig(os.path.join(img_dir, "{}_{}_sis.{}.png".format(metric, cat, file_datetime_string)))
            f_ccsi.clear()
            df.to_excel(xlWriter, sheet_name = "{}_{}_sis".format(metric, cat))



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
        f_spsi.savefig(os.path.join(img_dir, "sp_similarity_sis_{}.{}.png".format(k, file_datetime_string)))
        f_spsi.clear()
        dfspsi.to_excel(xlWriter, sheet_name = "sp_similarity_sis_{}".format(k))


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
        f_spsi.savefig(os.path.join(img_dir, "sp_len_similarity_sis_{}.{}.png".format(k, file_datetime_string)))
        f_spsi.clear()
        dfspsi.to_excel(xlWriter, sheet_name = "sp_len_similarity_sis_{}".format(k))


    print("\nCalculating sensitivity indices - Total route length")

    # Get array of parameter values and output values
    X = dfRouteLength.loc[:, problem['names']].values
    Y = dfRouteLength.loc[:, 'route_length_pp'].astype(float).values
    try:
        Sis = morris.analyze(problem, X, Y, num_resamples = 100, conf_level= 0.95, print_to_console = False, num_levels = num_levels, seed=random_seed)
    except ValueError as e:
        print(e)
        print(k)

    # Gather into a dataframe
    dfRLSis = pd.DataFrame(Sis).sort_values(by='mu_star', ascending=False)
    f_rlsi = morris_si_bar_figure(dfRLSis, "Mean Path Length Sensitivity", r"$\mathrm{\mu^*}$", dfRLSis['names'].replace(rename_dict))
    f_rlsi.savefig(os.path.join(img_dir, "route_length_pp_sis.{}.png".format(file_datetime_string)))
    f_rlsi.clear()
    dfRLSis.to_excel(xlWriter, sheet_name = "route_length_pp_sis_{}".format(k))


if 'sobol_si' in setting:

    print("\nCalculating sobol indices - Conflicts")
    title_dict = config["title_dict"]
    for cat in config["conflict_categories"]:
        dfC = conflicts_data[cat]
        for metric in config["conflict_metrics"]:
            X = dfC.loc[:, problem['names']].values
            Y = dfC.loc[:, metric].values.astype(float)
            if pd.Series(Y).isnull().any():
                print("Null values in utput for {} - {}, skipping".format(cat, metric))
                continue
                
            Sis = sobol.analyze(problem, Y, calc_second_order=calc_second_order, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)

            # Gather into a dataframe
            Sis['names'] = problem['names']
            Sis_filtered = {k:Sis[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf', 'names']}
            df = pd.DataFrame(Sis_filtered).sort_values(by='S1', ascending=False)
            df.to_excel(xlWriter, sheet_name="{}_{}_sobol1T".format(metric, cat))

            # Create figures
            f_si = sobol_si_bar_figure(df, "{} Sobol Indices - {} crossings".format(title_dict[metric], cat), df['names'].replace(rename_dict))
            f_si.savefig(os.path.join(img_dir, "{}_{}_sobol1T.{}.png".format(metric, cat, file_datetime_string)))
            f_si.clear()

            if calc_second_order==True:
                f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.001, rename_dict = rename_dict)
                f_si.savefig(os.path.join(img_dir, "{}_{}_sobol2T.{}.png".format(metric, cat, file_datetime_string)))
                f_si.clear()
                dfsi_second_order = save_second_order_sobol_indices(xlWriter, "{}_{}_sobol2T".format(metric, cat), Sis, problem)

                f_si = sobol_second_order_si_bar_figure(dfsi_second_order,  "{} Second Order Sobol Indices - {} crossings".format(title_dict[metric], cat), rename_dict)
                f_si.savefig(os.path.join(img_dir, "{}_{}_sobol2T_bar.{}.png".format(metric, cat, file_datetime_string)))


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
        df.to_excel(xlWriter, sheet_name = "sp_similarity_sobol_{}".format(k))

        # Plot
        f_si = sobol_si_bar_figure(df, "Shortest Path Similarity Sobol Indices", df['names'].replace(rename_dict))
        f_si.savefig(os.path.join(img_dir, "sp_similarity_sobol_{}.{}.png".format(k, file_datetime_string)))
        f_si.clear()

        if calc_second_order==True:
            f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.001, rename_dict = rename_dict)
            f_si.savefig(os.path.join(img_dir, "sp_similarity_sobol_2ndorder_{}.{}.png".format(k, file_datetime_string)))
            f_si.clear()
            dfsi_second_order = save_second_order_sobol_indices(xlWriter, "sp_similarity_sobol_2ndorder_{}".format(k), Sis, problem)

            f_si = sobol_second_order_si_bar_figure(dfsi_second_order,  "Shortest Path Similarity Second Order Sobol Indices", rename_dict)
            f_si.savefig(os.path.join(img_dir, "sp_similarity_sobol_2ndorder_bar_{}.{}.png".format(k, file_datetime_string)))

    print("Calculating Sobol indices - Total route length")

    X = dfRouteLength.loc[:, problem['names']].values
    Y = dfRouteLength.loc[:, 'route_length_pp'].astype(float).values
    try:
        Sis = sobol.analyze(problem, Y, calc_second_order=calc_second_order, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)
    except ValueError as e:
        print(e)
        print(k)


    # Gather into a dataframe
    Sis['names'] = problem['names']
    Sis_filtered = {k:Sis[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf', 'names']}
    df = pd.DataFrame(Sis_filtered).sort_values(by='S1', ascending=False)
    df.to_excel(xlWriter, sheet_name = "route_length_pp_sobol")

    # Plot
    f_si = sobol_si_bar_figure(df, "Mean Path Length Sensitivity", df['names'].replace(rename_dict))
    f_si.savefig(os.path.join(img_dir, "route_length_pp_sobol.{}.png".format(file_datetime_string)))
    f_si.clear()

    # If second_order then produce plot show interdependence of parameter sensitivity
    if calc_second_order==True:
        f_si = plot_sobol_indices(Sis, problem, criterion='ST', threshold=0.01, rename_dict = rename_dict)
        f_si.savefig(os.path.join(img_dir, "route_length_pp_sobol_2ndorder.{}.png".format(file_datetime_string)))
        f_si.clear()
        dfsi_second_order = save_second_order_sobol_indices(xlWriter, "route_length_pp_sobol_2ndorder", Sis, problem)

        f_si = sobol_second_order_si_bar_figure(dfsi_second_order,  "Mean Path Length Second Order Sensitivity", rename_dict)
        f_si.savefig(os.path.join(img_dir, "route_length_pp_sobol_2ndorder_bar.{}.png".format(file_datetime_string)))


#########################################
#
#
# Scatter plot of two variables, coloured by output metric
#
#
#########################################

if "epsilon_gamma_scatter" in setting:

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

    dfCrossAtTarget = dfCrossEventsCompleteJourney.loc[dfCrossEventsCompleteJourney['TacticalEdgeID'].isin(target_links)]
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

if 'variance_comparison' in setting:

    print("\nProducing single agents paths figure")

    with open("figure_config.json") as f:
        fig_config = json.load(f)

    gdfTopoVeh = gpd.read_file(vehicle_topographic_file)
    gdfTopoPed = gpd.read_file(pedestrian_topographic_file)
    gdfPedODs = gpd.read_file(ped_ods_file)


    assert dfSinglePedPaths['start_node'].unique().shape[0]==1
    assert dfSinglePedPaths['end_node'].unique().shape[0]==1

    start_node = dfSinglePedPaths['start_node'].unique()[0]
    end_node = dfSinglePedPaths['end_node'].unique()[0]

    study_area_rls = dfSinglePedPaths['FullStrategicPathString'].values[0]
    origin_id = gdfPaveNodes.loc[ gdfPaveNodes['fid']==start_node, 'juncNodeID'].values[0]
    dest_id = gdfPaveNodes.loc[ gdfPaveNodes['fid']==end_node, 'juncNodeID'].values[0]

    sp = study_area_rls
    tp_links = dfSinglePedPaths['edge_path'].unique()

    edgelist = []
    for edge_id in tp_links:
        for e in list(pavement_graph.edges(data=True)):
            if edge_id == e[-1]['fid']:
                edgelist.append(e)

    # Now get link counts to colour figure by
    # Aggregate single ped links to get edge data values
    edge_traverse_counts = dfSinglePedPaths['edge_path'].value_counts()
    edgedata = np.array([edge_traverse_counts[i] for i in tp_links])

    f_single_pad_paths = figure_single_ped_tactical_paths(study_area_rls, origin_id, dest_id, sp, gdfTopoVeh, gdfTopoPed, gdfORNodes, gdfORLinks, gdfPedODs, pavement_graph, dict_node_pos, edgelist, edgedata, plt.get_cmap('Reds'), [], fig_config)
    f_single_pad_paths.tight_layout()
    output_single_pad_paths = os.path.join(img_dir, "single_ped_paths.{}.png".format(file_datetime_string))
    f_single_pad_paths.savefig(output_single_pad_paths)

    
    # Not until here that we want to calculate tactical paths using alternative model
    alt_model_paths = []

    rng=np.random.RandomState(100)
    eps = rng.normal(0, 0.25, size=100)

    # Want to vary crossing cost on direct crossings and diagonal crossings separately.

    # To do that ideantify direct crossings that correspond to crossing infrastructure
    gdfPaveLinksDirectCrossings = gdfPaveLinks.loc[ gdfPaveLinks['linkType'] == 'direct_cross']
    gdfPaveLinksCAs = gpd.sjoin(gdfPaveLinks.loc[ gdfPaveLinks['linkType'] == 'direct_cross'], gdfCAs, op = 'within')
    direct_crossing_with_marked_cas = gdfPaveLinksCAs['fid'].unique()

    # Create alternative set of paths for each vehicle flow setting
    for run in dfRun.drop_duplicates(subset = 'addVehicleTicks')['run'].unique():
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

                    # Initialise cross cost as zero
                    dfLinkWeights['cross_cost'] = 0

                    # Then for road crossing link set crossing cost as a multiple of the average vehicle desnity on the link the crossing is on.
                    dfLinkWeights.loc[ dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'cross_cost'] = dfLinkWeights.loc[ dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'AvVehDen'] * k
                    dfLinkWeights.loc[ ~dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'cross_cost'] = dfLinkWeights.loc[ ~dfLinkWeights['fid'].isin(direct_crossing_with_marked_cas), 'AvVehDen'] * j
                    dfLinkWeights[weight_name] = dfLinkWeights['length'] + dfLinkWeights['cross_cost']

                # sense check the weights by comparing length to cross cost
                print(weight_name)
                print( (dfLinkWeights['cross_cost'] / dfLinkWeights['length']).describe())

                # Weight pavement network crossing links by average vehicle flow
                weight_attributes = dfLinkWeights.set_index( dfLinkWeights.apply(lambda row: (row['MNodeFID'], row['PNodeFID']), axis=1))[weight_name].to_dict()
                nx.set_edge_attributes(pavement_graph, weight_attributes, name = weight_name)

                alt_model_path = bd_utils.shortest_path_within_strategic_path(sp, gdfORLinks, gdfPaveNodes, pavement_graph, start_node, end_node, weight = weight_name)
                alt_model_paths.append(alt_model_path)

    # Create complementary figure for alternative model paths
    edgelist = [list(zip(i[:-1], i[1:])) for i in alt_model_paths]
    
    # Create series of all (u,v) tuple edges
    all_edges = np.concatenate(edgelist)
    s_edgelist = pd.Series(tuple(i) for i in all_edges)

    # Count number of times each edge is traversed
    edge_counts = s_edgelist.value_counts()

    f_single_alt_paths = figure_single_ped_tactical_paths(study_area_rls, origin_id, dest_id, sp, gdfTopoVeh, gdfTopoPed, gdfORNodes, gdfORLinks, gdfPedODs, pavement_graph, dict_node_pos, edge_counts.index, edge_counts.values, plt.get_cmap('Reds'), [], fig_config)
    f_single_alt_paths.tight_layout()
    output_single_pad_paths = os.path.join(img_dir, "single_ped_alt_model_paths.{}.png".format(file_datetime_string))
    f_single_alt_paths.savefig(output_single_pad_paths)

    # Now compare path length varaition
    clt_path_lengths = [nx.path_weight(pavement_graph, p, 'length') for p in dfPedRoutes.loc[ dfPedRoutes['ID']==ped_id_simple_paths, 'node_path']]
    alt_model_lengths = [nx.path_weight(pavement_graph, p, 'length') for p in alt_model_paths]

    s_clt = pd.Series(clt_path_lengths)
    s_alt = pd.Series(alt_model_lengths)
    print(s_clt.describe())
    print(s_alt.describe())


xlWriter.close()