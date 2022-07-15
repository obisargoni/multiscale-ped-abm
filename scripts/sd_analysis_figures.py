# SD analysis figures

import json
import os
import itertools
import pandas as pd
import numpy as np
import geopandas as gpd
import networkx as nx
import re
from datetime import datetime as dt
from sklearn.neighbors import KernelDensity

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
# Functions
#
#
#####################################
def sf(k):
    if k[1]=='always':
        return 0
    elif k[1]=='sometimes':
        return 1
    else:
        return 2

def pair_plot(dfDD, outcome_vars, policy_col, rename_dict, nbins, file_datetime_string):
    data = dfDD.loc[:, outcome_vars].rename(columns = rename_dict)
    data[rename_dict[policy_col]] = dfDD[policy_col]
    sns.pairplot(data, hue=rename_dict[policy_col], vars=[rename_dict[outcome_vars[0]], rename_dict[outcome_vars[1]]])
    plt.savefig(os.path.join(img_dir, 'pair_plot.{}-{}.{}bins.{}.png'.format(outcome_vars[0], outcome_vars[1], nbins, file_datetime_string)))

def kde_plot(ax, data, val_col, group_col, bandwidth=0.75, palette=['#66ff66', '#ffcc33', '#ff9966']):

    minv = data[val_col].min()
    maxv = data[val_col].max()
    rng = maxv-minv
    x = np.linspace(minv, maxv, 100*int(rng)).reshape(-1,1)

    groups = list(data[group_col].unique())
    groups.sort(key=sf)
    for i, g in enumerate(groups):
        df = data.loc[ data[group_col]==g]

        values = df[val_col].dropna().values.reshape(-1,1)
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(values)
        log_dens = kde.score_samples(x)

        ax.plot(x, np.exp(log_dens), color=palette[i], label=g)
    ax.legend()
    #ax.set_xlim(xmin=minv*0.9, xmax=maxv*1.1)
    ax.set_title(val_col)
    return ax

def hist_plot(ax, data, val_col, group_col, title, nbins=50, palette=['#66ff66', '#ffcc33', '#ff9966']):

    minv = data[val_col].min()
    maxv = data[val_col].max()
    rng = maxv-minv
    bins = np.linspace(minv, maxv, nbins)

    groups = list(data[group_col].unique())
    groups.sort(key=sf)
    for i, g in enumerate(groups):
        df = data.loc[ data[group_col]==g]

        ax.hist(df[val_col].dropna(), bins = bins, density=False, color=palette[i], alpha=1, label=g, histtype='step')

    ax.legend()
    #ax.set_xlim(xmin=minv*0.9, xmax=maxv*1.1)
    ax.set_title(title)
    return ax

def paths_heatmap_data(gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes):
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
    n_paths = dfPedRoutes['run'].unique().shape[0] * dfPedRoutes['ID'].unique().shape[0]
    edge_traverse_counts = dfPathORLinks['FullStrategicPathString'].value_counts()

    path_links = dfPathORLinks['FullStrategicPathString'].unique()
    edgedata = np.array([edge_traverse_counts[i] / n_paths for i in path_links])

    #tp_links = edge_traverse_counts.index()

    edgelist = []
    for edge_id in path_links:
        for e in list(or_graph.edges(data=True)):
            if edge_id == e[-1]['fid']:
                edgelist.append(e)

    return or_graph, dict_node_pos, edgelist, edgedata

def paths_heatmap_figure(fig, ax, or_graph, dict_node_pos, edgelist, edgedata, fig_config, cmap = 'Reds', title = 'Paths Heatmap', cbar_title = "Frequency of link traversal divided by total number of trips", title_font = {'fontsize':15}, labelsize=15, cbar_pad=0.03, label_pad = 20):
    # Get start and end nodes
    gdfStartNodes = gdfPaveNodes.loc[ gdfPaveNodes['fid'].isin(dfPedRoutes['start_node'])]
    gdfEndNodes = gdfPedODs.loc[gdfPedODs['fid'] =='od_0']

    ax = bd_utils.figure_rl_paths_heatmap(fig, ax, gdfORLinks, gdfStartNodes, gdfEndNodes, or_graph, dict_node_pos, edgelist, edgedata, plt.get_cmap(cmap), title, cbar_title, title_font, labelsize, fig_config, cbar_pad = cbar_pad)

    return ax

def combined_kde_path_heatmap_plot(dfDD, outcome_vars, policy_col, nbins, title_rename_dict, fig_config, or_graph, dict_node_pos, edgelist, edgedata, cmap = 'Reds', heatmap_title = 'Paths Heatmap', cbar_title = "Frequency of link traversal divided by total number of trips", title_font = {'fontsize':15}, labelsize=15, cbar_pad=0.03, label_pad = 20, figsize=(20,20)):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1, 2], width_ratios = [1,1])

    #  Varying density along a streamline
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])

    data = dfDD.loc[:, outcome_vars+[policy_col]]
    ax0 = hist_plot(ax0, data, outcome_vars[0], policy_col, title_rename_dict[outcome_vars[0]])
    ax1 = hist_plot(ax1, data, outcome_vars[1], policy_col, title_rename_dict[outcome_vars[1]])
    ax2 = paths_heatmap_figure(fig, ax2, or_graph, dict_node_pos, edgelist, edgedata, fig_config, cmap = cmap, title = heatmap_title, cbar_title = cbar_title, title_font = title_font, labelsize=labelsize, cbar_pad=cbar_pad, label_pad = label_pad)

    outpath = os.path.join(img_dir, 'kde_pathmap.{}-{}.{}bins.{}.png'.format(outcome_vars[0], outcome_vars[1], nbins, file_datetime_string))
    fig.savefig(outpath)

    return outpath

def sobol_si_figure(dfDD, problem, outcome_vars):


    # Group by policy setting and calculate sobol indices
    dfSIs = pd.DataFrame()
    for metric in outcome_vars:
        for pv in policy_values:
            X = dfDD.loc[dfDD[policy_param]==pv, problem['names']].values
            Y = dfDD.loc[dfDD[policy_param]==pv, metric].values.astype(float)
            if pd.Series(Y).isnull().any():
                print("Null values in utput for {} - {}, skipping".format(cat, metric))
                continue

            Sis = sobol.analyze(problem, Y, calc_second_order=False, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)
            df = pd.DataFrame(Sis)
            df['param'] = problem['names']
            df['metric']=metric
            df[policy_param]=pv
            dfSIs = pd.concat([dfSIs, df])

    policy_names = {"never":"Informal crossing prevented", "always": "Informal crossing allowed", "sometimes":"Informal crossing prevented on A-Roads"}

    f, axs = plt.subplots(3,2, figsize=(20,20), constrained_layout=True)
    axs = axs.reshape(3,2)
    ylims = [(-12, 40), (-12, 40), (-5, 5), (-5, 5)]

    grouped = dfSIs.groupby(['metric',policy_param])
    keys = list(grouped.groups.keys())
    keys.sort(key=sf)
    for i, (m, pv) in enumerate(keys):
        p = i//2
        r = i%2
        dfsi = grouped.get_group((m, pv))#.sort_values(by='S1', ascending=False)
        title = "{} - {}".format(m, policy_names[pv])
        ax = bd_utils.sobol_si_bar_subplot(axs[p,r], dfsi, title, dfsi['param'])
        #ax.set_ylim(ylims[i])

    outpath = os.path.join(img_dir,"sobol_si.{}-{}.{}.png".format(outcome_vars[0], outcome_vars[1], file_datetime_string))
    f.savefig(outpath)
    return outpath

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
gdfPaveNodes = gpd.read_file(pavement_nodes_file)
gdfORLinks = gpd.read_file(or_links_file)
gdfORNodes = gpd.read_file(or_nodes_file)
gdfPedODs = gpd.read_file(ped_ods_file)

weight_params = range(0, 100, 100)
dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(None, None, None, None, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)

# Model output data
dfDD = pd.read_csv(output_sd_data)

########################################
#
#
# Processing data to compare runs under different policies
#
# This involves matching up runs with the same parameter settings but under different policy conditions
#
#######################################
print("\nRunning Scenario Discovery Analysis")

# get policy parameter and split the data into groups for different policies
policy_param = list(policies.keys())[0]
policy_values = policies[policy_param]
scenario_param_cols =  [i for i in params if i!=policy_param]

# Now group by scenario and aggregate to find difference in outputs between policy conditions
'''
for c in scenario_param_cols:
    dfDD[c] = dfDD[c].astype(str) # Helps with grouping, makes matching doubles easier

dfPolicyDiff = dfDD.groupby(scenario_param_cols).agg(   PedDistDiff = pd.NamedAgg(column = "DistPAPed", aggfunc=lambda s: s.values[0] - s.values[1]),
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
#assert (dfPolicyDiff['CountRuns']==2).all()
'''

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

#assert 'latin' in config['setting'] # expect LH desig to be used when doing SD

with open("figure_config.json") as f:
    fig_config = json.load(f)


#
# Initial exploratory analysis of multiple outcomes
#
outcome_vars1 = ['route_length_pp','cross_entropy']
outcome_vars2 = ['DistPA','cross_entropy']
policy_col = 'informalCrossing'

title_rename_dict = {   "route_length_pp":"Average Pedestrian Route Length",
                        "DistPA": "Average Pedestrian Distance",
                        "cross_entropy":"Crossing Location Entropy", 
                        'informalCrossing':'Informal Crossing'}
#
# Create pairs plot
#
pair_plot(dfDD, outcome_vars1, policy_col, title_rename_dict, nbins, file_datetime_string)
pair_plot(dfDD, outcome_vars2, policy_col, title_rename_dict, nbins, file_datetime_string)

#
# KDE + paths heatmap plot
#
#plt.style.use('dark_background')
# get data for heatmap of paths
or_graph, dict_node_pos, edgelist, edgedata = paths_heatmap_data(gdfORLinks, gdfORNodes, dfPedRoutes, gdfPaveNodes)
combined_kde_path_heatmap_plot(dfDD, outcome_vars1, policy_col, nbins, title_rename_dict, fig_config, or_graph, dict_node_pos, edgelist, edgedata)
combined_kde_path_heatmap_plot(dfDD, outcome_vars2, policy_col, nbins, title_rename_dict, fig_config, or_graph, dict_node_pos, edgelist, edgedata)

plt.style.use('default')


#
# Sobol indices for each metric and policy setting
#
from SALib.analyze import sobol

policy_param = list(policies.keys())[0]
policy_values = policies[policy_param]
scenario_param_cols =  [i for i in params if i!=policy_param]

problem = init_problem(params)

sobol_si_figure(dfDD, problem, outcome_vars1)
sobol_si_figure(dfDD, problem, outcome_vars2)