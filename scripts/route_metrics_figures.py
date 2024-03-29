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

def hist_plot(ax, data, val_col, group_col, title, nhistbins = 25, palette=['#3d993d', '#cca328', '#af3f33']):

    minv = data[val_col].min()
    maxv = data[val_col].max()
    rng = maxv-minv
    bins = np.linspace(minv, maxv, nhistbins)

    groups = list(data[group_col].unique())
    groups.sort(key=sf)
    for i, g in enumerate(groups):
        df = data.loc[ data[group_col]==g]

        ax.hist(df[val_col].dropna(), bins = bins, density=False, color=palette[i], alpha=0.7, label=g, histtype='step', linewidth=6)

    ax.set_ylabel('count', fontsize = 15)
    #ax.set_xlim(xmin=minv*0.9, xmax=maxv*1.1)
    ax.set_title(title, fontsize = 28, pad=25)
    return ax

def multi_env_hist_plot(f, axs, dfDD, gdfORLinks, outcome_vars, env_col, title_rename_dict, fig_config, nhistbins = 25, palette = ['#3d993d', '#cca328', '#af3f33']):

    data = dfDD.loc[:, outcome_vars+[env_col]]
    for i in range(nvars):
        ax0 = hist_plot(axs[i], data, outcome_vars[i], env_col, title_rename_dict[outcome_vars[i]], nhistbins = nhistbins, palette = palette)

    # Add single legend
    #axs[-1].legend(fontsize = 20, bbox_to_anchor=(1,1), loc="upper left")

    return f

def get_metric_sis_envs(dfDD, problem, env_col, metric):
    dfMetricSIs = pd.DataFrame()
    envs = dfDD[env_col].unique()
    for env in envs:
            X = dfDD.loc[dfDD[env_col]==env, problem['names']].values
            Y = dfDD.loc[dfDD[env_col]==env, metric].values.astype(float)
            if pd.Series(Y).isnull().any():
                print("Null values in utput for {} - {}, skipping".format(cat, metric))
                continue

            Sis = sobol.analyze(problem, Y, calc_second_order=False, num_resamples=100, conf_level=0.95, print_to_console=False, parallel=False, n_processors=None, keep_resamples=False, seed=random_seed)
            df = pd.DataFrame(Sis)
            df['param'] = problem['names']
            df['metric']=metric
            df[env_col]=env
            dfMetricSIs = pd.concat([dfMetricSIs, df])
    return dfMetricSIs

def get_multiple_metrics_sis_envs(dfDD, problem, env_col, outcome_vars):
    # Group by policy setting and calculate sobol indices
    dfSIs = pd.DataFrame()
    for metric in outcome_vars:
        df = get_metric_sis_envs(dfDD, problem, env_col, metric)
        dfSIs = pd.concat([dfSIs, df])
    return dfSIs

def multi_env_sobol_si_plot(f, axs, dfSIs, gdfORLinks, env_col, env_values, outcome_vars, rename_dict, constrained_layout = True, colors = ['#3d993d', '#cca328', '#af3f33']):

    ylims = [(-12, 40), (-12, 40), (-5, 5), (-5, 5)]

    grouped = dfSIs.groupby(['metric'])
    keys = list(grouped.groups.keys())
    keys.sort(key=sf)
    for i, var in enumerate(outcome_vars):
        dfsi = dfSIs.loc[ dfSIs['metric']==var]

        nenvs = len(env_values)
        bar_width=0.8/nenvs
        xi = np.linspace((-bar_width/2) * (nenvs-1), (bar_width/2) * (nenvs-1), nenvs)

        # now loop through environments
        title = "{}".format(rename_dict[var])
        for j, env in enumerate(env_values):
            dfsip = dfsi.loc[ dfsi[env_col]==env]
            data = dfsip.set_index('param')[['ST','ST_conf']].sort_index()
        
            x_pos = np.arange(data.shape[0])+xi[j]
            axs[i].bar(x_pos, data['ST'], width=bar_width, yerr = data['ST_conf'], align='center', label=env, color = colors[j])

        axs[i].set_xticks(np.arange(data.shape[0]))
        axs[i].set_xticklabels([ rename_dict[i] for i in data.index], rotation=45, fontsize=20)
        #axs[i].legend()

        #axs[i].set_title(title, fontsize=24)
    axs[-1].legend(fontsize = 20, bbox_to_anchor=(-2,-0.32,2.5,0), loc="lower center", mode='expand', ncol=3)
    #bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2

    return f

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
nbins = None
bin_dist = 2
env_col = 'environment'

palette = ['#1b9e77', '#d95f02', '#7570b3']

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])
itn_links_file = os.path.join(gis_data_dir, config['mastermap_itn_processed_direction_file'])
itn_nodes_file = os.path.join(gis_data_dir, config['mastermap_node_processed_file'])
crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])

data_paths = bd_utils.get_data_paths(config['cc_results'], data_dir)
batch_file = data_paths["batch_file"]


ug_output_paths = bd_utils.get_ouput_paths(config['ug_results'], "", data_dir, nbins = bin_dist)
ug_output_path = ug_output_paths["output_route_data"]

qg_output_paths = bd_utils.get_ouput_paths(config['qg_results'], "", data_dir, nbins = bin_dist)
qg_output_path = qg_output_paths["output_route_data"]

cc_output_paths = bd_utils.get_ouput_paths(config['cc_results'], "", data_dir, nbins = bin_dist)
cc_output_path = cc_output_paths["output_route_data"]


output_paths = { 'Uniform Grid':    ug_output_path,
                'Quad Grid':        qg_output_path,
                'Clapham Common':   cc_output_path}

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
gdfORLinks = gpd.read_file(or_links_file)


'''
weight_params = range(0, 100, 100)
dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(None, None, None, None, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)
dfCrossEvents = bd_utils.load_and_clean_cross_events(gdfPaveLinks, cross_events_path = cross_events_file)

dfPedRoutesConsistentPeds = dfPedRoutes.loc[ ~dfPedRoutes['ID'].isin(dfPedRoutes_removedpeds['ID'])]
dfCrossEventsConsistentPeds = dfCrossEvents.loc[ ~dfCrossEvents['ID'].isin(dfPedRoutes_removedpeds['ID'])]

dfRouteLength = bd_utils.get_run_total_route_length(dfPedRoutesConsistentPeds, dfRun, pavement_graph, output_path = output_route_length_file)
dfSPSim = bd_utils.get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, range(0,100,100), distance_function = 'dice_dist', output_path = output_sp_similarity_path)
dfSPSimLen = bd_utils.get_shortest_path_similarity(dfPedRoutesConsistentPeds, dfRun, pavement_graph, dict_node_pos, range(0,100,100), distance_function = 'path_length', output_path = output_sp_similarity_length_path)
dfCrossCounts = dfCrossEventsConsistentPeds.merge(dfRun.reindex(columns = ['run','nPeds']), on='run').groupby('run').apply(lambda df: df.shape[0] / df['nPeds'].values[0]).reset_index().rename(columns = {0:'crossCountPP'})
'''

#
# Load data from multiple environments
#
dfDD = pd.DataFrame()
for env, data_path in output_paths.items():
    df = pd.read_csv(data_path)
    df[env_col] = env
    dfDD = pd.concat([dfDD, df])

dfDD  = dfDD.loc[ dfDD['informalCrossing']=='always']

##############################
#
#
# Produce figures
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


plot_types = ['histograph','si']
#outcome_vars = ['DistPA','crossCountPP','cross_entropy']
outcome_vars = ['route_length_pp', 'sp_sim', 'PostponeCountPP']

title_rename_dict = {   "route_length_pp":r"$\bar{L_r}$",
                        "DistPA": r"$\bar{D_r}$",
                        "crossCountPP":r"$\bar{C_r}$",
                        "cross_entropy":r"$CLE$", 
                        'informalCrossing':'Informal Crossing',
                        'PostponeCountPP': 'Postpone Crossing Per Ped',
                        'sp_sim': 'Fractional length difference from shortest path'}

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
                "timeThreshold": r"$\mathrm{\tau}$",
                "route_length_pp":r"$\bar{L_r}$",
                "DistPA": r"$\bar{D_r}$",
                "crossCountPP":r"$\bar{C_r}$",
                "cross_entropy":r"$CLE$", 
                'informalCrossing':'Informal Crossing'
                }

#
# Initialise figure
#
nplots = len(plot_types)
nvars = len(outcome_vars)
fig_width = 8
fig, axs = plt.subplots(nplots, nvars, figsize=(fig_width*len(outcome_vars),13))

#
# Histogram plots
#
#plt.style.use('dark_background')
fig = multi_env_hist_plot(fig, axs[0, :], dfDD, gdfORLinks, outcome_vars, env_col, title_rename_dict, fig_config, palette=palette)

#
# Sobol indices for each metric and policy setting
#
from SALib.analyze import sobol

policy_param = list(policies.keys())[0]
scenario_param_cols =  [i for i in params if i!=policy_param]

env_values = ['Uniform Grid','Quad Grid','Clapham Common']

problem = init_problem(params)

dfSIs = get_multiple_metrics_sis_envs(dfDD, problem, env_col, outcome_vars)
fig = multi_env_sobol_si_plot(fig, axs[1, :], dfSIs, gdfORLinks, env_col, env_values, outcome_vars, rename_dict, constrained_layout = False, colors = palette)

# annotate figure
texts = ['a)','b)','c)','d)','e)','f)']
for i, ax in enumerate(axs.reshape(1,-1)[0]):
    ax.text(-0.15, 0.98, texts[i], transform=ax.transAxes, fontsize=18)

fig.savefig(os.path.join(img_dir, "env_comparison_{}_{}_{}.png".format(config['ug_results'], config['qg_results'], config['cc_results'])))