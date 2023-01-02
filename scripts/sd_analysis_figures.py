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

def sf2(k):
    if k=='always':
        return 0
    elif k=='sometimes':
        return 1
    else:
        return 2

def pair_plot(dfDD, outcome_vars, policy_col, rename_dict, file_datetime_string):
    data = dfDD.loc[:, outcome_vars].rename(columns = rename_dict)
    data[rename_dict[policy_col]] = dfDD[policy_col]
    sns.pairplot(data, hue=rename_dict[policy_col], vars=[rename_dict[outcome_vars[0]], rename_dict[outcome_vars[1]]])
    plt.savefig(os.path.join(img_dir, 'pair_plot.{}-{}.{}.png'.format(outcome_vars[0], outcome_vars[1], file_datetime_string)))

def kde_plot(ax, data, val_col, group_col, bandwidth=0.75, palette=['#1b9e77', '#d95f02', '#7570b3']):

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

def hist_plot(ax, data, val_col, group_col, title, nhistbins = 25, palette=['#1b9e77', '#d95f02', '#7570b3']):

    minv = data[val_col].min()
    maxv = data[val_col].max()
    rng = maxv-minv
    bins = np.linspace(minv, maxv, nhistbins)

    groups = list(data[group_col].unique())
    groups.sort(key=sf)
    for i, g in enumerate(groups):
        df = data.loc[ data[group_col]==g]

        ax.hist(df[val_col].dropna(), bins = bins, density=False, color=palette[i], alpha=1, label=g, histtype='step', linewidth=6)

    ax.legend(fontsize = 20)
    ax.set_ylabel('count', fontsize = 20)
    #ax.set_xlim(xmin=minv*0.9, xmax=maxv*1.1)
    ax.set_title(title, fontsize = 24)
    return ax

def multi_hist_plot(dfDD, gdfORLinks, outcome_vars, policy_col, title_rename_dict, fig_config, inset_rec, nhistbins = 25, figsize=(20,10), ttc_threshold=1):
    nvars = len(outcome_vars)
    if nvars==4:
        fig, axs = plt.subplots(2, 2, figsize=figsize)
        axs = axs.reshape(1,-1)[0]
    else:
        fig, axs = plt.subplots(1, nvars, figsize=figsize)

    data = dfDD.loc[:, outcome_vars+[policy_col]]
    for i in range(nvars):
        ax0 = hist_plot(axs[i], data, outcome_vars[i], policy_col, title_rename_dict[outcome_vars[i]], nhistbins = nhistbins)

    # add inset showing the road network
    axins = fig.add_axes(inset_rec)
    gdfORLinks.plot(ax=axins, color='black')
    axins.set_axis_off()
    axins.set_title('Environment', y=-0.1)

    outpath = os.path.join(img_dir, 'hists_ttc{}.{}.png'.format(ttc_threshold,file_datetime_string))
    fig.savefig(outpath)

    return outpath

def get_metric_sis(dfDD, problem, policy_param, policy_values, metric):
    dfMetricSIs = pd.DataFrame()
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
            dfMetricSIs = pd.concat([dfMetricSIs, df])
    return dfMetricSIs

def get_multiple_metrics_sis(dfDD, problem, policy_param, policy_values, outcome_vars):
    # Group by policy setting and calculate sobol indices
    dfSIs = pd.DataFrame()
    for metric in outcome_vars:
        df = get_metric_sis(dfDD, problem, policy_param, policy_values, metric)
        dfSIs = pd.concat([dfSIs, df])
    return dfSIs

def sobol_si_figure(dfSIs, gdfORLinks, policy_param, policy_values, outcome_vars, rename_dict, inset_rec, constrained_layout = True, fig_width = 10, colors = ['#1b9e77', '#d95f02', '#7570b3'], ttc_threshold=1):
    nvars = len(outcome_vars)
    if nvars==4:
        f, axs = plt.subplots(2,2, figsize=(20,20), constrained_layout = constrained_layout)
        axs = axs.reshape(1,-1)[0]
    else:
        f, axs = plt.subplots(1,nvars, figsize=(fig_width*nvars,10), constrained_layout = constrained_layout)
    

    ylims = [(-12, 40), (-12, 40), (-5, 5), (-5, 5)]

    for i, m in enumerate(outcome_vars):
        dfsi = dfSIs.loc[ dfSIs['metric']==m]

        n_policies = len(policy_values)
        bar_width=0.8/n_policies
        xi = np.linspace((-bar_width/2) * (n_policies-1), (bar_width/2) * (n_policies-1), n_policies)

        # now loop through policy settings
        title = "{}".format(rename_dict[m])
        for j, pv in enumerate(policy_values):
            dfsip = dfsi.loc[ dfsi[policy_param]==pv]
            data = dfsip.set_index('param')[['ST','ST_conf']].sort_index()
        
            x_pos = np.arange(data.shape[0])+xi[j]
            axs[i].bar(x_pos, data['ST'], width=bar_width, yerr = data['ST_conf'], align='center', label=pv, color = colors[j])

        axs[i].set_xticks(np.arange(data.shape[0]))
        axs[i].set_xticklabels([ rename_dict[i] for i in data.index], rotation=45, fontsize=20)
        axs[i].legend()

        axs[i].set_title(title, fontsize=24)


    # add inset showing the road network
    axins = f.add_axes(inset_rec)
    gdfORLinks.plot(ax=axins, color='black')
    axins.set_axis_off()
    axins.set_title('Environment', y=-0.1)


    outpath = os.path.join(img_dir,"sobol_si_ttc{}.{}.png".format(ttc_threshold, file_datetime_string))
    f.savefig(outpath)
    return outpath


def interval_index_from_groups(groups):
    intervals = [pd.Interval(left=i, right=j, closed = 'right') for (i,j) in zip(groups[:-1], groups[1:])]
    interval_index = pd.IntervalIndex(intervals, verify_integrity=True)
    return interval_index



def agg_policy_comparison_figure(dfDD, gdfORLinks, group_param, policy_param, metric, rename_dict, inset_rec, title, colors = ['#1b9e77', '#d95f02', '#7570b3'], figsize = (15,7), quantile_groups = (0.25,0.75,1.0), quantile_labels = ("Bottom 25%", "Middle 50%", "Top 25%") ):

    cut_values = dfDD[group_param].drop_duplicates().quantile(quantile_groups).tolist()
    cut_values = [0] + cut_values
    cut_bins = interval_index_from_groups(cut_values)
    dfDD['group_level'] = pd.cut(dfDD[group_param], bins = cut_bins)
    dfDD['group_label'] = dfDD['group_level'].replace(dict(zip(cut_bins, quantile_labels)))

    # Get mean conflict counts for each flow level and 
    dfcomp = dfDD.groupby(['group_label',policy_param]).agg( av = pd.NamedAgg(column = metric, aggfunc=np.mean), err=pd.NamedAgg(column = metric, aggfunc=lambda x: np.std(x) / np.sqrt(x.shape[0]) ) ).reset_index()

    f, ax = plt.subplots(figsize = figsize)

    group_values = dfcomp['group_label'].unique()
    policy_values = dfcomp[policy_param].unique().tolist()
    policy_values.sort(key=sf2)

    x = np.arange(1, len(group_values)+1)


    for i, p in enumerate(policy_values):
        dfp = dfcomp.loc[ dfcomp[policy_param]==p]
        dx = -0.2+ (0.2*i)

        xi = x+dx

        ax.errorbar(xi, dfp['av'], yerr = dfp['err'], c = colors[i], fmt="x", markersize=15, label=p)

    ax.set_xticks(x)
    ax.set_xticklabels(group_values)
    ax.tick_params(axis='both', labelsize=15)

    ax.legend(fontsize=20)

    ax.set_title(title, fontsize=24)

    ax.set_ylabel(rename_dict[metric], fontsize=18)
    ax.set_xlabel(rename_dict[group_param], fontsize=18, labelpad=20)


    # add inset showing the road network
    axins = f.add_axes(inset_rec)
    gdfORLinks.plot(ax=axins, color='black')
    axins.set_axis_off()
    axins.set_title('Environment', y=-0.1)

    outpath = os.path.join(img_dir,"agg_comparison_{}.{}.png".format(metric,file_datetime_string))
    f.savefig(outpath)
    return f

def agg_policy_two_metric_comparison_figure(dfDD, gdfORLinks, group_param, policy_param, metrics, rename_dict, inset_rec, title, colors = ['#1b9e77', '#d95f02', '#7570b3'], figsize = (15,7), quantile_groups = (0.25,0.75,1.0), quantile_labels = ("Bottom 25%", "Middle 50%", "Top 25%"), ttc_threshold=1 ):

    cut_values = dfDD[group_param].drop_duplicates().quantile(quantile_groups).tolist()
    cut_values = [0] + cut_values
    cut_bins = interval_index_from_groups(cut_values)
    dfDD['group_level'] = pd.cut(dfDD[group_param], bins = cut_bins)
    dfDD['group_label'] = dfDD['group_level'].replace(dict(zip(cut_bins, quantile_labels)))

    f, axs = plt.subplots(2,1,figsize = figsize)

    handles = []
    labels = []

    for k, ax in enumerate(axs):
        
        # Get mean conflict counts for each flow level and policy
        dfcomp = dfDD.groupby(['group_label',policy_param]).agg( av = pd.NamedAgg(column = metrics[k], aggfunc=np.mean), err=pd.NamedAgg(column = metrics[k], aggfunc=lambda x: np.std(x) / np.sqrt(x.shape[0]) ) ).reset_index()

        group_values = dfcomp['group_label'].unique()
        policy_values = dfcomp[policy_param].unique().tolist()
        policy_values.sort(key=sf2)

        x = np.arange(1, len(group_values)+1)


        for i, p in enumerate(policy_values):
            dfp = dfcomp.loc[ dfcomp[policy_param]==p]
            dx = -0.25+ (0.25*i)

            xi = x+dx

            ebar = ax.errorbar(xi, dfp['av'], yerr = dfp['err'], c = colors[i], fmt='x', markersize=15, label=p)
            handles.append(ebar)
            labels.append(rename_dict[metrics[k]] + ", " + p)

        ax.set_ylabel(rename_dict[metrics[k]], fontsize=18)
        ax.set_xlabel(rename_dict[group_param], fontsize=18, labelpad=20)
        ax.tick_params(axis='both', labelsize=15)

        ax.set_xticks(x)
        ax.set_xticklabels(group_values)
    
    axs[0].legend(fontsize=20)#, loc='upper left', bbox_to_anchor=(-0.2, 0.8))
    axs[0].set_xticklabels([])
    axs[0].set_xlabel(None)

    f.suptitle(title, fontsize=24)

    # add inset showing the road network
    axins = f.add_axes(inset_rec)
    gdfORLinks.plot(ax=axins, color='black')
    axins.set_axis_off()
    axins.set_title('Environment', y=-0.15)

    outpath = os.path.join(img_dir,"agg_comparison_ttc{}_{}_{}.{}.png".format(ttc_threshold, metrics[0],metrics[1],file_datetime_string))
    f.savefig(outpath)
    return f

#####################################
#
#
# Choose data to analyze and sensitivity analysis setting
#
#
#####################################
file_datetime_string = config['file_datetime_string']
setting = config['setting']

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
nbins = None
bin_dist = 2
ttc_threshold = 1

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

output_paths = bd_utils.get_ouput_paths(file_datetime_string, data_dir, nbins = bin_dist, ttc_threshold=ttc_threshold)
output_sd_data = output_paths["output_sd_data"]

palette = ['#1b9e77', '#d95f02', '#7570b3']

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
#dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(None, None, None, None, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)

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
outcome_vars1 = ['route_length_pp','mean_link_cross_entropy']
outcome_vars2 = ['DistPA','mean_link_cross_entropy']
outcome_vars3 = ['route_length_pp', 'speedVeh','conflict_count','mean_link_cross_entropy']
policy_col = 'informalCrossing'

title_rename_dict = {   "route_length_pp":r"$\bar{L_r}$",
                        "DistPA": r"$\bar{D_r}$",
                        "crossCountPP":r"$\bar{C_r}$",
                        "conflict_count":r"$\bar{C_r}$",
                        "cross_entropy":r"$CLE$", 
                        "mean_link_cross_entropy":r"$CLE$", 
                        "speedVeh":r"$\bar{S^v_r}$",
                        'informalCrossing':'Informal Crossing'}
#
# Create pairs plot
#
#pair_plot(dfDD, outcome_vars1, policy_col, title_rename_dict, file_datetime_string)
#pair_plot(dfDD, outcome_vars2, policy_col, title_rename_dict, file_datetime_string)

#
# Histogram plots
#
#plt.style.use('dark_background')
inset_rec = [-0.02, 0.87, 0.13, 0.13]
multi_hist_plot(dfDD, gdfORLinks, outcome_vars3, policy_col, title_rename_dict, fig_config, inset_rec, figsize=(20,20), ttc_threshold=ttc_threshold)

plt.style.use('default')


#
# Sobol indices for each metric and policy setting
#
from SALib.analyze import sobol

rename_dict = { 'alpha':r"$\mathrm{\alpha}$",
                'lambda':r"$\mathrm{\lambda}$",
                "epsilon":r"$\mathrm{\epsilon}$",
                "gamma":r"$\mathrm{\gamma}$",
                "minCrossing": r"$\mathrm{MC}$",
                "tacticalPlanHorizon": r"$\mathrm{PH}$",
                "avNVehicles": r"$\bar{N^v}_r$",
                "addPedTicks": r"$\mathrm{T_{ped}}$",
                "pedSpeedSeed": r"$\mathrm{Seed_{pSpeed}}$",
                "pedMassSeed": r"$\mathrm{Seed_{pMass}}$",
                "caSampleSeed": r"$\mathrm{Seed_{CA}}$",
                "vehODSeed": r"$\mathrm{Seed_{veh}}$",
                "timeThreshold": r"$\mathrm{\tau}$",
                "route_length_pp":r"$\bar{L_r}$",
                "DistPA": r"$\bar{D_r}$",
                "crossCountPP":r"$\bar{C_r}$",
                "conflict_count":r"$\bar{C_r}$",
                "cross_entropy":r"$CLE$", 
                "mean_link_cross_entropy":r"$CLE$", 
                "speedVeh":r"$\bar{S^v_r}$",
                'informalCrossing':'Informal Crossing'
                }

policy_param = list(policies.keys())[0]
policy_values = policies[policy_param]
scenario_param_cols =  [i for i in params if i!=policy_param]

problem = init_problem(params)

dfSIs = get_multiple_metrics_sis(dfDD, problem, policy_param, policy_values, outcome_vars3)
sobol_si_figure(dfSIs, gdfORLinks, policy_param, policy_values, outcome_vars3, rename_dict, inset_rec, constrained_layout = False, fig_width = 9, colors = ['#1b9e77', '#d95f02', '#7570b3'], ttc_threshold=ttc_threshold)



#
# Looking at metrics for different levels of vehicle flow
#

'''
group_param = 'avNVehicles'
policy_param = 'informalCrossing'
metric = 'speedVeh'
title = 'Vehicle speed increases with crossing restrictions'
f = agg_policy_comparison_figure(dfDD, gdfORLinks, group_param, policy_param, metric, rename_dict, inset_rec, title, colors = ['#1b9e77', '#d95f02', '#7570b3'], figsize=(20,10))

group_param = 'avNVehicles'
policy_param = 'informalCrossing'
metric = 'conflict_count'
title = 'Conflicts decrease with crossing restrictions'
f = agg_policy_comparison_figure(dfDD, gdfORLinks, group_param, policy_param, metric, rename_dict, inset_rec, title, colors = ['#1b9e77', '#d95f02', '#7570b3'], figsize=(20,10))
'''

group_param = 'avNVehicles'
policy_param = 'informalCrossing'
metrics = ['speedVeh','conflict_count']
title = 'Comparing vehicle speed and conflicts between policies'
agg_policy_two_metric_comparison_figure(dfDD, gdfORLinks, group_param, policy_param, metrics, rename_dict, inset_rec, title, colors = ['#1b9e77', '#d95f02', '#7570b3'], figsize = (16,10), quantile_groups = (0.25,0.5,0.75,1.0), quantile_labels = ("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"), ttc_threshold=ttc_threshold )