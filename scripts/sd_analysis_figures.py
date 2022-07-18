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

def pair_plot(dfDD, outcome_vars, policy_col, rename_dict, file_datetime_string):
    data = dfDD.loc[:, outcome_vars].rename(columns = rename_dict)
    data[rename_dict[policy_col]] = dfDD[policy_col]
    sns.pairplot(data, hue=rename_dict[policy_col], vars=[rename_dict[outcome_vars[0]], rename_dict[outcome_vars[1]]])
    plt.savefig(os.path.join(img_dir, 'pair_plot.{}-{}.{}.png'.format(outcome_vars[0], outcome_vars[1], file_datetime_string)))

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

        ax.hist(df[val_col].dropna(), bins = bins, density=False, color=palette[i], alpha=1, label=g, histtype='step', linewidth=6)

    ax.legend(fontsize = 20)
    ax.set_ylabel('count', fontsize = 20)
    #ax.set_xlim(xmin=minv*0.9, xmax=maxv*1.1)
    ax.set_title(title, fontsize = 24)
    return ax

def multi_hist_plot(dfDD, gdfORLinks, outcome_vars, policy_col, title_rename_dict, fig_config, inset_rec, nhistbins = 25, figsize=(20,10)):
    nvars = len(outcome_vars)
    fig, axs = plt.subplots(1, nvars, figsize=figsize)

    data = dfDD.loc[:, outcome_vars+[policy_col]]
    for i in range(nvars):
        ax0 = hist_plot(axs[i], data, outcome_vars[i], policy_col, title_rename_dict[outcome_vars[i]], nhistbins = nhistbins)

    # add inset showing the road network
    axins = fig.add_axes(inset_rec)
    gdfORLinks.plot(ax=axins, color='black')
    axins.set_axis_off()
    axins.set_title('Environment', y=-0.2)

    outpath = os.path.join(img_dir, 'hists.{}.png'.format(file_datetime_string))
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
nbins = None
bin_dist = 2

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

output_paths = bd_utils.get_ouput_paths(file_datetime_string, vehicle_density_timestamp, data_dir, nbins = bin_dist)
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
outcome_vars3 = ['DistPA','crossCountPP','cross_entropy']
policy_col = 'informalCrossing'

title_rename_dict = {   "route_length_pp":r"$\bar{L_r}$",
                        "DistPA": r"$\bar{D_r}$",
                        "crossCountPP":r"$\bar{C_r}$",
                        "cross_entropy":r"$CLE$", 
                        'informalCrossing':'Informal Crossing'}
#
# Create pairs plot
#
pair_plot(dfDD, outcome_vars1, policy_col, title_rename_dict, file_datetime_string)
pair_plot(dfDD, outcome_vars2, policy_col, title_rename_dict, file_datetime_string)

#
# Histogram plots
#
#plt.style.use('dark_background')
inset_rec = [0, 0.85, 0.13, 0.13]
multi_hist_plot(dfDD, gdfORLinks, outcome_vars3, policy_col, title_rename_dict, fig_config, inset_rec, figsize=(10*len(outcome_vars3),10))

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