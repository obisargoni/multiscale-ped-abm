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

#from sklearn.model_selection import train_test_split # for splitting the data into train and test samples
from sklearn.metrics import classification_report # for model evaluation metrics
from sklearn import tree, linear_model # for decision tree models
from sklearn.feature_extraction import DictVectorizer

import seaborn as sns
import matplotlib.patches as mpatches

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

def hist_plot(ax, data, val_col, val_unit, group_col, title, nhistbins = 25, palette=['#3d993d', '#cca328', '#af3f33']):

    minv = data[val_col].min()
    maxv = data[val_col].max()
    rng = maxv-minv
    bins = np.linspace(minv, maxv, nhistbins)

    groups = list(data[group_col].unique())
    groups.sort(key=sf)
    for i, g in enumerate(groups):
        df = data.loc[ data[group_col]==g]

        ax.hist(df[val_col].dropna(), bins = bins, density=False, color=palette[i], alpha=0.7, label=g, histtype='step', linewidth=6)

    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylim(ymin=0, ymax=900)
    ax.set_xlabel(val_unit, fontsize = 23, labelpad=15)
    ax.set_title(title, fontsize = 40, pad=25)
    return ax

def multi_env_hist_plot(f, axs, dfDD, gdfORLinks, outcome_vars, outcome_units, env_col, title_rename_dict, fig_config, nhistbins = 25, palette = ['#3d993d', '#cca328', '#af3f33']):

    data = dfDD.loc[:, outcome_vars+[env_col]]
    for i in range(nvars):
        ax0 = hist_plot(axs[i], data, outcome_vars[i], outcome_units[i], env_col, title_rename_dict[outcome_vars[i]], nhistbins = nhistbins, palette = palette)
        axs[i].set_yticks(range(0, 1000, 200))
        axs[i].tick_params(axis='y', length=10)
        if i==0:
            axs[i].set_ylabel('count of runs', fontsize = 35, labelpad=17)
            axs[i].set_yticklabels(range(0, 1000, 200))
        axs[i].tick_params(axis = 'both', labelsize = 26)

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
                print("Null values in utput for {}, skipping".format(metric))
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

def multi_env_sobol_si_plot(f, axs, dfSIs, gdfORLinks, env_col, env_values, outcome_vars, rename_dict, params = ['alpha','lambda','epsilon','timeThreshold', 'addPedTicks','addVehicleTicks', 'tacticalPlanHorizon', 'minCrossing'], constrained_layout = True, colors = ['#3d993d', '#cca328', '#af3f33']):

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
            data = dfsip.set_index('param').loc[params][['ST','ST_conf']]
        
            x_pos = np.arange(data.shape[0])+xi[j]
            axs[i].bar(x_pos, data['ST'], width=bar_width, yerr = data['ST_conf'], align='center', label=env, color = colors[j])

        axs[i].set_yticks([])
        axs[i].set_yticklabels([])
        axs[i].set_ylim(ymin=0, ymax=1)
        axs[i].set_xticks(np.arange(data.shape[0]))
        axs[i].set_xticklabels([ rename_dict[i] for i in data.index], rotation=45)

        axs[i].set_yticks(np.round(np.linspace(0, 1, 6),1))
        axs[i].tick_params(axis='y', length=10)
        if i==0:
            axs[i].set_ylabel(r"$S_T$", fontsize=35, labelpad=17)
            axs[i].set_yticklabels( np.round(np.linspace(0, 1, 6),1))

        axs[i].tick_params(axis = 'y', labelsize = 26)
        axs[i].tick_params(axis = 'x', labelsize = 35)

        #axs[i].set_title(title, fontsize=24)
    axs[-1].legend(fontsize = 35, bbox_to_anchor=(-2.37,-0.45,2.5,0), loc="lower center", mode='expand', ncol=3) # use bbox_to_anchor = [-2,-0.32,2.5,0] for 3 variables
    #bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2

    return f

def fitting_tree(X, y, feature_names, criterion, splitter, mdepth, clweight, minleaf):

    # Fit the model
    model = tree.DecisionTreeRegressor(criterion=criterion, 
                                        splitter=splitter, 
                                        max_depth=mdepth,
                                        min_samples_leaf=minleaf, 
                                        random_state=0, 
                                  )
    clf = model.fit(X, y)

    # Predict class labels on training data
    pred_labels_tr = model.predict(X)

    # Tree summary and model evaluation metrics
    print('*************** Tree Summary ***************')
    print('Tree Depth: ', clf.tree_.max_depth)
    print('No. of leaves: ', clf.tree_.n_leaves)
    print('No. of features: ', clf.n_features_in_)
    print('--------------------------------------------------------')
    print("")
    
    print('*************** Evaluation on Training Data ***************')
    score_tr = model.score(X, y)
    print('Accuracy Score: ', score_tr)
    print('--------------------------------------------------------')
    
    # Use graphviz to plot the tree
    dot_data = tree.export_graphviz(clf, out_file=None, 
                                feature_names=feature_names, 
                                class_names=None,
                                filled=True, 
                                rounded=True, 
                                #rotate=True,
                               ) 
    graph = graphviz.Source(dot_data)
    
    # Return relevant data for chart plotting
    return X, y, clf, graph

def fitting_linear_regression(X, y, feature_names):

    # Fit the model
    model = linear_model.LinearRegression( fit_intercept=True)
    clf = model.fit(X, y)
    
    print('*************** Evaluation on Training Data ***************')
    score_tr = model.score(X, y)
    print('Accuracy Score: ', score_tr)
    print('--------------------------------------------------------')
    
    # Return relevant data for chart plotting
    return X, y, clf

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


ug_output_paths = bd_utils.get_ouput_paths(config['ug_results'], data_dir, nbins = bin_dist)
ug_output_path = ug_output_paths["output_route_data"]

qg_output_paths = bd_utils.get_ouput_paths(config['qg_results'], data_dir, nbins = bin_dist)
qg_output_path = qg_output_paths["output_route_data"]

cc_output_paths = bd_utils.get_ouput_paths(config['cc_results'], data_dir, nbins = bin_dist)
cc_output_path = cc_output_paths["output_route_data"]


output_paths = { 'Uniform Grid':    ug_output_path,
                'Quad Grid':        qg_output_path,
                'Clapham Common':   cc_output_path}

'''
oneway_output_paths = bd_utils.get_ouput_paths(config['one_way_results'], data_dir, nbins = bin_dist)
oneway_output_path = oneway_output_paths["output_route_data"]

twoway_output_paths = bd_utils.get_ouput_paths(config['two_way_results'], data_dir, nbins = bin_dist)
twoway_output_path = twoway_output_paths["output_route_data"]

output_paths = {'One Way Peds': oneway_output_path,
                'Two Way Peds': twoway_output_path}
'''

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
outcome_vars = ['route_length_pp', 'sp_sim_zerocount_pct', 'PostponeCountPP', 'pcntInfCross']
outcome_units = ['meters', '%', 'crossings', '%']

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
                'informalCrossing':'Informal Crossing',
                'sp_sim': r"$\bar{\Delta^{SP}_r}$",
                'sp_sim_dice': r"$\Delta SP_{dice}$",
                'sp_sim_zerocount': r"$N^{SP}_r$",
                'sp_sim_zerocount_pct':r'$SP_r',
                'PostponeCountPP': r"$\bar{P}_r$",
                'pcntInfCross': r"$I_r$"
                }

name_titles = { 'route_length_pp': 'Av. Route Length, '+ r"$\bar{L_r}$",
                'sp_sim_zerocount': 'N Shortest Paths, ' + r"$N^{SP}_r$",
                'sp_sim_zerocount_pct': 'Pcnt. Shortest Paths, ' + r"$SP_r$",
                'PostponeCountPP':'Av. Crossing Postponements, ' + r"$\bar{P}_r$",
                'pcntInfCross': 'Pcnt. Informal Crossings, ' + r"$I_r$"}

name_titles_pair = { 'route_length_pp': 'Av. Route\nLength,\n'+ r"$\bar{L_r}$",
                'sp_sim_zerocount': 'N Shortest\nPaths,\n' + r"$N^{SP}_r$",
                'sp_sim_zerocount_pct': 'Pcnt. Shortest\nPaths,\n' + r"$SP_r$",
                'PostponeCountPP':'Av. Crossing\nPostponements,\n' + r"$\bar{P}_r$",
                'pcntInfCross': 'Pcnt. Informal\nCrossings,\n' + r"$I_r$"}

# update the rename dict with name titles - doesn't look good.
for k,v in name_titles_pair.items():
    rename_dict[k] = v

#
# Initialise figure
#
nplots = len(plot_types)
nvars = len(outcome_vars)
fig_width = 10
fig, axs = plt.subplots(nplots, nvars, figsize=(fig_width*len(outcome_vars),17))
plt.subplots_adjust(left=0.05,
                    bottom=0.15, 
                    right=0.98, 
                    top=0.93, 
                    wspace=0.08, 
                    hspace=0.23)

#
# Histogram plots
#
#plt.style.use('dark_background')
fig = multi_env_hist_plot(fig, axs[0, :], dfDD, gdfORLinks, outcome_vars, outcome_units, env_col, name_titles, fig_config, palette=palette)

#
# Sobol indices for each metric and policy setting
#
from SALib.analyze import sobol

policy_param = list(policies.keys())[0]
scenario_param_cols =  [i for i in params if i!=policy_param]

env_values = ['Uniform Grid','Quad Grid','Clapham Common']
#env_values = ['One Way Peds','Two Way Peds']

problem = init_problem(params)

dfSIs = get_multiple_metrics_sis_envs(dfDD, problem, env_col, outcome_vars)
fig = multi_env_sobol_si_plot(fig, axs[1, :], dfSIs, gdfORLinks, env_col, env_values, outcome_vars, rename_dict, constrained_layout = False, colors = palette)

# annotate figure
texts = ['i','ii','iii','iv','v','vi', 'vii','viii','ix']
for i, ax in enumerate(axs.reshape(1,-1)[0]):
    ax.text(0.945, 0.935, texts[i], transform=ax.transAxes, fontsize=20, fontweight='bold')

fig.savefig(os.path.join(img_dir, "env_comparison_{}_{}_{}.png".format(config['ug_results'], config['qg_results'], config['cc_results'])))


#
# Regression Tree plots
#
'''
import graphviz # for plotting decision tree graphs

# Tree settings
criterion = 'squared_error'
splitter = 'best'
mdepth = 4
clweight = None

minleaf = 10
model_results = {}
grouped = dfDD.groupby(env_col)
for i, (env) in enumerate(grouped.groups.keys()):
    model_results[env] = {}
    for var in outcome_vars:
        df_env = grouped.get_group((env))
        
        data_dict = df_env.loc[:, problem['names']].to_dict('records')
        vec = DictVectorizer()  # create the DictVectorizer object
        vec_array = vec.fit_transform(data_dict).toarray()  # execute process on the record dictionaries and transform the result into a numpy array object
        
        X, y, clf, graph = fitting_tree(vec_array, df_env[var].values, vec.get_feature_names_out(), criterion, splitter, mdepth, clweight, minleaf)
        model_results[env][var] = {'X':X, 'y':y, 'clf':clf, 'graph':graph}

        # save figure
        graph.format='png'
        graph.render(filename="tree.{}.{}".format(env, var), directory=img_dir)
'''

#
# Pair plot figure
#
dfDDPair = dfDD.rename(columns = rename_dict)
x_vars = [rename_dict[i] for i in ['alpha','lambda','epsilon','addVehicleTicks', 'tacticalPlanHorizon', 'minCrossing']]
y_vars = [rename_dict[i] for i in outcome_vars]
grid = sns.pairplot(dfDDPair, hue=env_col, palette = palette, x_vars = x_vars, y_vars = y_vars, kind = 'reg', diag_kind = 'hist', height = 2.5, plot_kws=dict(scatter_kws=dict(s=0.8, alpha=0.1)))

# remove default legend
grid._legend.remove()

# Create patches for the legend
ug_patch = mpatches.Patch(color=palette[0], label='Uniform Grid')
qg_patch = mpatches.Patch(color=palette[1], label='Quad Grid')
cc_patch = mpatches.Patch(color=palette[2], label='Clapham Common')
grid.axes[-1,-1].legend(handles = [ug_patch, qg_patch, cc_patch], fontsize = 16, bbox_to_anchor=(-4,-0.55,3.8,0), loc="lower center", mode='expand', ncol=3)

# Set font sizes
for ax in grid.axes.reshape(1,-1)[0]:
    if len(ax.get_xticklabels()) != 0:
        ax.set_xticklabels(ax.get_xticklabels(), fontdict = dict(fontsize=11))
        ax.set_xlabel(ax.get_xlabel(), fontdict = dict(fontsize=15))

    if len(ax.get_yticklabels()) != 0:
        ax.set_yticklabels(ax.get_yticklabels(), fontdict = dict(fontsize=11))
        ax.set_ylabel(ax.get_ylabel(), fontdict = dict(fontsize=15))

# Use annotations to label plots with metrics
for i, ax in enumerate(grid.axes[:,0]):
    metric_title = ax.yaxis.get_label_text()
    ax.text(-0.65, 0.5, metric_title, transform=ax.transAxes, fontsize=15, horizontalalignment='center',verticalalignment='center')
    ax.yaxis.set_label_text(outcome_units[i], fontdict = dict(fontsize=11))

# Edit x ticks for a couple of plots

# Column of MC axes
for i, ax in enumerate(grid.axes[:,-1]):
    ax.set_xticks([0,1])
    if i==3:
        ax.set_xticklabels(['false', 'true'])

# Column of T_veh axes
for i, ax in enumerate(grid.axes[:,-3]):
    if i==3:
        t = [0,1,2,3,4]
        l = [5 * (2**i) for i in t]
        ax.set_xticks(t)
        ax.set_xticklabels(l)

# Change lambda tick labels
grid.axes[-1,1].set_xticks([0,1.0,2.0])
grid.axes[-1,1].set_xticklabels(['0','1.0','2.0'])

grid.axes[-1,0].set_xticks([0,0.5,1.0])
grid.axes[-1,0].set_xticklabels(['0','0.5','1.0'])


grid.savefig(os.path.join(img_dir, 'pair_plot.{}.{}.{}.png'.format(config['ug_results'], config['qg_results'], config['cc_results'])))

#
# Print some useful summary stats
#

print("\nDescribe percentage of peds following shortest path")
print(dfDD.groupby(env_col)['sp_sim_zerocount_pct'].describe())

print("\nDescribe different in length to shortest path")
print(dfDD.groupby(env_col)['sp_sim'].describe())

print("\nDescribe number of crossings")
print(dfDD.groupby(env_col)['crossCountPP'].describe())

print("\nDescribe percent informal crossing")
print(dfDD.groupby(env_col)['pcntInfCross'].describe())

print("\nDescribe postpone crossing - 70% quantile")
print(dfDD.groupby(env_col)['PostponeCountPP'].quantile(0.7))

