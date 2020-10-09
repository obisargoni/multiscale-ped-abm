from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from shapely.geometry import LineString
from shapely.geometry import Point
import csv
import matplotlib.pyplot as plt

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
    return filtered_files[0]



#####################################
#
# Globals
#
#####################################

data_dir = "..\\output\\batch\\model_run_data\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")
output_directory = "..\\output\\processed_data\\"

geo_data_dir = "..\\data\\model_gis_data\\"
ped_od_file = "OD_pedestrian_nodes.shp"
vehicle_polygons = "topographicAreaVehicle-withPriority-RoadLinkGrouping-Crossing.shp"


# Output paths
img_dir = "..\\output\\img\\"
path_deviation_fig = "path_deviation.png"
path_deviation_bar_fig = "path_deviation_bar.png"
crossing_percent_fig = "crossing_percent.png"

output_data_dir = "..\\output\\processed_data\\"

project_crs = {'init': 'epsg:27700'}

file_datetime_string = "2020.Aug.31.13_37_31"
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")

####################################
#
# Get pedestrian crossing choices
#
####################################

file_re = get_file_regex("pedestrian_locations", file_datetime = None)
ped_locations_file = most_recent_directory_file(data_dir, file_re)


df_ped_loc = pd.read_csv(os.path.join(data_dir, ped_locations_file))

# Select just the crossing choice column
desired_columns = ['ID', 'CrossingChoice','tick','run']
df_cc = df_ped_loc.reindex(columns = desired_columns)



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

df_ped_cc = df_cc.groupby(['run','ID'], group_keys=True)['CrossingChoice'].apply(get_peds_crossing_choice).reset_index()

# Get count of peds each run
ser_run_ped_counts = df_ped_cc.groupby('run')['ID'].apply(lambda s: s.unique().shape[0])
df_ped_counts = pd.DataFrame({'run_npeds':ser_run_ped_counts})

# Now group by run to get count of each crossing type
df_cc_count = df_ped_cc.groupby(['run','CrossingChoice']).count().unstack()
df_cc_count.columns = [c[1] for c in df_cc_count.columns]
df_cc_count.rename(columns = {'none':'undecided'}, inplace=True)

# Join to df of npeds per run and calculate percentages
df_cc_count = pd.merge(df_cc_count, df_ped_counts, left_index = True, right_index = True)

data_cols = [c for c in ['undecided','unmarked','unsignalised'] if c in df_cc_count.columns]
for c in data_cols:
    df_cc_count[c] = (df_cc_count[c] / df_cc_count['run_npeds']) * 100.0


###################################
#
# Now join data with run dataframe
#
###################################

file_re = get_file_regex("pedestrian_locations", file_datetime = None, suffix = 'batch_param_map')
batch_file = most_recent_directory_file(data_dir, file_re)
df_run = pd.read_csv(os.path.join(data_dir, batch_file))

'''
run_selection_dict =    {
                        "lambda":           [0.1,1],
                        "alpha":            [0.1,0.9],
                        "addVehicleTicks":  [10,50],
                        "epsilon":          [2, 6],
                        "gamma":            [0.1,0.9]
                        }

for col in run_selection_dict.keys():
    df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]
'''


# Join this to the data
df_cc_count = pd.merge(df_cc_count, df_run, left_index = True, right_on = 'run', how = 'inner')
df_cc_count.fillna(0, inplace = True)


#####################################
#
# Functions to make figures
#
#####################################
def plot_group_heatmap(df_group, row, col, value_col = 'cratio', alpha_col = 'uratio', cmap = plt.cm.viridis, ax = None, cbarlabel="Ratio of unmarked crossing to unsignalised crossing"):
    # get the data
    rgba, row_labels, col_labels = heatmap_rgba_data(df_group, row, col, value_col = value_col, alpha_col = alpha_col, cmap = cmap)

    im, cbar = heatmap(rgba, row_labels, col_labels, ax = ax, x_label = col, y_label = row, cbar_kw = dict(shrink = 0.7), cbarlabel=cbarlabel)

    return im, cbar

def heatmap_rgba_data(df_group, row, col, value_col = 'cratio', alpha_col = None, cmap = plt.cm.viridis):

    # Values used for pixel colours
    colour_data = df_group.reindex(columns = [row, col, value_col]).set_index([row, col]).unstack()

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

def heatmap(data, row_labels, col_labels, ax=None, x_label = None, y_label = None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    if x_label is not None:
        ax.set_xlabel(x_label)

    if y_label is not None:
        ax.set_ylabel(y_label)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.xaxis.set_label_position('top') 

    # Rotate the tick labels and set their alignment.
    #plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def batch_run_heatmap(df_data, groupby_columns, parameter_sweep_columns, rename_dict):

    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of 'addVehicleTicks'
    p = len(df_run[groupby_columns[0]].unique())
    q = len(df_run[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    # Set fixed indices for q and r axes
    qi = 0

    f,axs = plt.subplots(p, q,figsize=(20,10), sharey=False, sharex = False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        key_index = key_indices[pi, qi]
        group_key = keys[key_index]
        df_group = grouped.get_group(group_key)

        # Select the corresponding axis
        ax = axs[pi, qi]

        plot_group_heatmap(df_group, parameter_sweep_columns[0], parameter_sweep_columns[1], ax = ax)

    # Now add text annotations to indicate the scenario
    for i in range(p):
        ki = key_indices[i, 0]
        group_key = keys[ki]
        ax = axs[i, 0]

        s = None
        if group_key[0] == 10:
            s = "High\nVehicle\nFlow"
        elif group_key[0] == 50:
            s = "Low\nVehicle\nFlow"
        plt.text(-0.4,0.5, s, fontsize = 11, transform = ax.transAxes)
    '''
    for j in range(q):
        ki = key_indices[-1, j]
        group_key = keys[ki]

        ax = axs[-1, j]

        s = "{}".format(rename_dict[group_key[1]])
        plt.text(1.05,1.1, s, fontsize = 11, transform = ax.transAxes)
    '''

    return f, axs


#####################################
#
#
# Plot specific data prep
#
#
#####################################
rename_dict = {'addVehicleTicks':"Ticks\nBetween\nVehicle\nAddition",'alpha':r"$\mathrm{\alpha}$",'lambda':r"$\mathrm{\lambda}$"}

# Values of interest are the ratio of unmarked crossings to unsignalised crossings and the fractio of pedestrians that are undecided
def crossing_ratio(row, c1 = 'unmarked', c2 = 'unsignalised'):
    if row[c2] == 0:
        return 0.0
    else:
        return row[c1] / row[c2]

df_cc_count['cratio'] = df_cc_count.apply(lambda r: crossing_ratio(r), axis = 1)
df_cc_count['uratio'] = 1 - df_cc_count['undecided'] / 100.0

# Groups by the variables I want to keep constant in eac plot
groupby_columns = ['addVehicleTicks', 'epsilon']
parameter_sweep_columns = ['alpha', 'lambda']

f, axs = batch_run_heatmap(df_cc_count, groupby_columns, parameter_sweep_columns, rename_dict)

f.show()