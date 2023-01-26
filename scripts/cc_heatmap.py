import os
import json
import re
from datetime import datetime as dt
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from colorspacious import cspace_convert

from importlib import reload
import batch_data_utils
batch_data_utils = reload(batch_data_utils)
from batch_data_utils import get_processed_crossing_locations_data, get_data_paths, load_and_clean_cross_events


#####################################
#
# Functions to make figures
#
#####################################
def heatmap_rgba_data(df_group, row, col, value_col = 'unmarked_pcnt', alpha_col = None, cmap = plt.cm.viridis):

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

def heatmap(data, row_labels, col_labels, ax = None, x_label = None, y_label = None, **kwargs):
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

    return im

def annotate_heatmap(im, data=None, value_data = None, valfmt="{x:.0f}", textcolors=["white","black"], threshold=None, exclude = [], **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    exclude
        Values to exclude from annotating. Optional
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i,j] in exclude:
                continue
            else:
                kw.update(color=textcolors[int(30 < value_data[i,j] < 70)])
                text = im.axes.text(j, i, 'u', **kw)
                texts.append(text)

    return texts

def batch_run_heatmap(df_data, groupby_columns, parameter_sweep_columns, value_col, alpha_col, annotate_col, rename_dict, cmap, title = None, cbarlabel = None, output_path = None):

    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    # Want to get separate array of data for each value of 'addVehicleTicks'
    p = len(df_data[groupby_columns[0]].unique())
    q = len(df_data[groupby_columns[1]].unique())

    key_indices = np.reshape(np.arange(len(keys)), (p,q))

    f,axs = plt.subplots(p, q, figsize=(13,10), sharey=False, sharex=False)

    # Make sure axes array in shame that matches the layout
    axs = np.reshape(axs, (p, q))

    # Select data to work with and corresponding axis
    for pi in range(p):
        for qi in range(q):
            key_index = key_indices[pi, qi]
            group_key = keys[key_index]
            df_group = grouped.get_group(group_key)

            # Select the corresponding axis
            ax = axs[pi, qi]

            # get the data and plot image
            df_group[value_col] = df_group[value_col].fillna(0.0)
            rgba, row_labels, col_labels = heatmap_rgba_data(df_group, parameter_sweep_columns[0], parameter_sweep_columns[1], value_col = value_col, alpha_col = alpha_col, cmap = cmap)
            im = heatmap(rgba, row_labels, col_labels, ax = ax,  y_label = rename_dict[parameter_sweep_columns[0]], x_label = rename_dict[parameter_sweep_columns[1]], cmap = cmap)

            # Annotate heatmap
            if annotate_col is not None:
                annotate_data = df_group.reindex(columns = [parameter_sweep_columns[0], parameter_sweep_columns[1], annotate_col]).set_index([parameter_sweep_columns[0], parameter_sweep_columns[1]]).unstack().values
                colour_data = df_group.reindex(columns = [parameter_sweep_columns[0], parameter_sweep_columns[1], value_col]).set_index([parameter_sweep_columns[0], parameter_sweep_columns[1]]).unstack().values
                texts = annotate_heatmap(im, data=annotate_data, value_data = colour_data, exclude = [0, 40], fontweight='bold')


    # Adjust the plot to make space for the colourbar axis
    plt.subplots_adjust(right=0.8, wspace = 0.1)

    # Create new axis at far right of plot - [left, bottom, width, height]
    cax = f.add_axes([0.82, 0.2, 0.03, 0.6])
    # Create colorbar
    cbar = f.colorbar(im, cax=cax, anchor = (0,0.7), cmap = cmap)

    # Now add text annotations to indicate the scenario
    for i in range(p):
        ki = key_indices[i, 0]
        group_key = keys[ki]
        ax = axs[i, 0]

        s = "{}".format(rename_dict[group_key[0]])
        plt.text(-0.4,0.5, s, fontsize = 11, transform = ax.transAxes)


    for j in range(q):
        ki = key_indices[-1, j]
        group_key = keys[ki]

        ax = axs[-1, j]

        s = "{}".format(rename_dict[group_key[1]])
        plt.text(0.3,-0.1, s, fontsize = 11, transform = ax.transAxes)

    if title is not None:
        f.suptitle(title, fontsize=16, y = 1)
    if cbarlabel is not None:
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    if output_path is not None:
        plt.savefig(output_path)

    return f, axs


#####################################
#
#
# Get data from alpha and lambda param sweep runs for each configuration
#
#
#####################################
with open(".//config.json") as f:
    config = json.load(f)

rename_dict = { "avNVehicles": r"$\bar{N^v}_r$",
                'alpha':r"$\mathrm{\alpha}$",
                'lambda':r"$\mathrm{\lambda}$",
                "epsilon":r"$\mathrm{\epsilon}$",
                "gamma":r"$\mathrm{\gamma}$",
                "between": "Between Configuration",
                "beyond":"Beyond Configuration",
                15:"High\nVehicle\nFlow",
                1:"Low\nVehicle\nFlow"
                }

configuration_datetime_strings = {
                                    "between":  dt.strptime(config['al_sweep_between'], "%Y.%b.%d.%H_%M_%S"),
                                    "beyond":   dt.strptime(config['al_sweep_beyond'], "%Y.%b.%d.%H_%M_%S")
                                }

configuration_datetime_strings = {
                                    "between":  config['al_sweep_between'],
                                    "beyond":   config['al_sweep_beyond']
                                }


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
crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])

# Model output data
between_data_paths = get_data_paths(configuration_datetime_strings['between'], data_dir)
between_cross_events_file = between_data_paths["cross_events"]
between_batch_file = between_data_paths["batch_file"]

beyond_data_paths = get_data_paths(configuration_datetime_strings['beyond'], data_dir)
beyond_cross_events_file = beyond_data_paths["cross_events"]
beyond_batch_file = beyond_data_paths["batch_file"]

dfRunBetween = pd.read_csv(os.path.join(data_dir, between_batch_file))
dfRunBeyond = pd.read_csv(os.path.join(data_dir, beyond_batch_file))
assert (dfRunBetween==dfRunBeyond).all().all()


# GIS Data
gdfPaveLinks = gpd.read_file(pavement_links_file)
gdfPaveNodes = gpd.read_file(pavement_nodes_file)
gdfCAs = gpd.read_file(crossing_alternatives_file)

dfCrossEventsBetween = load_and_clean_cross_events(gdfPaveLinks, cross_events_path = between_cross_events_file)
dfCrossEventsBeyond = load_and_clean_cross_events(gdfPaveLinks, cross_events_path = beyond_cross_events_file)

# Aggregate to get percent informal or formal crossing
dfCrossTypesBetween = dfCrossEventsBetween.groupby('run')['CrossingType'].apply(lambda s: (s.value_counts() / s.value_counts().sum()) *100 ).unstack().reset_index()
dfCrossTypesBeyond = dfCrossEventsBeyond.groupby('run')['CrossingType'].apply(lambda s: (s.value_counts() / s.value_counts().sum()) *100 ).unstack().reset_index()

 # Merge to batch run params files
dfCrossTypesBetween = pd.merge(dfRunBetween, dfCrossTypesBetween, on='run')
dfCrossTypesBeyond = pd.merge(dfRunBetween, dfCrossTypesBeyond, on='run')
dfCrossTypesBetween['configuration'] = 'between'
dfCrossTypesBeyond['configuration'] = 'beyond'
df_cc_count_al = pd.concat([dfCrossTypesBetween, dfCrossTypesBeyond])


# This old method accounted for agents being undecided. Also based on old data outputs
'''
btwn_ped_cc = get_processed_crossing_locations_data(data_dir, "pedestrian_locations", configuration_datetime_strings['between'])
btwn_ped_cc["configuration"] = "between"

bynd_ped_cc = get_processed_crossing_locations_data(data_dir, "pedestrian_locations", configuration_datetime_strings['beyond'])
bynd_ped_cc["configuration"] = "beyond"

df_cc_count_al = pd.concat([btwn_ped_cc, bynd_ped_cc])
'''

# Groups by the variables I want to keep constant in eac plot
groupby_columns = ['avNVehicles', 'configuration']
parameter_sweep_columns = ['alpha', 'lambda']

fig_title = "Crossing Choices\n{} and {} parameter sweep".format(r"$\mathrm{\alpha}$", r"$\mathrm{\lambda}$")
fig_path = os.path.join(img_dir, "al_crossing_heatmap_{}.png".format(configuration_datetime_strings['between']))

f, axs = batch_run_heatmap(df_cc_count_al, groupby_columns, parameter_sweep_columns, 'unmarked', None, None, rename_dict, title = fig_title, cbarlabel = "Proportion choosing informal crossings", cmap = plt.cm.coolwarm_r, output_path = fig_path)
f.show()


#####################################
#
#
# Get data from epsilon and gamma parameter sweep runs for each configuration
#
#
#####################################
configuration_datetime_strings = {
                                    "between":  dt.strptime(config['eg_sweep_between'], "%Y.%b.%d.%H_%M_%S"),
                                    "beyond":   dt.strptime(config['eg_sweep_beyond'], "%Y.%b.%d.%H_%M_%S")
                                }

btwn_ped_cc = get_processed_crossing_locations_data(data_dir, "pedestrian_locations", configuration_datetime_strings['between'])
btwn_ped_cc["configuration"] = "between"

bynd_ped_cc = get_processed_crossing_locations_data(data_dir, "pedestrian_locations", configuration_datetime_strings['beyond'])
bynd_ped_cc["configuration"] = "beyond"

df_cc_count_eg = pd.concat([btwn_ped_cc, bynd_ped_cc])

# If zero pedestrians have crossed this gives nan. Replace with zero. Need to use alpha instead to represent missing data (or mask array...?)
df_cc_count_eg.fillna(0, inplace = True)

# Create col that is zero if all peds undecided and one otherwise. Use this to set opacity of the plot
df_cc_count_eg['opacity'] = 1 - np.floor(df_cc_count_eg['undecided_frac'])

# Groups by the variables I want to keep constant in eac plot
groupby_columns = ['addVehicleTicks', 'configuration']
parameter_sweep_columns = ['epsilon', 'gamma']

fig_title = "Crossing Choices\n{} and {} parameter sweep".format(r"$\mathrm{\epsilon}$", "$\mathrm{\gamma}$")
fig_path = os.path.join(img_dir, "eg_crossing_heatmap_{}.png".format(configuration_datetime_strings['between'].strftime("%Y.%b.%d.%H_%M_%S")))

# 'inverse_undecided_frac'

f, axs = batch_run_heatmap(df_cc_count_eg, groupby_columns, parameter_sweep_columns, 'unmarked_pcnt', 'opacity', 'undecided', rename_dict, title = fig_title, cbarlabel = "Proportion choosing informal crossings", cmap = plt.cm.coolwarm_r, output_path = fig_path)
f.show()
