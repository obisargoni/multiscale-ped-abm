from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from shapely.geometry import LineString
from shapely.geometry import Point
import csv

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

file_datetime_string = "2020.Feb.20.14_58_33"
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")

####################################
#
# Get pedestrian crossing choices
#
####################################

file_re = get_file_regex("pedestrian_locations")
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


# Now group by run to get count of each crossing type
df_cc_count = df_ped_cc.groupby(['run','CrossingChoice']).count().unstack()
df_cc_count.columns = [c[1] for c in df_cc_count.columns]
df_cc_count.rename(columns = {'none':'undecided'}, inplace=True)

# Multiply by 10 so values are percentages
df_cc_count = df_cc_count*10


###################################
#
# Now join data with run dataframe
#
###################################

file_re = get_file_regex("pedestrian_locations", suffix = 'batch_param_map')
batch_file = most_recent_directory_file(data_dir, file_re)
df_run = pd.read_csv(os.path.join(data_dir, batch_file))

selection_columns = ['lambda', 'alpha', 'addVehicleTicks']
selction_values = [ [0.1,1],
                    [0.1,0.9],
                    [5,20]
                    ]

run_selection_dict = {selection_columns[i]:selction_values[i] for i in range(len(selection_columns))}

for col in selection_columns:
    df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]


# Join this to the data
df_cc_count = pd.merge(df_cc_count, df_run, left_index = True, right_on = 'run', how = 'inner')



#####################################
#
# Make fgures
#
#####################################
import matplotlib.pyplot as plt


def batch_run_map(df_data, run_col, rename_dict, title, output_path):

    groupby_columns = ['addVehicleTicks','alpha','lambda']
    data_columns = ['undecided','unmarked','unsignalised']
    colors = ['grey','blue', 'red']

    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    p = len(run_selection_dict[groupby_columns[0]])
    q = len(run_selection_dict[groupby_columns[1]])
    r = len(run_selection_dict[groupby_columns[2]])
    key_indices = np.reshape(np.arange(len(keys)), (p,q*r))

    fig_indices = np.reshape(key_indices, (p,q*r))
    f,axs = plt.subplots(p, q*r,figsize=(20,10), sharey=True, sharex = True)
    for ki in range(len(keys)):
        group_key = keys[ki]
        data = grouped.get_group(group_key)

        # Check have got a single runs data
        assert data[run_col].unique().shape[0] == 1

        # if one of the data columns not present initialise it as 0
        for c in [i for i in data_columns if i not in data.columns]:
            data[c] = 0

        # Get axis
        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        # Plot bar chart
        tick_labels = data_columns
        x = range(len(tick_labels))
        y = data[data_columns].fillna(0).iloc[0].values

        ax = ax.bar(x, y, tick_label = tick_labels, color = colors)


    # Add in row group info to the last axis in the row
    for i in range(p):
        ki = key_indices[i, q*r - 1]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}:\n{}".format(rename_dict[groupby_columns[0]], group_key[0])
        plt.text(1.1,0.5, s, fontsize = 15, transform = ax.transAxes)

    for j in range(q):
        ki = key_indices[0, j*r]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}: {}".format(rename_dict[groupby_columns[1]], group_key[1])
        plt.text(1.05,1.1, s, fontsize = 11, transform = ax.transAxes)

        # Add some explanitory text
        t = ""
        if group_key[1] == 0.1:
            t = "Sensitive to traffic"
        elif group_key[1] == 0.9:
            t = "Sensitive to journey time"
        plt.text(0.95,1.2, t, fontsize = 15, transform = ax.transAxes)

    for k in range(q*r):
        ki = key_indices[1, k]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}: {}".format(rename_dict[groupby_columns[2]], group_key[2])
        plt.text(0.45, -0.15, s, fontsize = 11, transform = ax.transAxes)

        # Add some explanitory text
        t = ""
        if group_key[2] == 0.1:
            t = "Plans ahead"
        elif group_key[2] == 1:
            t = "Considers nearby\nalternatives more"
        plt.text(0.35,-0.3, t, fontsize = 15, transform = ax.transAxes)


    f.suptitle(title, fontsize=16, y = 1)
    f.show()
    plt.savefig(output_path)

map_output_path = "..\\output\\img\\crossing_choice_bar.png"
rename_dict = {'addVehicleTicks':"Ticks\nBetween\nVehicle\nAddition",'alpha':r"$\mathrm{\alpha}$",'lambda':r"$\mathrm{\lambda}$"}
batch_run_map(df_cc_count, 'run', rename_dict, "Crossing Choices Percentages", map_output_path)