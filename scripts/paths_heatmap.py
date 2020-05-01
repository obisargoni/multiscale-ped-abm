#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
from datetime import datetime as dt

import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

#from cartoframes.auth import set_default_credentials
#from cartoframes.viz import Map, Layer, animation_widget, animation_style, basemaps

#set_default_credentials('cartoframes')


# In[2]:


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

def most_recent_directory_file(directory, file_regex):
    files = os.listdir(directory)
    filtered_files = [f for f in files if file_regex.search(f) is not None]
    filtered_files.sort(key = lambda x: dt_from_file_name(x, file_regex), reverse=True)
    return filtered_files[0]

def dt_from_file_name(file_name, regex):
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt


###############################
#
#
# Read in pedestrian locations data
#
#
###############################
data_dir = "..\\output\\batch\\model_run_data\\"
outpath = "C:\\Users\\obisargoni\\eclipse-workspace\\repastInterSim\\output\\hex_bin_trajectories\\hex_bin_trajectories.shp"

project_crs = {'init': 'epsg:27700'}
wsg_crs = {'init':'epsg:4326'}

file_re = get_file_regex("pedestrian_locations")
ped_locations_file = most_recent_directory_file(data_dir, file_re)

df_ped_loc = pd.read_csv(os.path.join(data_dir, ped_locations_file))

# Now split df into one for locations, one for route coords
lo_cols = [i for i in df_ped_loc.columns if 'Loc' in i]

df_loc = df_ped_loc.dropna(subset = lo_cols)
for c in lo_cols:
    df_loc[c] = df_loc[c].map(float)


gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.LocXString, df_loc.LocYString))
gdf_loc.crs = project_crs

gdf_loc = gdf_loc.to_crs(wsg_crs)


###############################
#
#
# Read in batch runs metadata
#
#
################################
file_re = get_file_regex("pedestrian_initial_route", suffix = 'batch_param_map')
batch_file = most_recent_directory_file(data_dir, file_re)
df_run = pd.read_csv(os.path.join(data_dir, batch_file))

selection_columns = ['vehiclePriorityCostRatio', 'addVehicleTicks', 'cellCostUpdate']
selction_values = [ [1,100],
                    [5, 50],
                    [0, 100]
                    ]

run_selection_dict = {selection_columns[i]:selction_values[i] for i in range(len(selection_columns))}

for col in selection_columns:
    df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]


#################################
#
#
# Load the hexagonal bins and join with pedestrian trajectories to get the count per bin. 
# Once processed save this as a new version of the hexagonal bins
#
#
#################################
'''

# Use hexagonal tiles to bin pedestrian locations, create heat map of trajectories
hex_polys_file = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\hexgrid.shp"
gdf_hex = gpd.read_file(hex_polys_file)
print(gdf_hex.crs)

gdf_hex = gdf_hex.to_crs(wsg_crs)
gdf_hex['hex_id'] = np.arange(gdf_hex.shape[0])
gdf_hex = gdf_hex.reindex(columns = ['hex_id','geometry'])
print(gdf_hex.head())

# spatial join hex polys to the pedestrian locations
gdf_hex_traj = gpd.sjoin(gdf_hex, gdf_loc, op = 'contains')

# Get point count per hex per run
ser_hex_counts = gdf_hex_traj.groupby(["hex_id", "run"])['index_right'].apply(lambda df: df.shape[0])
df_hex_counts = pd.DataFrame(ser_hex_counts).reset_index().rename(columns = {'index_right':'loc_count'})

# Join back with the polygons and with df run to get parameter info

df_hex_counts = pd.merge(df_hex_counts, gdf_hex, on = 'hex_id', how = 'left')
df_hex_counts = pd.merge(df_hex_counts, df_run, on = 'run', how = 'inner') # Only include the runs we are interested in

gdf_hex_counts = gpd.GeoDataFrame(df_hex_counts)
gdf_hex_counts.crs = gdf_hex.crs

# Save the data
gdf_hex_counts.to_file(outpath)

'''
#################################
#
#
# Load pre-prepared hex bin data
#
#
################################

gdf_hex_counts = gpd.read_file(outpath)
gdf_hex_counts.rename(columns = {"addPedTick":"addPedTicks", "vehiclePri":"vehiclePriorityCostRatio",
                                    "addVehicle":"addVehicleTicks", "cellCostUp":"cellCostUpdate"}, inplace = True)
print(gdf_hex_counts.columns)


################################
#
#
# Create figure
#
#
################################

#gdf_hex_counts = gdf_hex_counts.to_crs(project_crs)


import matplotlib.pyplot as plt
import contextily as cx

def batch_run_map(df_data, data_col, run_col, rename_dict, title, output_path):

    groupby_columns = ['addVehicleTicks','vehiclePriorityCostRatio','cellCostUpdate']
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

        # Get axis
        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        bb = data.total_bounds
        im, bounds = cx.bounds2img(*bb, ll = True, source=cx.providers.CartoDB.Positron, zoom = 19)

        ax.set_axis_off()
        ax.imshow(im, extent = bounds)
        ax = data.to_crs(epsg=3857).plot(column = data_col, alpha = 0.7, cmap = plt.cm.viridis, edgecolor = 'none', ax = ax, legend = False)

    # Add in row group info to the last axis in the row
    for i in range(p):
        ki = key_indices[i, q*r - 1]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}:\n{}".format(rename_dict[groupby_columns[0]], group_key[0])
        plt.text(1.1,0.5, s, fontsize = 12, transform = ax.transAxes)

    for j in range(q):
        ki = key_indices[0, j*r]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}: {}".format(rename_dict[groupby_columns[1]], group_key[1])
        plt.text(1.05,1.1, s, fontsize = 12, transform = ax.transAxes)

    for k in range(q*r):
        ki = key_indices[1, k]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}: {}".format(rename_dict[groupby_columns[2]], group_key[2])
        plt.text(0.45, -0.1, s, fontsize = 12, transform = ax.transAxes)


    f.suptitle(title, fontsize=16, y = 1)
    f.show()
    plt.savefig(output_path)

map_output_path = "..\\output\\img\\binned_trajectories_w_background.png"
rename_dict = {'addVehicleTicks':"Vehicle\nAddition\nFrequency:",'vehiclePriorityCostRatio':r"$\mathrm{V}$",'cellCostUpdate':r"$\mathrm{\beta}$"}
batch_run_map(gdf_hex_counts, 'loc_count', 'run', rename_dict, "Paths Heatmap", map_output_path)

'''
gdf_hex_run1 = gdf_hex_counts.loc[ ~gdf_hex_counts['loc_count_'].isnull()]
bb = gdf_hex_run1.total_bounds


# In[90]:

im, bounds = cx.bounds2img(*bb, ll = True, source=cx.providers.CartoDB.Positron, zoom = 19)



# In[91]:


f, ax = plt.subplots(1, figsize = (15,15))
ax.set_axis_off()
ax.imshow(im, extent = bounds)
ax = gdf_hex_run1.to_crs(epsg=3857).plot(column='loc_count_',alpha=1, k=7, cmap=plt.cm.viridis, edgecolor='w', linewidth=0.1, ax = ax)


# In[92]:


# Make a map for each run

f, axs = plt.subplots(2,4, figsize = (20,10), sharex= True, sharey = True)
cols = [c for c in gdf_hex_counts.columns if 'loc_coun' in c]
for i in range(len(cols)):
    r = int(i/4)
    c = i % 4
    ax = axs[r,c]
    ax.set_axis_off()
    gdf_hex_run = gdf_hex_counts.loc[ ~gdf_hex_counts[cols[i]].isnull()]
    ax.imshow(im, extent = source)
    ax = gdf_hex_run.to_crs(epsg=3857).plot(column = cols[i], alpha = 0.7, cmap = plt.cm.viridis, edgecolor = 'none', ax = ax, legend = False)


# In[93]:


f.savefig("..\\output\\img\\binned_trajectories_w_background.png")

'''