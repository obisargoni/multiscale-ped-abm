#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
from datetime import datetime as dt
import json
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon

import batch_data_utils as bd_utils

import matplotlib.pyplot as plt
import contextily as cx


###############################
#
#
# Read in pedestrian locations data
#
#
###############################
with open(".//config.json") as f:
    config = json.load(f)

file_datetime_string = config['file_datetime_string']
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\Clapham Common\\"
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

project_crs = {'init': 'epsg:27700'}
wsg_crs = {'init':'epsg:4326'}

hex_polys_file = os.path.join(gis_data_dir, "simple_pedestrian_trips",  "hexgrid1m.shp")

file_re = bd_utils.get_file_regex("pedestrian_locations", file_datetime = file_datetime)
ped_locations_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_locations", file_datetime = file_datetime, suffix = 'batch_param_map')
batch_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

# output paths
output_trajectories_path = os.path.join(data_dir, "hex_bin_trajectories.{}.gpkg".format(file_datetime_string)) 
map_output_path = os.path.join(img_dir, "binned_trajectories_w_background.{}.png".format(file_datetime_string))


#################################
#
#
# Process data
#
#
#################################

df_ped_loc = pd.read_csv(ped_locations_file)

# Now split df into one for locations, one for route coords
lo_cols = [i for i in df_ped_loc.columns if 'Loc' in i]

df_loc = df_ped_loc.dropna(subset = lo_cols)
for c in lo_cols:
    df_loc[c] = df_loc[c].map(float)


gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.LocXString, df_loc.LocYString))
gdf_loc.crs = project_crs

gdf_loc = gdf_loc.to_crs(wsg_crs)


# Read in batch runs metadata
df_run = pd.read_csv(os.path.join(data_dir, batch_file))


selection_columns = ['lambda', 'alpha', 'avNVehicles']
selction_values = [ [0.4,1.6],
                    [0.1,0.9],
                    [15,1]
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
# Use hexagonal tiles to bin pedestrian locations, create heat map of trajectories
gdf_hex = gpd.read_file(hex_polys_file)
print(gdf_hex.crs)

gdf_hex = gdf_hex.to_crs(wsg_crs)
gdf_hex['hex_id'] = np.arange(gdf_hex.shape[0])
gdf_hex = gdf_hex.reindex(columns = ['hex_id','geometry'])

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
gdf_hex_counts.to_file(output_trajectories_path, drive='GPKG')
#################################
#
#
# Load pre-prepared hex bin data
#
#
################################

gdf_hex_counts = gpd.read_file(output_trajectories_path)
gdf_hex_counts.rename(columns = {'avNVehicle':'avNVehicles'}, inplace=True)


################################
#
#
# Create figure
#
#
################################

#gdf_hex_counts = gdf_hex_counts.to_crs(project_crs)


def batch_run_map(df_data, run_selection_dict, data_col, run_col, rename_dict, title, output_path):

    global tbounds

    groupby_columns = ['avNVehicles','alpha','lambda']
    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    p = len(set(run_selection_dict[groupby_columns[0]]))
    q = len(set(run_selection_dict[groupby_columns[1]]))
    r = len(set(run_selection_dict[groupby_columns[2]]))
    key_indices = np.reshape(np.arange(len(keys)), (p,q*r))

    fig_indices = np.reshape(key_indices, (p,q*r))
    f,axs = plt.subplots(p, q*r,figsize=(20,10), sharey=True, sharex = True)

    # Set bounds for map
    map_bounds = [-0.13535773, 51.46418611, -0.13410014, 51.46472197]

    for ki in range(len(keys)):
        group_key = keys[ki]
        data = grouped.get_group(group_key)

        # Check have got a single runs data
        assert data[run_col].unique().shape[0] == 1

        # Get axis
        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        im, bounds = cx.bounds2img(*map_bounds, ll = True, source=cx.providers.Thunderforest.Transport(apikey="ebf9c5aef5b546ab9ea40180032937b5"))
        tbounds = bounds

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

        s = "{}".format(rename_dict[group_key[0]])
        plt.text(1.1,0.5, s, fontsize = 11, transform = ax.transAxes)

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
        if group_key[1] == 0.2:
            t = "Sensitive to traffic"
        elif group_key[1] == 0.8:
            t = "Sensitive to journey time"
        plt.text(0.95,1.2, t, fontsize = 15, transform = ax.transAxes)

    for k in range(q*r):
        ki = key_indices[1, k]
        group_key = keys[ki]

        i,j= np.where(fig_indices == ki)
        assert len(i) == len(j) == 1
        ax = axs[i[0], j[0]]

        s = "{}: {}".format(rename_dict[groupby_columns[2]], group_key[2])
        plt.text(0.45, -0.1, s, fontsize = 11, transform = ax.transAxes)

        # Add some explanitory text
        t = ""
        if group_key[2] == 0.4:
            t = "Plans ahead"
        elif group_key[2] == 1.6:
            t = "Considers nearby\nalternatives more"
        plt.text(0.35,-0.35, t, fontsize = 15, transform = ax.transAxes)

    if title is not None:
        f.suptitle(title, fontsize=16, y = 0.93)
    f.show()
    plt.savefig(output_path)

def batch_run_map_single(df_data, data_col, run_col, rename_dict, title, output_path):

    global tbounds

    groupby_columns = ['avNVehicles']
    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    key_indices = np.arange(len(keys))
    f,axs = plt.subplots(1, len(key_indices),figsize=(20,10), sharey=True, sharex = True)

    # Set bounds for map
    map_bounds = [-0.13525118, 51.46425201, -0.13335044, 51.46488333]

    for ki in range(len(keys)):
        group_key = keys[ki]
        data = grouped.get_group(group_key)

        # Check have got a single runs data
        assert data[run_col].unique().shape[0] == 1

        # Get axis
        ax = axs[ki]

        im, bounds = cx.bounds2img(*map_bounds, ll = True, source=cx.providers.Thunderforest.Transport(apikey="ebf9c5aef5b546ab9ea40180032937b5"))
        tbounds = bounds

        ax.set_axis_off()
        ax.imshow(im, extent = bounds)
        ax = data.to_crs(epsg=3857).plot(column = data_col, alpha = 0.7, cmap = plt.cm.viridis, edgecolor = 'none', ax = ax, legend = False)

        ax.set_title(rename_dict[group_key])

    if title is not None:
        f.suptitle(title, fontsize=16, y = 1)
    f.show()
    plt.savefig(output_path)

rename_dict = {15:"High Vehicle Flow", 1:"Low Vehicle Flow", 'alpha':r"$\mathrm{\alpha}$",'lambda':r"$\mathrm{\lambda}$"}

batch_run_map(gdf_hex_counts, run_selection_dict, 'loc_count', 'run', rename_dict, "Beyond Configuration", map_output_path)
#batch_run_map_single(gdf_hex_counts, 'loc_count', 'run', rename_dict, None, map_output_path)



#####################################
#
#
# Animation
#
#
#####################################

from matplotlib.animation import FuncAnimation, FFMpegWriter
plt.rcParams['animation.ffmpeg_path'] = "C:\\Anaconda3\\bin\\ffmpeg"


# Select batch runs for animation
# Read in batch runs metadata
df_run = pd.read_csv(os.path.join(data_dir, batch_file))


selection_columns = ['lambda', 'alpha', 'avNVehicles']
selction_values = [ [1.0,1.0],
                    [0.5,0.5],
                    [15,1]
                    ]

run_selection_dict = {selection_columns[i]:selction_values[i] for i in range(len(selection_columns))}

for col in selection_columns:
    df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]


# Get crossing locations and origin and destination
gdfCA = gpd.read_file(os.path.join(gis_data_dir,  "CrossingAlternatives.shp"))
gdfOD = gpd.read_file(os.path.join(gis_data_dir, "OD_pedestrian_nodes.shp"))

gdfCA = gdfCA.to_crs(epsg=3857)
gdfOD = gdfOD.to_crs(epsg=3857)

# Get x and y coords to plot at each tick
gdf_loc = gdf_loc.to_crs(epsg=3857)

gdf_loc_display = gdf_loc.loc[ gdf_loc['run'].isin(df_run['run'])]

high_flow_run = df_run.loc[ df_run['avNVehicles']==df_run['avNVehicles'].max(), 'run'].values[0]
low_flow_run = df_run.loc[ df_run['avNVehicles']==df_run['avNVehicles'].min(), 'run'].values[0]

gdf_loc_display['x'] = gdf_loc_display['geometry'].map(lambda g: g.coords[0][0])
gdf_loc_display['y'] = gdf_loc_display['geometry'].map(lambda g: g.coords[0][1])
gdf_loc_display['c'] = gdf_loc_display['run'].replace({high_flow_run:0, low_flow_run:1})

points = gdf_loc_display.groupby('tick').apply(lambda df: df.loc[:, ['x', 'y', 'c']].values).values

# Fine starting tick
starting_tick = int( max( gdf_loc_display.loc[gdf_loc_display['run']==high_flow_run, 'tick'].min(), gdf_loc_display.loc[gdf_loc_display['run']==low_flow_run, 'tick'].min()) )

points = points[starting_tick:]

# filter points to speed up animation
points_filter = []
for i, p in enumerate(points):
    if i%3==0:
        points_filter.append(p)

# Initialise figure
fig, ax = plt.subplots(figsize=(10,10))
xdata, ydata = [], []
scat = ax.scatter([], [], s=20, vmin=0, vmax=1, cmap=plt.cm.bwr, alpha=0.7)
scat.cmap = plt.cm.bwr

def init():
    #map_bounds = [-0.13525118, 51.46425201, -0.13335044, 51.46488333]
    map_bounds = [-0.1351, 51.4643, -0.13425, 51.46465]
    im, bounds = cx.bounds2img(*map_bounds, ll = True, source=cx.providers.Thunderforest.Transport(apikey="ebf9c5aef5b546ab9ea40180032937b5"))
    ax.set_axis_off()
    ax.imshow(im, extent = bounds)

    gdfCA.plot(ax=ax, color = 'green', linewidth = 1.5)
    gdfOD.plot(ax=ax, color='black', linewidth=2)
    return scat,

def update(frame_points):
    # Set x and y data...
    scat.set_offsets(frame_points[:, :2])
    
    # Set sizes...
    #sizes = np.array([200]*len(frame_points))
    #scat.set_sizes(sizes)
    
    # Set colors..
    scat.set_array(frame_points[:, 2])

    # We need to return the updated artist for FuncAnimation to draw..
    # Note that it expects a sequence of artists, thus the trailing comma.
    return scat,

ani = FuncAnimation(fig, update, frames = points_filter[:1000], init_func=init, blit=True)

FFwriter = FFMpegWriter()
ani.save(os.path.join(img_dir, 'lower_level_crossing_animation.{}.mp4'.format(file_datetime_string)), writer = FFwriter)
