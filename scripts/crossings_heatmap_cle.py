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
import networkx as nx
from shapely.geometry import Point, Polygon

import batch_data_utils as bd_utils

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

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\Clapham Common\\processed_gis_data"
data_dir = config['batch_data_dir']
img_dir = "..\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

project_crs = {'init': 'epsg:27700'}
wsg_crs = {'init':'epsg:4326'}

nbins = None
bin_dist = 2

hex_polys_file = os.path.join("S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\Clapham Common\\simple_pedestrian_trips", "hexgrid1m.shp")

data_paths = bd_utils.get_data_paths(file_datetime_string, data_dir)
cross_events_file = data_paths["cross_events"]
ped_routes_file = data_paths['pedestrian_routes']
batch_file = data_paths["batch_file"]

output_paths = bd_utils.get_ouput_paths(file_datetime_string, data_dir, nbins = bin_dist, ttc_threshold=1.5)
output_link_cross_entropy = output_paths["output_link_cross_entropy"]


pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])


# output paths
output_trajectories_path = os.path.join(data_dir, "hex_bin_crossings.{}.gpkg".format(file_datetime_string)) 
map_output_path = os.path.join(img_dir, "binned_crossings_w_background.{}.png".format(file_datetime_string))

gdfPaveLinks = gpd.read_file(pavement_links_file)
gdfPaveNodes = gpd.read_file(pavement_nodes_file)
gdfORLinks = gpd.read_file(or_links_file)

# Get networkx graph and node positions
pavement_graph = nx.Graph()
gdfPaveLinks['length'] = gdfPaveLinks['geometry'].length.map(lambda x: int(x)).replace(0,1) # Repalce 0 with 1 to prevent links with 0 weight. Only affects 4 links. (gdfPaveLinks['geometry'].length < 1).value_counts()
gdfPaveLinks['edge_data'] = gdfPaveLinks.apply(lambda row: {'length':row['length'], 'fid':row['fid']}, axis=1)
edges = gdfPaveLinks.loc[:, ['MNodeFID', 'PNodeFID', 'edge_data']].values
pavement_graph.add_edges_from(edges)


#################################
#
#
# Get road crosisng data
#
#
#################################
dfRun = pd.read_csv(os.path.join(data_dir, batch_file))


weight_params = range(0,100,101)
dfPedRoutes, dfPedRoutes_removedpeds = bd_utils.load_and_clean_ped_routes(gdfPaveLinks, gdfORLinks, gdfPaveNodes, pavement_graph, weight_params, ped_routes_path = ped_routes_file, strategic_path_filter=False)
dfCrossEvents = bd_utils.load_and_clean_cross_events(gdfPaveLinks, cross_events_path = cross_events_file)
dfCrossEvents['cross_linestring'] = dfCrossEvents['CrossingCoordinatesString'].map(lambda x: bd_utils.linestring_from_crossing_coord_string(x))

dfCLE = bd_utils.calculate_average_link_level_crossing_location_entropy(dfCrossEvents, dfPedRoutes, gdfPaveLinks, gdfPaveNodes, gdfORLinks, dfRun, nbins = None, bin_dist = 2, output_path = output_link_cross_entropy)

gdfCrossLines = gpd.GeoDataFrame(dfCrossEvents.reindex(columns = ['run','ID','cross_linestring']), geometry = 'cross_linestring', crs = project_crs)
gdfCrossLines = gdfCrossLines.to_crs(wsg_crs)

#################################
#
#
# Load the hexagonal bins and join with crossing links to highlight crossing locations 
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
gdf_hex_traj = gpd.sjoin(gdf_hex, gdfCrossLines, op = 'intersects')

# Get point count per hex per run
ser_hex_counts = gdf_hex_traj.groupby(["hex_id", "run"])['index_right'].apply(lambda df: df.shape[0])
df_hex_counts = pd.DataFrame(ser_hex_counts).reset_index().rename(columns = {'index_right':'loc_count'})
df_hex_counts['loc_prop'] = df_hex_counts.groupby('run')['loc_count'].transform(lambda s: s / 200) # divide by number of ped agents

# Join back with the polygons and with df run to get parameter info

df_hex_counts = pd.merge(df_hex_counts, gdf_hex, on = 'hex_id', how = 'left')
df_hex_counts = pd.merge(df_hex_counts, dfRun, on = 'run', how = 'inner') # Only include the runs we are interested in

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


import matplotlib.pyplot as plt
import contextily as cx

def batch_run_cross_map_single(df_data, data_col, run_col, rename_dict, title, output_path, map_bounds =  [-0.1351, 51.4643, -0.13425, 51.46465]):

    global tbounds

    groupby_columns = ['run']
    grouped = df_data.groupby(groupby_columns)
    keys = list(grouped.groups.keys())

    key_indices = np.arange(len(keys))
    f,axs = plt.subplots(1, len(key_indices),figsize=(15,7), sharey=True, sharex = True)

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
        ax = data.to_crs(epsg=3857).plot(column = data_col, alpha = 0.7, cmap = plt.cm.viridis, vmin=0.0, vmax=1.0, edgecolor = 'none', ax = ax, legend = False)

        # Add annotation to show CLE value
        s = r"$CLE=$" + str(np.round(data['mean_link_cross_entropy'].unique()[0], 3))
        plt.text(0.45,0.9, s, fontsize = 11, transform = ax.transAxes)

    # add colour bar
    smap = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=0, vmax=1.0))
    cbar = f.colorbar(smap, ax=axs[-1], fraction=0.1, shrink = 0.5)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label("% of road crossings", x=0, rotation=270)
    cbar.set_ticks([0,0.2,0.4,0.6,0.8,1.0])
    cbar.set_ticklabels([0,20,40,60,80,100])


    if title is not None:
        f.suptitle(title, fontsize=16, y = 0.8)
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0.07)

rename_dict = {15:"High Vehicle Flow", 1:"Low Vehicle Flow", 'alpha':r"$\mathrm{\alpha}$",'lambda':r"$\mathrm{\lambda}$"}

# Select two runs to compare the crossing location entropy
dfCLE['cle_quintile'] = pd.qcut(dfCLE['mean_link_cross_entropy'], 5)
dfCLE['cle_quintile_num'] = dfCLE['cle_quintile'].replace(dict(enumerate(dfCLE['cle_quintile'].value_counts().index)))

# select runs to plots data for 
dfCLEPlot = dfCLE.loc[ dfCLE['cle_quintile_num'].isin([1,3])].groupby('cle_quintile_num').apply(lambda df: df.iloc[10:11])

gdf_hex_counts_plot = pd.merge(gdf_hex_counts, dfCLEPlot.reindex(columns = ['run','mean_link_cross_entropy', 'cle_quintile_num']), on = 'run', how = 'inner')

batch_run_cross_map_single(gdf_hex_counts_plot, 'loc_prop', 'run', rename_dict, 'Crossing Location Entropy', map_output_path, map_bounds =  [-0.1351, 51.4643, -0.13425, 51.46465])




#
# Animation of road crossing - sigspatial replicate
#


from matplotlib.animation import FuncAnimation, FFMpegWriter
plt.rcParams['animation.ffmpeg_path'] = "C:\\Anaconda3\\bin\\ffmpeg"

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\Clapham Common\\sigspatial_replicate\\geodata"


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
gdfCA = gpd.read_file(os.path.join(gis_data_dir,  "beyond_crossing.shp"))
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
fig, ax = plt.subplots(figsize=(7,7))
xdata, ydata = [], []
scat = ax.scatter([], [], s=30, vmin=0, vmax=1, cmap=plt.cm.bwr, alpha=0.7)
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