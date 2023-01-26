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
        s = r"$CLE_r=$" + str(np.round(data['mean_link_cross_entropy'].unique()[0], 3))
        plt.text(0.45,0.9, s, fontsize = 11, transform = ax.transAxes)

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