from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from shapely.geometry import LineString
from shapely.geometry import Point

######################################
#
# Functions
#
######################################

def dt_from_file_name(file_name, regex):
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt

def most_recent_directory_file(directory, file_regex):
	files = os.listdir(directory)
	filtered_files = [f for f in files if file_regex.search(f) is not None]
	filtered_files.sort(key = lambda x: dt_from_file_name(x, file_regex), reverse=True)
	return filtered_files[0]

def linestring_from_coord_string(coord_string):
    points = []
    for c in coord_string.split("),("):
        x,y = l_re.search(c).groups()
        p = Point(float(x), float(y))
        points.append(p)
    l = LineString(points)
    return l

#####################################
#
# Globals
#
#####################################

data_dir = "..\\data\\model_run_exports\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")
output_directory = "..\\output\\processed_data\\"

geo_data_dir = "..\\data\\single_road_data\\"
ped_od_file = "OD_pedestrian_nodes.shp"
vehicle_polygons = "topographicAreaVehicle-withPriority-RoadLinkGrouping-Crossing.shp"


####################################
#
# Get pedestrian locations and route as seperate dataframes
#
####################################
    
file_prefix = "pedestrian_locations"
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})")
ped_locations_file = most_recent_directory_file(data_dir, file_re)

# read the data
df_ped_loc = pd.read_csv(os.path.join(data_dir, ped_locations_file))

loc_cols = ['Location', 'RouteCoordinate']
for c in loc_cols:
    for i in [0,1]:
        df_ped_loc[c+str(i)] = df_ped_loc[c].map(lambda x: float(l_re.search(x).groups()[i]))

# Now split df into one for locations, one for route coords
route_cols = [i for i in df_ped_loc.columns if 'RouteCoordinate' in i]
lo_cols = [i for i in df_ped_loc.columns if 'Location' in i]

df_loc = df_ped_loc.drop(route_cols, axis=1)
df_route = df_ped_loc.drop(lo_cols, axis=1)


##################################
#
# Get pedestrian trajectories as linestrings
#
##################################

# Create geopandas data frame in order to crete coordiante object from coords
gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.Location0, df_loc.Location1))
gdf_traj = gdf_loc.groupby(['run', 'ID'])['geometry'].apply(lambda x: LineString(x.tolist())).reset_index()
gdf_traj['traj_length'] = gdf_traj['geometry'].map(lambda l: l.length)


#################################
#
# Get pedestrian primary route coordinates and direct route linestring
#
#################################

# Work flow to get primary route coordinates
file_prefix = 'pedestrian_primary_route'
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})")
ped_primary_route_file = most_recent_directory_file(data_dir, file_re)

df_prim = pd.read_csv(os.path.join(data_dir, ped_primary_route_file))


# Primary route coords are removed as agent progresses, so select only the first time step for each agent
df_prim.sort_values(by = ['run','ID','tick'], inplace=True)
df_prim_first = df_prim.groupby(['run','ID']).apply(lambda df: df.iloc[0])
df_prim_first.index = np.arange(df_prim_first.shape[0])
df_prim_first['geometry'] = df_prim_first['PrimaryRouteCoordinatesString'].map(lambda x: linestring_from_coord_string(x))

gdf_prim = gpd.GeoDataFrame(df_prim_first)
gdf_prim['length'] = gdf_prim['geometry'].map(lambda l: l.length)


################################
#
# Join dataframes to compare walked route to direct primary route
#
################################

gdf_comb = pd.merge(gdf_traj, gdf_prim, on = ['run', 'ID'], how='outer')
gdf_comb['deviation'] = gdf_comb['traj_length'] - gdf_comb['length']
gdf_comb['pct_deviation'] = gdf_comb['deviation'] / gdf_comb['length']


###############################
#
# Load pedestrian ODs and Road Polygons data and use to calculate proportion of pedestrian journeys that use/don't use crossing
#
###############################
gdf_ped_od = gpd.read_file(os.path.join(geo_data_dir, ped_od_file))
gdf_veh_poly = gpd.read_file(os.path.join(geo_data_dir, vehicle_polygons))

# Match up ped ODs to ped trajectories

###############################
#
# Save processed data
#
###############################
