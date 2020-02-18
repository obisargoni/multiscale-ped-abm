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

data_dir = "..\\data\\model_run_exports"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")
output_directory = "..\\output\\processed_data\\"

geo_data_dir = "..\\data\\single_road_data\\"
ped_od_file = "OD_pedestrian_nodes.shp"
vehicle_polygons = "topographicAreaVehicle-withPriority-RoadLinkGrouping-Crossing.shp"


# Output paths
img_dir = "..\\output\\img\\"
path_deviation_fig = "path_deviation.png"
crossing_percent_fig = "crossing_percent.png"

output_data_dir = "..\\output\\processed_data\\"


####################################
#
# Get pedestrian locations and route as seperate dataframes
#
####################################
    
file_prefix = "pedestrian_locations"
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})\.csv")
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
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})\.csv")
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

# Buffer ODs slightly
gdf_ped_od.rename(columns = {"id":"od_id"}, inplace=True)
gdf_ped_od['geometry'] = gdf_ped_od['geometry'].buffer(0.1)

gdf_prim['start'] = gdf_prim['geometry'].map(lambda l: Point(l.coords[0]))
gdf_prim['end'] = gdf_prim['geometry'].map(lambda l: Point(l.coords[-1]))

# Spatial join
gdf_prim.set_geometry("start")
gdf_prim = gpd.sjoin(gdf_prim, gdf_ped_od, how='left', op='intersects')
gdf_prim.drop('index_right', axis=1, inplace=True)
gdf_prim.rename(columns = {"od_id":"od_id_start"}, inplace=True)

gdf_prim.set_geometry("end")
gdf_prim = gpd.sjoin(gdf_prim, gdf_ped_od, how='left', op='intersects')
gdf_prim.drop('index_right', axis=1, inplace=True)
gdf_prim.rename(columns = {"od_id":"od_id_end"}, inplace=True)

gdf_prim.set_geometry("geometry")

# Select ped routes that could have made use of a crossing
df_crossing_ped = gdf_prim.reindex(columns = ['run','ID', 'od_id_start', 'od_id_end'])
df_crossing_ped = df_crossing_ped.loc[	((df_crossing_ped['od_id_start']==2) & (df_crossing_ped['od_id_end']==1)) |
										((df_crossing_ped['od_id_start']==4) & (df_crossing_ped['od_id_end']==1)) |
										((df_crossing_ped['od_id_start']==3) & (df_crossing_ped['od_id_end']==2)) |

										((df_crossing_ped['od_id_start']==1) & (df_crossing_ped['od_id_end']==2)) |
										((df_crossing_ped['od_id_start']==1) & (df_crossing_ped['od_id_end']==4)) |
										((df_crossing_ped['od_id_start']==2) & (df_crossing_ped['od_id_end']==3)) ]

# Filter routes by those that could pass through crossing
gdf_cross_traj = pd.merge(gdf_traj, df_crossing_ped, on = ['run','ID'], how = 'inner')
gdf_cross_traj = gpd.GeoDataFrame(gdf_cross_traj)

# Get crossing polygons
gdf_crossing_polys = gdf_veh_poly.loc[gdf_veh_poly['priority'] == 'pedestrian'].reindex(columns = ['fid','priority','geometry'])

gdf_cross_traj = gpd.sjoin(gdf_cross_traj, gdf_crossing_polys, how = 'left')
gdf_cross_traj['use_crossing'] = ~gdf_cross_traj['fid'].isnull()
gdf_cross_traj.drop(list(gdf_crossing_polys.columns) + ['index_right'], axis=1, inplace=True)

df_cross_pct = gdf_cross_traj.groupby('run').apply(lambda df: df.loc[df['use_crossing']==True].count() / df.count())
df_cross_pct = df_cross_pct.rename(columns = {'run':'cross_pct'}).reindex(columns = ['cross_pct']).reset_index()

###############################
#
# Save processed data
#
###############################



###############################
#
# Make figures
#
###############################
import matplotlib.pyplot as plt
import seaborn as sn

runs = gdf_comb['run'].unique()
n = len(runs)

# Histogram of deviation from straight line distance
f,axs = plt.subplots(1,n,figsize=(10,20), sharey=True, sharex = True)
for i in range(n):
	data = gdf_comb.loc[gdf_comb['run']==runs[i], 'pct_deviation']
	f.get_axes()[i].hist(data, bins = 20)

f.show()
plt.savefig(os.path.join(img_dir, path_deviation_fig))


f,ax = plt.subplots(1,1,figsize=(10,20), sharey=True, sharex = True)
ax.bar(df_cross_pct['run'], df_cross_pct['cross_pct'])
f.show()
plt.savefig(os.path.join(img_dir, crossing_percent_fig))