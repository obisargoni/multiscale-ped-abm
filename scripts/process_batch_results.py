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
        found = l_re.search(c)
        if found is None:
            return

        x,y = found.groups()
        p = Point(float(x), float(y))
        points.append(p)
    l = LineString(points)
    return l

def exctract_coords(string, i, regex):
    found = regex.search(string)
    if found is None:
        return
    else:
        try:
            coord = float(found.groups()[i])
        except TypeError as e:
            return None
        return coord

#####################################
#
# Globals
#
#####################################

data_dir = "..\\output\\batch\\model_run_data\\"
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
#ped_locations_file = "manual_pedestrian_locations_batch4_6.csv"

# read the data
'''
f = open(os.path.join(data_dir, ped_locations_file))
csv_reader = csv.reader(f)
data = []
for row in csv_reader:
    if len(row) != 5:
        continue
    else:
        data.append(row)
        
df_ped_loc = pd.DataFrame(data)
data = None
df_ped_loc.columns = df_ped_loc.iloc[0]
df_ped_loc = df_ped_loc.iloc[1:]
'''
df_ped_loc = pd.read_csv(os.path.join(data_dir, ped_locations_file))

# Now split df into one for locations, one for route coords
lo_cols = [i for i in df_ped_loc.columns if 'Loc' in i]

df_loc = df_ped_loc.dropna(subset = lo_cols)
for c in lo_cols:
    df_loc[c] = df_loc[c].map(float)

##################################
#
# Get pedestrian trajectories as linestrings
#
##################################

# Create geopandas data frame in order to crete coordiante object from coords
gdf_traj = gdf_loc.groupby(['run', 'ID'])['geometry'].apply(lambda x: LineString(x.tolist())).reset_index()
gdf_traj['traj_length'] = gdf_traj['geometry'].map(lambda l: l.length)
gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.LocXString, df_loc.LocYString))


#################################
#
# Get pedestrian primary route coordinates and direct route linestring
#
#################################

# Work flow to get primary route coordinates
file_prefix = 'pedestrian_initial_route'
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})\.csv")
ped_primary_route_file = most_recent_directory_file(data_dir, file_re)
#ped_primary_route_file = "manual_pedestrian_primary_route_batch4_6.csv"

# read the data
'''
f = open(os.path.join(data_dir, ped_primary_route_file))
csv_reader = csv.reader(f)
data = []
for row in csv_reader:
    if len(row) != 4:
        continue
    else:
        data.append(row)
        
df_prim = pd.DataFrame(data)
data = None
df_prim.columns = df_prim.iloc[0]
df_prim = df_prim.iloc[1:]
'''
df_prim = pd.read_csv(os.path.join(data_dir, ped_primary_route_file))


# Primary route coords are removed as agent progresses, so select only the first time step for each agent
df_prim.sort_values(by = ['run','ID','tick'], inplace=True)
df_prim_first = df_prim.groupby(['run','ID']).apply(lambda df: df.iloc[0])
df_prim_first.index = np.arange(df_prim_first.shape[0])
df_prim_first['geometry'] = df_prim_first['InitialRouteCoordinatesString'].map(lambda x: linestring_from_coord_string(x))
df_prim_first.dropna(subset = ['geometry'], inplace = True)

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
gdf_comb = gdf_comb.reindex(columns = ['run', 'ID', 'traj_length', 'length', 'deviation', 'pct_deviation'])
gdf_comb = gdf_comb.dropna(subset = ['pct_deviation'])


###############################
#
# Load pedestrian ODs and Road Polygons data and use to calculate proportion of pedestrian journeys that use/don't use crossing
#
###############################
gdf_ped_od = gpd.read_file(os.path.join(geo_data_dir, ped_od_file))
gdf_veh_poly = gpd.read_file(os.path.join(geo_data_dir, vehicle_polygons))

# Buffer ODs slightly
gdf_ped_od.rename(columns = {"id":"od_id"}, inplace=True)
gdf_ped_od['geometry'] = gdf_ped_od['geometry'].buffer(0.5)

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
df_crossing_ped = df_crossing_ped.loc[  ((df_crossing_ped['od_id_start']==2) & (df_crossing_ped['od_id_end']==1)) |
                                        ((df_crossing_ped['od_id_start']==4) & (df_crossing_ped['od_id_end']==1)) |
                                        ((df_crossing_ped['od_id_start']==3) & (df_crossing_ped['od_id_end']==2)) |

                                        ((df_crossing_ped['od_id_start']==1) & (df_crossing_ped['od_id_end']==2)) |
                                        ((df_crossing_ped['od_id_start']==1) & (df_crossing_ped['od_id_end']==4)) |
                                        ((df_crossing_ped['od_id_start']==2) & (df_crossing_ped['od_id_end']==3)) ]

# Filter routes by those that could pass through crossing
gdf_cross_traj = pd.merge(gdf_traj, df_crossing_ped, on = ['run','ID'], how = 'inner')
gdf_cross_traj = gpd.GeoDataFrame(gdf_cross_traj)

# Get crossing polygons
gdf_crossing_polys = gdf_veh_poly.loc[gdf_veh_poly['priority'] == 'pedestrian_crossing'].reindex(columns = ['fid','priority','geometry'])

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
# Filter data by desired batch runs
#
###############################

# Read in batch runs metadata
file_prefix = 'pedestrian_initial_route'
file_suffix = 'batch_param_map'
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})\."+file_suffix + "\.csv")
batch_file = most_recent_directory_file(data_dir, file_re)
df_run = pd.read_csv(os.path.join(data_dir, batch_file))

selection_columns = ['vehiclePriorityCostRatio', 'addVehicleTicks', 'cellCostUpdate']
selction_values = [ [1,10],
                    [5,10],
                    [0,10]
                    ]

run_selection_dict = {selection_columns[i]:selction_values[i] for i in range(len(selection_columns))}

for col in selection_columns:
    df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]

# Join this to the data
gdf_comb = pd.merge(gdf_comb, df_run, left_on = 'run', right_on = 'run', how = 'inner')
df_cross_pct = pd.merge(df_cross_pct, df_run, left_on = 'run', right_on = 'run', how = 'inner')

###############################
#
# Make figures
#
###############################
import matplotlib.pyplot as plt
import seaborn as sn


# Histogram of deviation from straight line distance
groupby_columns = ['addVehicleTicks','vehiclePriorityCostRatio','cellCostUpdate']
grouped = gdf_comb.groupby(groupby_columns)
keys = list(grouped.groups.keys())

x = len(run_selection_dict[groupby_columns[0]])
y = len(run_selection_dict[groupby_columns[1]])
z = len(run_selection_dict[groupby_columns[2]])
key_indices = np.reshape(np.arange(len(keys)), (x,y,z))

fig_indices = np.reshape(key_indices, (x, y*z))

f,axs = plt.subplots(x, y*z,figsize=(10,20), sharey=True, sharex = True)
for ki in range(len(keys)):
    group_key = keys[ki]
    data = grouped.get_group(group_key)
    assert data['run'].unique().shape[0] == 1

    stats = data['pct_deviation'].describe().map(lambda x: round(x,4))
    results_string = "mean = {}\nstd = {}".format(stats['mean'], stats['std'])

    i,j= np.where(fig_indices == ki)
    assert len(i) == len(j) == 1
    ax = axs[i[0], j[0]]

    ax.hist(data['pct_deviation'])
    plt.text(0.1,0.9, results_string, fontsize = 9, transform = ax.transAxes)
    #ax.set_title("{},{},{}".format(*group_key), fontsize = 9)

# Add in row group legend values
for i in range(x):
    ki = key_indices[i, y-1, z-1]
    group_key = keys[ki]

    i,j= np.where(fig_indices == ki)
    assert len(i) == len(j) == 1
    ax = axs[i[0], j[0]]

    s = "{}:\n{}".format(groupby_columns[0], group_key[0])
    plt.text(1.1,0.5, s, fontsize = 9, transform = ax.transAxes)

for j in range(y):
    ki = key_indices[0, j, 0]
    group_key = keys[ki]

    i,j= np.where(fig_indices == ki)
    assert len(i) == len(j) == 1
    ax = axs[i[0], j[0]]

    s = "{}: {}".format(groupby_columns[1], group_key[1])
    plt.text(0.75,1.1, s, fontsize = 9, transform = ax.transAxes)

for ki in np.nditer(key_indices):
    group_key = keys[ki]

    i,j= np.where(fig_indices == ki)
    assert len(i) == len(j) == 1
    ax = axs[i[0], j[0]]

    ax.set_title("{}:\n{}".format(groupby_columns[2], group_key[2]), fontsize = 9)

f.show()
plt.savefig(os.path.join(img_dir, path_deviation_fig))

df_cross_pct = df_cross_pct.dropna(subset=['cross_pct'])
f,ax = plt.subplots(1,1,figsize=(10,20), sharey=True, sharex = True)
ax.bar(df_cross_pct['run'], df_cross_pct['cross_pct'])
f.show()
plt.savefig(os.path.join(img_dir, crossing_percent_fig))