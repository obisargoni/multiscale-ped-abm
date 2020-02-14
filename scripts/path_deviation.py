from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from shapely.geometry import LineString
from shapely.geometry import Point

file_list = os.listdir(data_dir)
file_list[0]

data_dir = "..\\data\\model_run_exports\\"
file_prefix = "pedestrian_locations"

file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})")
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

def dt_from_file_name(file_name, regex):
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt
    
def dt_from_file_name(file_name, regex = file_re):
    s = regex.search(file_name)
    file_dt = dt.strptime(s.group(), "%Y%b%d%H%M%S")
    return file_dt

def dt_from_file_name(file_name):
    regex = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})")
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt
    

filtered_files = [f for f in files if file_re.search(f) is not None]
filtered_files.sort(key = dt_from_file_name, reverse=True)
ped_locations_file = filtered_files[0]

# read the data
df_ped_loc = pd.read_csv(os.path.join(data_dir, ped_locations_file))


loc_cols = ['Location', 'RouteCoordinate']
for c in loc_cols:
    for i in [0,1]:
        df_ped_loc[c+str(i)] = df_ped_loc[c].map(lambda x: l_re.search(x).groups()[i])
        
df_ped_loc.head()

# Now split df into one for locations, one for route coords
route_cols = [i for i in df_ped_loc.columns if 'RouteCoordinate' in i]
lo_cols = [i for i in df_ped_loc.columns if 'Location' in i]

df_loc = df_ped_loc.drop(route_cols, axis=1)
df_route = df_ped_loc.drop(lo_cols, axis=1)

# Create geopandas data frame in order to crete coordiante object from coords
gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.Location0, df_loc.Location1))
for c in loc_cols:
    for i in [0,1]:
        df_ped_loc[c+str(i)] = df_ped_loc[c].map(lambda x: float(l_re.search(x).groups()[i]))
        
df_loc = df_ped_loc.drop(route_cols, axis=1)
df_route = df_ped_loc.drop(lo_cols, axis=1)

gdf_loc = gpd.GeoDataFrame(df_loc, geometry=gpd.points_from_xy(df_loc.Location0, df_loc.Location1))
gdf_traj = gdf_loc.groupby(['run', 'ID'])['geometry'].apply(lambda x: LineString(x.tolist())).reset_index()
gdf_traj['traj_length'] = gdf_traj['geometry'].map(lambda l: l.length)


# Work flow to get primary route coordinates
prefix = 'pedestrian_primary_route'
file_re = re.compile(file_prefix+r"\.(\d{4})\.([a-zA-Z]{3})\.(\d{2})\.(\d{2})_(\d{2})_(\d{2})")

filtered_files = [f for f in files if file_re.search(f) is not None]
filtered_files.sort(key = dt_from_file_name, reverse=True)
ped_locations_file = filtered_files[0]


df_prim = pd.read_csv(os.path.join(data_dir, ped_route_file))

def linestring_from_coord_string(coord_string):
    points = []
    for c in coord_string.split("),("):
        x,y = l_re.search(c).groups()
        p = Point(float(x), float(y))
        points.append(p)
    l = LineString(points)
    return l


# Primary route coords are removed as agent progresses
# So Filter out all time steps other than 1
df_prim = df_prim.loc[ df_prim['tick']==1]
df_prim['primary_route_line'] = df_prim['PrimaryRouteCoordiantesString'].map(lambda x: linestring_from_coord_string(x))
gdf_prim = gpd.GeoDataFrame(df_prim)
gdf_prim.rename(columns = {'primary_route_line':'geometry'}, inplace=True)
gdf_prim['length'] = gdf_prim['geometry'].map(lambda l: l.length)


# Join geodataframes together and calculate % difference from primary route path
gdf_comb = pd.merge(gdf_traj, gdf_prim, on = ['run', 'ID'], how='outer')
gdf_comb['deviation'] = gdf_comb['traj_length'] - gdf_comb['length']
gdf_comb['pct_deviation'] = gdf_comb['deviation'] / gdf_comb['length']
