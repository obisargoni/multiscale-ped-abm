# Check pedestrian perception of vehicle dominance/vehicle priority cost
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import random

gis_data_dir = "..\\data\\single_road_data\\"
road_link_path = os.path.join(gis_data_dir, "mastermap-itn RoadLink Intersect Within with orientation.shp")

gdf_road_link = gpd.read_file(road_link_path)
gdf_road_link['length'] = gdf_road_link['geometry'].length

df_ped_est = pd.DataFrame({"length":gdf_road_link['length']})
df_vehicle_numbers = pd.DataFrame({"vehicle_number":np.arange(10)})

df_ped_est['key'] = 0
df_vehicle_numbers['key'] = 0
df_ped_est = pd.merge(df_ped_est, df_vehicle_numbers)
df_ped_est.drop('key', axis = 1, inplace=True)

# Now define constants
vehicle_width = 1.8
vehicle_length = 4.0
vehicle_mph = 10
vehicle_speed = vehicle_mph * 1609.34 * (1.0/(60*60)) # 20mph in ms-1
vehicle_reaction_time = 0.85 # Taken from trl report https://trl.co.uk/sites/default/files/PPR313_new.pdf
pedestrian_speed = 0.8
lane_width = 3.65


def vehicle_space(vn, rl, lw = lane_width, vw = vehicle_width, vl = vehicle_length):
    vehicle_density = vn / (rl * lw)

    return vehicle_density * vw * vl

def vehicle_separation_space(vn, rl, lw = lane_width, vw = vehicle_width, vs = vehicle_speed, vrt = vehicle_reaction_time):
    vehicle_density = vn / (rl * lw)

    return vehicle_density * vw * vs * vrt

def pedestrian_gap_space(vn, rl, lw = lane_width, vw = vehicle_width, vs = vehicle_speed, ps = pedestrian_speed):
    vehicle_density = vn / (rl * lw)

    # Accounting for space ahead of vehicle required for pedestrian to cross only non zero if there are vehicles present
    if vehicle_density == 0:
        pedestrian_gap_space = 0
    else:
        pedestrian_gap_space = vw * vs * (lw/ps) * (1/rl)
    return pedestrian_gap_space


def vehicle_road_proportion(vn, rl, lw = lane_width, vw = vehicle_width, vl = vehicle_length, vs = vehicle_speed, vrt = vehicle_reaction_time, ps = pedestrian_speed):
    v_space = vehicle_space(vn, rl, lw, vw, vl)

    v_sep_space = vehicle_separation_space(vn, rl, lw, vw, vs, vrt)

    ped_gap_space = pedestrian_gap_space(vn, rl, lw, vw, vs, ps)

    vrp = v_space + v_sep_space + ped_gap_space

    # Calculation can go over 1, in practive vehicles should not be able to occupy more than 100% of road space
    return min(vrp, 1.0)


df_ped_est['vehicle_space'] = df_ped_est.apply(lambda row: vehicle_space(row['vehicle_number'], row['length']), axis = 1)
df_ped_est['vehicle_separation_space'] = df_ped_est.apply(lambda row: vehicle_separation_space(row['vehicle_number'], row['length']), axis = 1)
df_ped_est['pedestrian_gap_space'] = df_ped_est.apply(lambda row: pedestrian_gap_space(row['vehicle_number'], row['length']), axis = 1)
df_ped_est["vehicle_road_proportion"] = df_ped_est.apply(lambda row: vehicle_road_proportion(row['vehicle_number'], row['length']), axis = 1)