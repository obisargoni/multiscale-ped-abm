import pandas as pd
import numpy as np
import os
import geopandas as gpd
import json

from shapely.geometry import Point

##################
#
# Config
#
##################

with open("config.json") as f:
    config = json.load(f)

# Proportion of ITN nodes to use as vehicle ODs
prop_random_ODs = 0.3

gis_data_dir = config['gis_data_dir']
processed_gis_dir = os.path.join(gis_data_dir, "processed_gis_data")

pavement_nodes_file = os.path.join(processed_gis_dir, config["pavement_nodes_file"])

pedestrian_od_flows = os.path.join(processed_gis_dir, config['pedestrian_od_flows'])
pedestrian_od_file = os.path.join(processed_gis_dir, config['pedestrian_od_file'])

poi_file = os.path.join(gis_data_dir, config["poi_file"])
centre_poi_ref = config["centre_poi_ref"]
dist_from_centre_threshold = 50

#################
#
# Load data and select pedestrian ODs
#
#################

# Select Origin pavement nodes based on POIs
'''
gdf_pois = gpd.read_file(poi_file)
gdf_pave_node = gpd.read_file(pavement_nodes_file)

centre_poi_geom = gdf_pois.loc[ gdf_pois['ref_no'] == centre_poi_ref, 'geometry'].values[0]

gdf_pave_node['dist_to_centre'] = gdf_pave_node['geometry'].distance(centre_poi_geom)

Os = gdf_pave_node.sort_values(by = 'dist_to_centre', ascending = True)['fid'].values[0:1]

# Select destination nodes randomly
candidates = gdf_pave_node.loc[ gdf_pave_node['dist_to_centre'] > dist_from_centre_threshold, 'fid'].values

nDs = int(prop_random_ODs * len(candidates))

Ds = np.random.choice(candidates, nDs)

ODs = np.concatenate([Os, Ds])

gdfODs = gdf_pave_node.loc[ gdf_pave_node['fid'].isin(ODs), ['fid', 'geometry']]

gdfODs.to_file(pedestrian_od_file)
'''

def displace_point(p, d, bearing):
	x = p.x + d*np.sin(bearing)
	y = p.y + d*np.cos(bearing)
	return Point([x,y])


# Load ped ODs to get number of origins/destinations
gdfODs = gpd.read_file(pedestrian_od_file)

points = gdfODs['geometry'].values
ds = [0.1]*gdfODs.shape[0]
random_bearings = np.random.rand(gdfODs.shape[0])
new_points = list(map(displace_point, points, ds, random_bearings))

gdfODs['geometry'] = new_points

Os = gdfODs.loc[ gdfODs['fid'] == 'pave_node_9', 'fid'].values
Ds = gdfODs.loc[ gdfODs['fid'] != 'pave_node_9', 'fid'].values
ODs = np.concatenate([Os, Ds])

#################
#
# Generate Flows
#
#################

# Initialise OD flows matrix
flows = np.zeros([len(ODs), len(ODs)])
dfFlows = pd.DataFrame(flows, columns = ODs, index = ODs)

# Set random flows between origins and destinations
for o in Os:
	for d in Ds:
		if (o == d):
			continue

		f_ij = np.random.rand()
		dfFlows.loc[o, d] = f_ij

		f_ji = np.random.rand()
		dfFlows.loc[d, o] = f_ji


dfFlows.to_csv(pedestrian_od_flows, index=False)