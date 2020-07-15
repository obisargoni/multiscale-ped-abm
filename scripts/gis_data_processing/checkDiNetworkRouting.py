
# Load OD nodes
# Load itn nodes
# Spatial join so that I have the IDs of the OD nodes
# Load RLnodes and RRI
# Create network using these dataframes
# Query network to identify possible routes between ODs

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import networkx as nx


#############################
#
#
# Initialise paths to inputs and outputs
#
#
#############################
gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\processed_gis_data"

OD_shapefile = "OD_vehicle_nodes_intersect_within.shp"

route_info_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\itn_route_info"
itn_network_edge_data_file = os.path.join(route_info_dir, "itn_edge_list.csv")

output_flows_path = os.path.join(gis_data_dir, "vehicleODFlows.csv")


############################
#
#
# Load the data
#
#
###########################

# Load the vehicle OD nodes and the road link data
gdfOD = gpd.read_file(os.path.join(gis_data_dir, OD_shapefile))

dfEdgeList = pd.read_csv(itn_network_edge_data_file)
graphRL = nx.from_pandas_edgelist(dfEdgeList, source='start_node', target = 'end_node', edge_attr = ['RoadLinkFID','weight'], create_using=nx.DiGraph())
assert graphRL.is_directed()

print(gdfOD.head())

# Now iterate between all OD pairs and find which paths are possible and which are not
dfPossibleFlows = pd.DataFrame(columns = ['O','D', 'flowPossible'])
excludeFIDs = ['osgb4000000029970447'] # Nodes that shouldn't be considered as ODs
ODfids = gdfOD['fid']
for o in ODfids:
	for d in ODfids:
		if (o == d):
			continue
		elif (o in [excludeFIDs]) | (d in excludeFIDs):
			# Record which can't be
			dfPossibleFlows = dfPossibleFlows.append({'O':o,'D':d,'flowPossible':0}, ignore_index=True)
		else:
			try:
				path = nx.shortest_path(graphRL,o,d)
				# record which OD can be traversed
				dfPossibleFlows = dfPossibleFlows.append({'O':o,'D':d,'flowPossible':1}, ignore_index=True)
			except Exception as err:
				# Record which can't be
				dfPossibleFlows = dfPossibleFlows.append({'O':o,'D':d,'flowPossible':0}, ignore_index=True)

# Use this dataframe to create flows
dfFlows = dfPossibleFlows.rename(columns = {'flowPossible':'flow'})
dfFlows['flow'] = dfFlows['flow'] * np.random.rand(dfFlows.shape[0])
dfFlows = dfFlows.set_index(['O','D']).unstack().fillna(0)

# Get rid of multiindex
dfFlows.columns = [c[1] for c in dfFlows.columns]

# Save this dataframe and use as the flows matrix
dfFlows.to_csv(output_flows_path, index=False, header=False)