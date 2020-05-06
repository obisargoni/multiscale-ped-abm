
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
itn_link_shapefile = "mastermap-itn RoadLink Intersect Within with orientation.shp" # The shapefile of the road network used in the model

route_info_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\itn_route_info"
road_route_info_path = os.path.join(route_info_dir, "extracted_RRI.csv")
road_node_info_path = os.path.join(route_info_dir, "extracted_RLNodes.csv")

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
gdfLink = gpd.read_file(os.path.join(gis_data_dir, itn_link_shapefile))

# Calculate the linestring lengths, to use as network edge weights
gdfLink['weight'] = gdfLink['geometry'].map(lambda g: g.length)

# No need to join together because the OD nodes file has the fids in
dfRLNode = pd.read_csv(road_node_info_path)
dfRRI = pd.read_csv(road_route_info_path)

# Filter the nodes to be just those in the ITN Link data used for the model
dfRLNode = dfRLNode.loc[ dfRLNode["RoadLinkFID"].isin(gdfLink.fid)]

# Check all road links are present and that thy all have 2 nodes
assert len(dfRLNode.RoadLinkFID.unique()) == len(gdfLink.fid.unique())
assert dfRLNode['PlusNodeFID'].isnull().any() == False
assert dfRLNode['MinusNodeFID'].isnull().any() == False

# Merge with data on which links are directed
dfRLNode = pd.merge(dfRLNode, dfRRI, left_on = "RoadLinkFID", right_on = "DirectedLinkFID", how = "left", indicator=True)
print(dfRLNode["_merge"].value_counts())
dfRLNode.drop("_merge", axis = 1, inplace = True)

# Merge with link weight and id data
dfRLNode = pd.merge(dfRLNode, gdfLink, left_on = "RoadLinkFID", right_on = "fid", how = "outer", indicator = True)

# Now build the edge list
dfEdgeList = pd.DataFrame(columns = ['start_node', 'end_node'])
print(dfRLNode["_merge"].value_counts())
dfRLNode.drop("_merge", axis = 1, inplace = True)

# df1 are teh links that go from minus node to plus node, orientation = +
df1 = dfRLNode.loc[dfRLNode['DirectedLinkOrientation'] == "+"].reindex(columns = ['MinusNodeFID','PlusNodeFID', "RoadLinkFID", "weight"])
# df2 are the nodes that go from plus to minus
df2 = dfRLNode.loc[dfRLNode['DirectedLinkOrientation'] == "-"].reindex(columns = ['PlusNodeFID','MinusNodeFID', "RoadLinkFID", "weight"])
# df3 are the non directed link, part 1
df3 = dfRLNode.loc[dfRLNode['DirectedLinkOrientation'].isnull()].reindex(columns = ['MinusNodeFID','PlusNodeFID', "RoadLinkFID", "weight"])
# df4 are the non directed links part 2
df4 = dfRLNode.loc[dfRLNode['DirectedLinkOrientation'].isnull()].reindex(columns = ['PlusNodeFID','MinusNodeFID', "RoadLinkFID", "weight"])

df1.rename(columns = {'MinusNodeFID':'start_node', 'PlusNodeFID':'end_node'}, inplace=True)
df2.rename(columns = {'PlusNodeFID':'start_node', 'MinusNodeFID':'end_node'}, inplace=True)
df3.rename(columns = {'PlusNodeFID':'start_node', 'MinusNodeFID':'end_node'}, inplace=True)
df4.rename(columns = {'MinusNodeFID':'start_node', 'PlusNodeFID':'end_node'}, inplace=True)

dfEdgeList = pd.concat([df1,df2,df3,df4], join = "outer")
dfEdgeList = dfEdgeList.reindex(columns=['start_node','end_node', 'RoadLinkFID','weight'])


# Now create the network
graphRL = nx.from_pandas_edgelist(dfEdgeList, source='start_node', target = 'end_node', edge_attr = ['RoadLinkFID','weight'], create_using=nx.DiGraph())
print(graphRL.is_directed())

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