# Script for producing the pedestrian network

import json
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import os
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString
import itertools
import re


################################
#
#
# Load data
#
#
# Road network data used to cut pedestrian and vehicle polugons.
# Pedestrian and vehicle polygons
#
################################

projectCRS = {'init' :'epsg:27700'}


with open("config.json") as f:
    config = json.load(f)


gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
output_directory = os.path.join(gis_data_dir, "processed_gis_data")


gdfORLink = gpd.read_file(os.path.join(output_directory, config["openroads_link_processed_file"]))
gdfORNode = gpd.read_file(os.path.join(output_directory, config["openroads_node_processed_file"]))
gdfORLink.crs = projectCRS
gdfORNode.crs = projectCRS

gdfTopoVeh = gpd.read_file(os.path.join(output_directory, config["topo_vehicle_processed_file"]))
gdfTopoPed = gpd.read_file(os.path.join(output_directory, config["topo_pedestrian_processed_file"]))
gdfTopoVeh.crs = projectCRS
gdfTopoPed.crs = projectCRS

# Load boundary data - used to identify traffic island pedestrian polygons
gdfBoundary = gpd.read_file(os.path.join(output_directory, config["boundary_file"]))
gdfBoundary.crs = projectCRS

output_ped_nodes_file = os.path.join(output_directory, "pedNetworkNodes.shp")
output_ped_links_file = os.path.join(output_directory, "pedNetworkLinks.shp")

#################################
#
#
# Functions
#
#
#################################
def find_and_remove_edge_connect_ends(graph, edge_attribute, edge_value):
	edge_to_remove = None
	for u,v,d in graph.edges(data=True):
		if d[edge_attribute] == edge_value:
			edge_to_remove = (u,v)
			break
	return remove_edge_connect_ends(graph, edge_to_remove)

def remove_edge_connect_ends(graph, edge):
	'''
	'''
	edges_to_add = []
	edges_to_remove = [edge]

	u = edge[0]
	v = edge[1]

	## Find nodes linked to the 'from' node of this edge
	## Create new edges that link these directly to the 'to' node
	for node in graph[v]:
		if node != u:
			data = graph.get_edge_data(v, node)

			# overwrite this edge by removing it and adding replacement than connects to u node
			# id of new edge no longer matches geographic representation of link
			edges_to_remove.append((v,node))
			edges_to_add.append((u, node, data))
	
	# Add these edges to the graph
	graph.add_edges_from(edges_to_add)

	# Remove the edge - need to later clean up the orphan nodes
	graph.remove_edges_from(edges_to_remove)
	return graph

def remove_multiple_edges(graph, edge_attribute, edge_values):
	for edge_value in edge_values:
		graph = find_and_remove_edge_connect_ends(graph, edge_attribute, edge_value)
	return graph

def save_geometries(*geoms, path):
	gdf = gpd.GeoDataFrame({'geometry':geoms})
	gdf.to_file(path)
	return gdf

def intersection_of_multiple_geometries(geoms):
	'''Finds intersection between multiple geometries
	'''
	i = geoms.pop()
	for g in geoms:
		i = i.intersection(g)
	return i

def find_single_road_node_pedestrian_nodes(graph, road_node_id, gdfVehPolys, gdfPedPolys):

	# Method for getting ped nodes for a single junctions
	edges = graph.edges(road_node_id)

	# Get the road link IDs connected to this junction node
	rl_fids = [graph.get_edge_data(*e)['fid'] for e in edges]

	gdfPedNodes = gpd.GeoDataFrame()

	# Iterate over pairs of road link IDs in order to find the pedestrian nodes that lie between these two road links
	for fid_pair in itertools.combinations(rl_fids, 2):
		# Initialise data to go into the geodataframe
		row = {"juncNodeID":road_node_id}

		veh_polys_indices = gdfVehPolys.loc[ gdfVehPolys['roadLinkID'].isin(fid_pair)].index
		ped_polys1_indices = gdfPedPolys.loc[ gdfPedPolys['roadLinkID'] == fid_pair[0]].index
		ped_polys2_indices = gdfPedPolys.loc[ gdfPedPolys['roadLinkID'] == fid_pair[1]].index


		# Then form pairs of vehicle polygons and ped polygons
		veh_pairs = list(itertools.combinations(veh_polys_indices, 2))
		ped_pairs = list(itertools.product(ped_polys1_indices, ped_polys2_indices))

		# loop through all combinations of these pairs and find intersection between all 4 polygons. This gives a node in the pedestrian network.
		intersections = []
		for iv1,iv2 in veh_pairs:
			for ip1,ip2 in ped_pairs:
				row['v1pID'] = gdfVehPolys.loc[iv1, 'polyID']
				row['v1rlID'] = gdfVehPolys.loc[iv1, 'roadLinkID']

				row['v2pID'] = gdfVehPolys.loc[iv2, 'polyID']
				row['v2rlID'] = gdfVehPolys.loc[iv2, 'roadLinkID']

				row['p1pID'] = gdfPedPolys.loc[ip1, 'polyID']
				row['p1rlID'] = gdfPedPolys.loc[ip1, 'roadLinkID']
				
				row['p2pID'] = gdfPedPolys.loc[ip2, 'polyID']
				row['p2rlID'] = gdfPedPolys.loc[ip2, 'roadLinkID']

				v1_g = gdfVehPolys.loc[iv1, 'geometry']
				v2_g = gdfVehPolys.loc[iv2, 'geometry']

				p1_g = gdfPedPolys.loc[ip1, 'geometry']
				p2_g = gdfPedPolys.loc[ip2, 'geometry']

				i = intersection_of_multiple_geometries([v1_g,v2_g,p1_g,p2_g])
				if i.is_empty == False:
					row['geometry'] = i
					gdfPedNodes = gdfPedNodes.append(row, ignore_index = True)

	# Filter dataframe to exclude empty geometries
	return gdfPedNodes

def find_multiple_road_node_pedestrian_nodes(graph, road_node_ids, gdfVehPolys, gdfPedPolys):
	gdfPedNodes = gpd.GeoDataFrame()

	for road_node_id in road_node_ids:
		gdf = find_single_road_node_pedestrian_nodes(graph, road_node_id, gdfVehPolys, gdfPedPolys)
		gdfPedNodes = pd.concat([gdfPedNodes, gdf])
		'''
		except Exception as err:
			print(road_node_id)
			print(err)
		'''
	gdfPedNodes.index = np.arange(gdfPedNodes.shape[0])
	return gdfPedNodes

def connect_pavement_ped_nodes(gdfPN, gdfPedPolys, gdfLink, road_graph):

	ped_node_edges = []

	for rl_id in gdfLink['fid'].values:

		# Get start and end node for this road link
		node_records = gdfLink.loc[ gdfLink['fid'] == rl_id, ['MNodeFID','PNodeFID']].to_dict(orient='records')
		assert len(node_records) == 1
		u = node_records[0]['MNodeFID']
		v = node_records[0]['PNodeFID']

		# Get small section of road network connected to this road link
		neighbour_links = [edge_data['fid'] for e, edge_data in road_graph[u].items()]
		neighbour_links += [edge_data['fid'] for e, edge_data in road_graph[v].items()]

		# Get the geometries for these road links
		gdfLinkSub = gdfLink.loc[ gdfLink['fid'].isin(neighbour_links)]

		# Get pairs of ped nodes
		gdfPedNodesSub = gdfPN.loc[(gdfPN['v1rlID']==rl_id) | (gdfPN['v2rlID']==rl_id)]
		ped_node_pairs = itertools.combinations(gdfPedNodesSub['fid'].values, 2)

		for ped_u, ped_v in ped_node_pairs:
			# Create linestring to join ped nodes
			g_u = gdfPedNodesSub.loc[ gdfPedNodesSub['fid'] == ped_u, 'geometry'].values[0]
			g_v = gdfPedNodesSub.loc[ gdfPedNodesSub['fid'] == ped_v, 'geometry'].values[0]

			l = LineString([g_u, g_v])

			# Now check if linestring intersects any of the road link geometries
			if gdfLinkSub['geometry'].map(lambda g: g.intersects(l)).any():
				continue
			else:
				edge_data = {'road_link':None, 'ped_poly':None}

				# Need to identify which pedestrian polygon(s) this edge corresponds to
				candidates = gdfPedPolys.loc[ gdfPedPolys['roadLinkID'] == rl_id, 'polyID'].unique()
				nested_ped_node_polys = gdfPedNodesSub.loc[ gdfPedNodesSub['fid'].isin([ped_u, ped_v]), ['p1pID','p2pID']].to_dict(orient = 'split')['data']
				ped_node_polys = [item for sublist in nested_ped_node_polys for item in sublist]
				edge_data['ped_poly'] = " ".join([p for p in candidates if p in ped_node_polys])

				ped_node_edges.append((ped_u, ped_v, edge_data))

	return ped_node_edges


######################################
#
#
# Create nx network frrom road link data
#
# Graph is not an exact representation of the road link data.
# Some road links do not have any pedestrian polygons associated to them.
# These links should not form part of this network as they will not be traversed by the pedestrians,
# although they represent important geographic information.
#
# To exclude these from the network the output node of these links is replaced by the input node for all
# edges.
#
# Then, network better represents decisions pedestrian agent must make. Need to use this version of the network in pedestrian
# strategic route.
#
#######################################
G = nx.Graph()
gdfORLink['fid_dict'] = gdfORLink['fid'].map(lambda x: {"fid":x})
edges = gdfORLink.loc[:,['MNodeFID','PNodeFID', 'fid_dict']].to_records(index=False)
G.add_edges_from(edges)

# Get fids of links which don't have any pedestrian polygons.
fids_no_ped_polys = gdfORLink.loc[~gdfORLink['fid'].isin(gdfTopoPed['roadLinkID']), 'fid'].values

# Don't remove those that are of a *significant* length - raise assertion error if this is going to happen
gdfORLink['length'] = gdfORLink['geometry'].length
assert gdfORLink.loc[gdfORLink['fid'].isin(fids_no_ped_polys) & (gdfORLink['length'] > 25)].shape[0] == 0

# Sinplify road network such that edges that don't have ped polygons associated to them are removed, and connecting nodes directly linked to source node
G_clean = remove_multiple_edges(G.copy(), 'fid', fids_no_ped_polys)


#################################
#
#
# Get nodes of the pedestrian network
#
#
# Identify junctions in the road network
# Get the pedestrian and vehicle polygons linked to that junctions
# Get intersection between vehicle and ped polygonsfor each road network link of the junction
#
##################################

# Nodes with degree > 2 are junctions. Get these
node_degrees = G_clean.degree()
df_node_degree = pd.DataFrame(node_degrees)
df_node_degree.columns = ['nodeID', 'nodeDegree']

# Identify traffic island pedestrian polygons, want to exclude these nodes from the network
gdfTopoPedBoundary = gpd.sjoin(gdfTopoPed, gdfBoundary, op = 'intersects')
non_boundary_ped_polys = gdfTopoPed.loc[ ~(gdfTopoPed['polyID'].isin(gdfTopoPedBoundary['polyID']))]

island_poly_ids = []

# Find the connected clusters that correspond to connected clusters
# initialise df to recod the neighbours for each ped-vehicle space polygon. use to identify clusters
dfPedNeighbours = pd.DataFrame()
for index, row in gdfTopoPed.iterrows():
    # Touches identifies polygons with at least one point in common but interiors don't intersect. So will work as long as none of my topographic polygons intersect
    neighborFIDs = gdfTopoPed[gdfTopoPed.geometry.touches(row['geometry'])]['polyID'].tolist()

    # Polygons without neighbours need to be attached to themselves so they can be identified as a cluster of one
    if len(neighborFIDs) == 0:
    	neighborFIDs.append(row['polyID'])
    #neighborFIDs = neighborFIDs.remove(row['fid']) polygons don't touch them selves because they intersect
    df = pd.DataFrame({'neighbourid':neighborFIDs})
    df['polyID'] = row['polyID']
    dfPedNeighbours = pd.concat([dfPedNeighbours, df])

g = nx.from_pandas_edgelist(dfPedNeighbours, source='polyID', target = 'neighbourid')
ccs = list(nx.connected_components(g))

# If nodes of connected component all do not intersect a bounday identify these nodes as belonging to a traffic island
for cc in ccs:
	if np.isin(np.array(list(cc)), non_boundary_ped_polys['polyID'].unique()).all():
		island_poly_ids += list(cc)

island_polys = gdfTopoPed.loc[ gdfTopoPed['polyID'].isin(island_poly_ids)]
#island_polys.to_file("island_ped_polys.shp")

# Exclude island polys from the ped network
gdfTopoPed = gdfTopoPed.loc[~(gdfTopoPed['polyID'].isin(island_poly_ids))]
gdfPedNodes = find_multiple_road_node_pedestrian_nodes(G_clean, df_node_degree['nodeID'].values, gdfTopoVeh, gdfTopoPed)

# Remove the multipoints - seems like these are traffic islands - might want to think about including in future
gdfPedNodes = gdfPedNodes.loc[ gdfPedNodes['geometry'].type != 'MultiPoint']

# Create id for each ped node
gdfPedNodes['fid'] = ['pave_node_{}'.format(i) for i in gdfPedNodes.index]
gdfPedNodes.crs = projectCRS
gdfPedNodes.to_file(output_ped_nodes_file)



##############################
#
#
# Create network by connecting nodes
#
# First connect pedestrian nodes that are on the same side of the road.
# - for each road link get ped nodes
# - make pairwise connections if link does not intersect road link (or section of road network +- 1 links from road link
# - handle cases where connection made between nearests this side ped nodes and farther away this side ped node by checking for multiple paths between node pairs.
#
#
##############################


ped_poly_edges = connect_pavement_ped_nodes(gdfPedNodes, gdfTopoPed, gdfORLink, G)

# Now when networkx graph is created duplicated edges will be ignored (since not a MultiGraph)
G_ped_poly = nx.Graph()
G_ped_poly.add_edges_from(ped_poly_edges)


#############################
#
# Connect ped nodes across road links
#
#############################

# Geoup by junctions node id
grouped = gdfPedNodes.groupby('juncNodeID')
group_names = list(grouped.groups.keys())
group_sizes = grouped.apply(lambda df: df.shape[0])
group_sizes = grouped.transform(lambda df: df.shape[0])


# Write process for connecting nodes in a group
# join ped nodes that have the same road link id
pcol_re = re.compile(r'p\dpID.*')

def check_ped_poly_id_repreated(row, ped_poly_col_regex = pcol_re):
	ped_poly_cols = [i for i in row.index if pcol_re.search(i) is not None]
	return row[ped_poly_cols].dropna().duplicated().any()

def connect_junction_ped_nodes(df, ped_node_col, v1_poly_col, v2_poly_col):

	# Connect nodes if they share a road link ID since this means they lie on the opposite side of the same road
	# This process also ends up connecting nodes to themselves which is corrected for afterwards
	junc_edges1 = pd.merge(df, df, on = v1_poly_col, how = 'inner', suffixes = ('_to', '_from'))
	junc_edges1['road_link'] = junc_edges1[v1_poly_col]

	junc_edges2 = pd.merge(df, df, left_on = v1_poly_col, right_on = v2_poly_col, how = 'inner', suffixes = ('_to', '_from'))
	junc_edges2['road_link'] = junc_edges2[v1_poly_col+'_to']

	junc_edges3 = pd.merge(df, df, on = v2_poly_col, how = 'inner', suffixes = ('_to', '_from'))
	junc_edges3['road_link'] = junc_edges3[v2_poly_col]

	junc_edges = pd.concat([junc_edges1, junc_edges2, junc_edges3])

	junc_edges['edge_data'] = junc_edges.apply(lambda row: {"road_link":row["road_link"],"ped_poly":None}, axis=1)

	return junc_edges

junc_edges = gdfPedNodes.groupby('juncNodeID').apply(connect_junction_ped_nodes, 'fid','v1rlID', 'v2rlID')

# Drop duplicates, drop self loops and recreate index
junc_edges.drop_duplicates(subset=['fid_from','fid_to'], inplace = True)
junc_edges = junc_edges.loc[ junc_edges['fid_from'] != junc_edges['fid_to']]
junc_edges.index = np.arange(len(junc_edges))

# Where the link connects nodes on the same ped poly, remove reference to road link as in these cases link does not cross road link
junc_edges['share_ped_poly'] = junc_edges.apply(check_ped_poly_id_repreated, axis = 1)

# Commented out because I realised that it might cross the link at one end of the ped polygon, even if it doesn't cross at another end.
'''
for i in junc_edges.loc[ junc_edges['share_ped_poly'] == True].index:
	row = junc_edges.loc[i]
	row['edge_data'] = {'road_link':None, 'ped_poly':None}
	junc_edges.loc[i] = row
'''

junc_edges['edge'] = junc_edges.apply(lambda row: (row['fid_from'], row['fid_to'], row['edge_data']), axis=1)


# Add these edges to the network
G_junc = nx.Graph()
G_junc.add_edges_from(junc_edges['edge'].values)
G_ped_poly.add_edges_from(junc_edges['edge'].values)


###################################
#
#
# LineString gdf of ped network
#
#
####################################
edge_data = []
for u,v,data in list(G_ped_poly.edges(data=True)):
	row = {'from_node':u, 'to_node':v,'road_link':data['road_link'], 'ped_poly':data['ped_poly']}
	edge_data.append(row)
dfPedNetwork = pd.DataFrame(edge_data)

# Now join with node coordinates
dfPedNetwork = pd.merge(dfPedNetwork, gdfPedNodes.reindex(columns = ['fid','geometry']), left_on = 'from_node', right_on = 'fid', how = 'left', indicator = True)
assert dfPedNetwork.loc[ dfPedNetwork['_merge'] != 'both'].shape[0]==0
dfPedNetwork.drop(['_merge','fid'], axis = 1, inplace = True)
dfPedNetwork.rename(columns = {'geometry':'geometry_from'}, inplace = True)

dfPedNetwork = pd.merge(dfPedNetwork, gdfPedNodes.reindex(columns = ['fid','geometry']), left_on = 'to_node', right_on = 'fid', how = 'left', suffixes = ('','_to'), indicator = True)
assert dfPedNetwork.loc[ dfPedNetwork['_merge'] != 'both'].shape[0]==0
dfPedNetwork.drop(['_merge','fid'], axis = 1, inplace = True)
dfPedNetwork.rename(columns = {'geometry':'geometry_to'}, inplace = True)

# Create LineString connecting edge nodes
dfPedNetwork['geometry'] = dfPedNetwork.apply(lambda row: LineString([row['geometry_from'], row['geometry_to']]), axis = 1)
dfPedNetwork = dfPedNetwork.reindex(columns = ['from_node','to_node','road_link','ped_poly','geometry'])
gdfPedNetwork = gpd.GeoDataFrame(dfPedNetwork, geometry = 'geometry')

# Rename node columns to match other network data columns and create fid field
gdfPedNetwork.rename(columns = {"from_node":"MNodeFID", "to_node":"PNodeFID", "road_link":"pedRLID", "ped_poly":"pedRoadID"}, inplace=True)
gdfPedNetwork['fid'] = "pave_link" + gdfPedNetwork['MNodeFID'].str.replace("pave_node","") + gdfPedNetwork['PNodeFID'].str.replace("pave_node","")

gdfPedNetwork.crs = projectCRS
gdfPedNetwork.to_file(output_ped_links_file)

# Also select those edges that connect ped nodes that share a polygon from junc edges
gdfPedNetwork['edge'] = gdfPedNetwork.apply(lambda row: (row['MNodeFID'],row['PNodeFID']), axis=1)
junc_edges['edgea'] = junc_edges.apply(lambda row: (row['fid_from'],row['fid_to']), axis=1)
junc_edges['edgeb'] = junc_edges.apply(lambda row: (row['fid_to'],row['fid_from']), axis=1)


# Check those edges that connect nodes around a junction that are on the same ped poly - just for checking
share_poly_edges = list(junc_edges.loc[(junc_edges['share_ped_poly'] == True) & (junc_edges['fid_from'] != junc_edges['fid_to']), 'edgea'].values) + list(junc_edges.loc[(junc_edges['share_ped_poly'] == True) & (junc_edges['fid_from'] != junc_edges['fid_to']), 'edgeb'].values)
gdfPedNetSharePoly = gdfPedNetwork.loc[ gdfPedNetwork['edge'].isin(share_poly_edges)]
#gdfPedNetSharePoly.drop(['edge'], axis=1).to_file("pedNetworkSharePoly.shp")
