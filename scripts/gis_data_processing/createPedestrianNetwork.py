# Script for producing the pedestrian network

import json
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import os
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString
import itertools


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


# Create nx network frrom road link data
G = nx.Graph()
gdfORLink['fid_dict'] = gdfORLink['fid'].map(lambda x: {"fid":x})
edges = gdfORLink.loc[:,['startNode','endNode', 'fid_dict']].to_records(index=False)
G.add_edges_from(edges)


#################################
#
#
# Functions
#
#
#################################
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
		row = {"junc_node":road_node_id}

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
				row['v1_polyID'] = gdfVehPolys.loc[iv1, 'polyID']
				row['v1_roadLinkID'] = gdfVehPolys.loc[iv1, 'roadLinkID']

				row['v2_polyID'] = gdfVehPolys.loc[iv2, 'polyID']
				row['v2_roadLinkID'] = gdfVehPolys.loc[iv2, 'roadLinkID']

				row['p1_polyID'] = gdfPedPolys.loc[ip1, 'polyID']
				row['p1_roadLinkID'] = gdfPedPolys.loc[ip1, 'roadLinkID']
				
				row['p2_polyID'] = gdfPedPolys.loc[ip2, 'polyID']
				row['p2_roadLinkID'] = gdfPedPolys.loc[ip2, 'roadLinkID']

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

	return gdfPedNodes



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
node_degrees = G.degree()
df_node_degree = pd.DataFrame(node_degrees)
df_node_degree.columns = ['nodeID', 'nodeDegree']
nodes_twopl = df_node_degree.loc[df_node_degree['nodeDegree'] > 2, 'nodeID']
gdfORNodeJuncs = gdfORNode.loc[gdfORNode['node_fid'].isin(nodes_twopl)]

junc_node = nodes_twopl.values[0]

gdfPedNodes1 = find_junction_pedestrian_nodes(junc_node, gdfTopoVeh, gdfTopoPed)
gdfPedNodes1.crs = projectCRS
gdfPedNodes1.to_file("pedNodes1.shp")