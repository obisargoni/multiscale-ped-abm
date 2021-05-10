# Script for producing the pedestrian network using a road network and walkable bounday

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


gis_data_dir = config['gis_data_dir']
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

output_ped_nodes_file = os.path.join(output_directory, config["pavement_nodes_file"])
output_ped_links_file = os.path.join(output_directory, config["pavement_links_file"])


#################################
#
#
# Functions
#
#
#################################
def road_node_pedestrian_nodes_metadata(graph, road_node_id):

	# Method for getting ped nodes for a single junctions
	edges = graph.edges(road_node_id, data=True)

	dfPedNodes = pd.DataFrame()

	# Iterate over pairs of road link IDs in order to find the pedestrian nodes that lie between these two road links
	for e1, e2 in itertools.combinations(edges, 2):
		# Initialise data to go into the geodataframe
		row = {"juncNodeID":road_node_id, "rlID1":e1[-1]['fid'], "rlID2":e2[-1]['fid']}

		dfPedNodes = dfPedNodes.append(row, ignore_index = True)
		
	# Filter dataframe to exclude empty geometries
	return dfPedNodes

def multiple_road_node_pedestrian_nodes_metadata(graph, road_node_ids):
	dfPedNodes = pd.DataFrame()

	for road_node_id in road_node_ids:
		df = road_node_pedestrian_nodes_metadata(graph, road_node_id)
		dfPedNodes = pd.concat([dfPedNodes, df])

	gdfPedNodes.index = np.arange(gdfPedNodes.shape[0])
	return gdfPedNodes

def assign_boundary_coordinates_to_ped_nodes(dfPN, gdfRoadLinks, gdfBounds):
	"""Identify coordinates for ped nodes based on the bounday.
	"""

	# Loop through nodes, get corresponding road links and the boundary between them
	for ix, row in dfPN.iterrows():

		rlID1 = row['rlID1']
		rlID2 = row['rlID2']

		rlg1 = gdfRoadLinks.loc[ gdfRoadLinks['fid'] == rlID1]
		rlg2 = gdfRoadLinks.loc[ gdfRoadLinks['fid'] == rlID2]

		road_node = rlg1.intersection(rlg2)
		assert isinstance(road_node, Point)
		road_node = np.array(road_node.coords)

		if road_node.coords == rlg1.coords[0]:
			rlp1 = np.array(rlg1.coords[-1])
		else:
			rlp1 = np.array(rlg1.coords[0])

		if road_node.coords == rlg2.coords[0]:
			rlp2 = np.array(rlg2.coords[-1])
		else:
			rlp2 = np.array(rlg2.coords[0])

		v1 = rlp1 - road_node
		v2 = rlp2 - road_node

		# Find boundary between using ray casting
		a1 = angle_between_north_and_vector(v1)
		a2 = angle_between_north_and_vector(v2)

		intersecting_boundary_geom_ids = np.array([])
		for l in rays_between_angles(a1, a2):
			ids = gdfBounds.loc[ gdfBounds.intersects(l)].index
			intersecting_boundary_geom_ids = np.concat([intersecting_boundary_geom_ids, ids])


		# Now find nearest boundary coordinate from intersecting boundaries
		min_dist = sys.maxsize
		ped_node_geom = None
		for row_id in intersecting_boundary_geom_ids:
			geom = gdfBounds.loc[row_id, 'geometry']

			for c in geom.coords:
				rn = Point(road_node)
				p = Point(c)
				d = rn.distance(p)
				if d < min_dist:
					min_dist = d
					ped_node_geom = p

		row['geometry'] = ped_node_geom


def rays_between_angles(a1, a2, sample_res = 10, ray_length = 50):
	sample_res = (2*np.pi) * (sample_res/360.0) # 10 degrees in rad
	sampled_angles = range(a1,a2,sample_res)
	for sa in sampled_angles:
		p1 = Point(road_node)
		p2 = Point([p1.x + ray_length*np.sin(sa), p1.y + ray_length*np.cos(sa)])
		l = LineString([p1,p2])
		yield l

def unit_vector(v):
	magV = np.linalg.norm(v)
	return v / magV

def angle_between_north_and_unit_vector(u):
	n = (0,1)

	signX = np.sign(u[0])
	signY = np.sign(u[1])

	dp = np.dot(n, u)
	a = np.acos(dp)

    # Dot product gives angle between vectors. We need angle clockwise from north
    if (signX == 1):
    	return a
    elif (signX == -1):
    	return 2*np.pi - a
    else:
    	return a

def angle_between_north_and_vector(v):
 	unit_v = unit_vector(v)
 	return angle_between_north_and_unit_vector(unit_v)
			

def ang(lineA, lineB):
    # Get nicer vector form
    lACoords = lineA.coords
    lBCoords = lineB.coords
    vA = [(lACoords[0][0]-lACoords[1][0]), (lACoords[0][1]-lACoords[1][1])]
    vB = [(lBCoords[0][0]-lBCoords[1][0]), (lBCoords[0][1]-lBCoords[1][1])]
    # Get dot prod
    dot_prod = np.dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = round(dot_prod/magA/magB, 4)
    # Get angle in radians and then convert to degrees
    angle = math.acos(cos_)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        # As in if statement
        return 360 - ang_deg
    else: 

        return ang_deg
################################
#
#
# Produce nodes metadata
#
# For every node in the road network, create pedestrian nodes in the regions between the links in that network
#
################################

# Load the Open Roads road network as a nx graph
G = nx.MultiGraph()
gdfORLink['fid_dict'] = gdfORLink['fid'].map(lambda x: {"fid":x})
edges = gdfORLink.loc[:,['MNodeFID','PNodeFID', 'fid_dict']].to_records(index=False)
G.add_edges_from(edges)


# Node metadata
dfPedNodes = multiple_road_node_pedestrian_nodes_metadata(G, df_node_degree['nodeID'].values)


##################################
#
#
# Assign coordinates to nodes
#
# Do I need to consider a connected component of walkable surfaces?
#
#
##################################

gdfPedNodes = assign_boundary_coordinates_to_ped_nodes(dfPedNodes, gdfORLink, gdfBoundary)

# alternatively could have assign_ped_veh_polygon_coordaintes_to_ped_nodes()

# Remove the multipoints - seems like these are traffic islands - might want to think about including in future
gdfPedNodes = gdfPedNodes.loc[ gdfPedNodes['geometry'].type != 'MultiPoint']

# Create id for each ped node
gdfPedNodes['fid'] = ['pave_node_{}'.format(i) for i in gdfPedNodes.index]
gdfPedNodes.crs = projectCRS

gdfPedNodes.to_file(output_ped_nodes_file)







