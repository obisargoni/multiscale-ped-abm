# Script for producing the pedestrian network using a road network and walkable bounday
import sys
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

# Load boundary data - used to identify pavement nodes
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
def unit_vector(v):
    magV = np.linalg.norm(v)
    return v / magV

def angle_between_north_and_unit_vector(u):
    n = (0,1)

    signX = np.sign(u[0])
    signY = np.sign(u[1])

    dp = np.dot(n, u)
    a = np.arccos(dp)

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
    angle = math.arccos(cos_)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        # As in if statement
        return 360 - ang_deg
    else: 
        return ang_deg


def sample_angles(a1, a2, sample_res):
    if a1 < a2:
        sampled_angles = np.arange(a1+sample_res, a2,sample_res)
    else:
        sampled_angles = []
        ang = a1+sample_res
        while ang < 2*np.pi:
            sampled_angles.append(ang)
            ang+=sample_res
        
        ang-=2*np.pi
        
        while ang < a2:
            sampled_angles.append(ang)
            ang+=sample_res

        sampled_angles = np.array(sampled_angles)

    return sampled_angles

def in_angle_range(ang, a1, a2):
    if a1 < a2:
        return (ang>a1) & (ang<a2)
    else:
        b1 = (ang>a1) & (ang<0)
        b2 = (ang>=0) & (ang<a2)
        return b1 | b2

def rays_between_angles(a1, a2, p1, sample_res = 10, ray_length = 50):
    sample_res = (2*np.pi) * (sample_res/360.0) # 10 degrees in rad
    sampled_angles = sample_angles(a1, a2, sample_res)
    for sa in sampled_angles:
        p2 = Point([p1.x + ray_length*np.sin(sa), p1.y + ray_length*np.cos(sa)])
        l = LineString([p1,p2])
        yield l

def linestring_bearing(l, start_point):
    if start_point.coords[0] == l.coords[0]:
        end_coord = np.array(l.coords[-1])
    elif start_point.coords[0] == l.coords[-1]:
        end_coord = np.array(l.coords[0])
    else:
        return None

    start_coord = np.array(start_point.coords[0])

    v = end_coord - start_coord

    return angle_between_north_and_vector(v)

def road_node_pedestrian_nodes_metadata(graph, road_node_geom, road_node_id):

    # Method for getting ped nodes for a single junctions
    edges = list(graph.edges(road_node_id, data=True))

    # Find neighbouring edges based on the edge geometry bearing
    edge_bearings = [linestring_bearing(e[-1]['geometry'], road_node_geom) for e in edges]

    edges_w_bearing = list(zip(edges, edge_bearings))
    edges_w_bearing.sort(key = lambda e: e[1])
    edges, edge_bearings = zip(*edges_w_bearing)
    edge_pairs = zip(edges, edges[1:] + edges[0:1])
    
    # Iterate over pairs of road link IDs in order to find the pedestrian nodes that lie between these two road links
    dfPedNodes = pd.DataFrame()
    for e1, e2 in edge_pairs:
        # Initialise data to go into the geodataframe
        row = {"juncNodeID":road_node_id, "rlID1":e1[-1]['fid'], "rlID2":e2[-1]['fid']}

        dfPedNodes = dfPedNodes.append(row, ignore_index = True)
        
    # Filter dataframe to exclude empty geometries
    return dfPedNodes

def multiple_road_node_pedestrian_nodes_metadata(graph, gdfRoadNodes):
    dfPedNodes = pd.DataFrame()

    for road_node_id, road_node_geom in gdfRoadNodes.loc[:, ['node_fid', 'geometry']].values:
        df = road_node_pedestrian_nodes_metadata(graph, road_node_geom, road_node_id)
        dfPedNodes = pd.concat([dfPedNodes, df])

    dfPedNodes.index = np.arange(dfPedNodes.shape[0])
    return dfPedNodes

def nearest_geometry_point_between_angles(a1, a2, start_point, seriesGeoms):
        intersecting_boundary_geom_ids = np.array([])
        for l in rays_between_angles(a1, a2, start_point):
            ids = seriesGeoms.loc[ seriesGeoms.intersects(l)].index
            intersecting_boundary_geom_ids = np.concatenate([intersecting_boundary_geom_ids, ids])


        # Now find nearest boundary coordinate from intersecting boundaries
        min_dist = sys.maxsize
        nearest_point = None
        for row_id in set(intersecting_boundary_geom_ids):
            geom = seriesGeoms.loc[row_id]

            for c in geom.exterior.coords:
                p = Point(c)
                d = start_point.distance(p)

                # ensure that point lies in direction between input and output angles
                a = linestring_bearing(LineString([start_point, p]), start_point)

                if (d < min_dist) & (in_angle_range(a, a1, a2)):
                    min_dist = d
                    nearest_point = p
        
        return nearest_point

def assign_boundary_coordinates_to_ped_nodes(dfPN, gdfRoadLinks, serBounds, crs = projectCRS):
    """Identify coordinates for ped nodes based on the bounday.
    """

    # Initialise output
    gdfPN = gpd.GeoDataFrame(dfPN)
    gdfPN['geometry'] = None
    gdfPN.set_geometry('geometry', inplace=True, crs = crs)

    # Loop through nodes, get corresponding road links and the boundary between them
    for ix, row in gdfPN.iterrows():

        rlID1 = row['rlID1']
        rlID2 = row['rlID2']

        rlg1 = gdfRoadLinks.loc[ gdfRoadLinks['fid'] == rlID1, 'geometry'].values[0]
        rlg2 = gdfRoadLinks.loc[ gdfRoadLinks['fid'] == rlID2, 'geometry'].values[0]

        road_node = rlg1.intersection(rlg2)
        if isinstance(road_node, Point)==False:
            print("Non-point intersection. {} {} {}".format(ix, rlID1, rlID2))
            print(road_node.type)
            print("\n") 
            continue
        else:
            pass

        # Find boundary between using ray casting
        a1 = linestring_bearing(rlg1, road_node)
        a2 = linestring_bearing(rlg2, road_node)

        ped_node_geom = nearest_geometry_point_between_angles(a1, a2, road_node, serBounds)

        gdfPN.at[ix, 'geometry'] = ped_node_geom

    return gdfPN


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
gdfORLink['fid_dict'] = gdfORLink.apply(lambda x: {"fid":x['fid'],'geometry':x['geometry']}, axis=1)
edges = gdfORLink.loc[:,['MNodeFID','PNodeFID', 'fid_dict']].to_records(index=False)
G.add_edges_from(edges)

# Node metadata
dfPedNodes = multiple_road_node_pedestrian_nodes_metadata(G, gdfORNode)


##################################
#
#
# Assign coordinates to nodes
#
# Do I need to consider a connected component of walkable surfaces?
#
#
##################################

# Buffer the boundary so that nodes are located slightly away from walls
gdfBoundaryBuff = gdfBoundary.buffer(0.5)

gdfPedNodes = assign_boundary_coordinates_to_ped_nodes(dfPedNodes, gdfORLink, gdfBoundaryBuff)

# alternatively could have assign_ped_veh_polygon_coordaintes_to_ped_nodes()

# Remove the multipoints - seems like these are traffic islands - might want to think about including in future
'''
gdfPedNodes = gdfPedNodes.loc[ gdfPedNodes['geometry'].type != 'MultiPoint']

# Create id for each ped node
gdfPedNodes['fid'] = ['pave_node_{}'.format(i) for i in gdfPedNodes.index]
gdfPedNodes.crs = projectCRS

gdfPedNodes.to_file(output_ped_nodes_file)
'''