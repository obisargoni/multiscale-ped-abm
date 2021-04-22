'''
This is a script for editing OS's ITN Road Network shapefiles so that:
- the attributeds include an id for the 'to' and 'from' nodes
- line strings are duplicated along links that are bidirectional

Purpose of this script is to use the extracted orientation information from the gml data to edit the roads linesting shape file such that 
two way road are duplicated with one linestring orientated in one direction and the other orientated in another direction.
'''

import json
import numpy as np
import pandas as pd
import geopandas as gpd
import os
import networkx as nx
import osmnx
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString

######################
#
#
# Functions
#
#
######################
import math

def make_linestring_coords_2d(l):
    new_coords = []
    for c in l.coords:
        if len(c) > 1:
            new_coords.append((c[0],c[1]))
    return LineString(new_coords)

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]
def ang(lineA, lineB):
    # Get nicer vector form
    lACoords = lineA.coords
    lBCoords = lineB.coords
    vA = [(lACoords[0][0]-lACoords[1][0]), (lACoords[0][1]-lACoords[1][1])]
    vB = [(lBCoords[0][0]-lBCoords[1][0]), (lBCoords[0][1]-lBCoords[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
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

def group(lst, n, overlapp = 1):
    '''Group a list in to groups of size n, with overlap number of iterms overlapping
    '''
    for i in range(0, len(lst), n - overlapp):
        val = lst[i:i+n]
        if len(val) == n:
            yield tuple(val)

def split_line(l, npoints = 2):
    '''Takes a list string geometry and splits it into subline strings made up of npoints.
    Useful for simplifying a linestring geometry
    '''
    coord_pairs = group(l.coords, npoints, overlapp = 1)
    for cpair in coord_pairs:
        yield LineString(cpair)

def simplify_line_angle(l, angle_threshold = 10):
    '''Break a linestring into components such that the anglular distance along each component 
    is below the threshold. Also simplify the linestrings to be two coordinates only. Cleans road network so that each line string
    acts as maximal angular deviation unit of angular distance
    
    Might want to separate these two cleaning operations
    '''
    simplified_lines = []
    
    split_lines = split_line(l, 2)

    la = next(split_lines)
    lb = None
    c1 = la.coords[0]
    angle = 0.0

    for lb in split_lines:
        angle+=abs(ang(la,lb))

        # If angle has reached threshold create simplified linestring
        if angle >= angle_threshold:
            c2 = la.coords[-1]
            simplified_lines.append(LineString((c1,c2)))
            c1 = c2
            angle = 0
        la = lb

    c2 = lb.coords[-1]
    simplified_lines.append(LineString((c1,c2)))

    return simplified_lines

# Disolved geometries are multi polygons, explode to single polygons
def simplify_line_gdf_by_angle(indf, angle_threshold, id_col, new_id_col):
    outdf = gpd.GeoDataFrame(columns=indf.columns)
    outdf[new_id_col] = np.nan
    for idx, row in indf.iterrows():
        l = row['geometry']
        if len(l.coords) == 2:
            row[new_id_col] = row[id_col]
            outdf = outdf.append(row,ignore_index=True)
        else:
            simplified_lines = simplify_line_angle(l, angle_threshold)
            multdf = gpd.GeoDataFrame(columns=indf.columns)
            nlines = len(simplified_lines)
            multdf = multdf.append([row]*nlines,ignore_index=True)
            for i in range(nlines):
                multdf.loc[i,'geometry'] = simplified_lines[i]
                multdf.loc[i, new_id_col] = multdf.loc[i, id_col] + "_{}".format(i)
            outdf = outdf.append(multdf,ignore_index=True)
    return outdf

def largest_connected_component_nodes_within_dist(G, source_node, dist, weight):
    lccNodes = max(nx.connected_components(G), key=len)

    lccG = G.subgraph(lccNodes).copy()
    
    shortest_paths = nx.single_source_dijkstra(lccG, source_node, target=None, cutoff=dist, weight=weight)   
    reachable_nodes = shortest_paths[0].keys()

    return reachable_nodes

def graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True):
    """
    Convert a Undirected Graph to node and/or edge GeoDataFrames.
    This function is the inverse of `graph_from_gdfs`.
    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    nodes : bool
        if True, convert graph nodes to a GeoDataFrame and return it
    edges : bool
        if True, convert graph edges to a GeoDataFrame and return it
    node_geometry : bool
        if True, create a geometry column from node x and y attributes
    fill_edge_geometry : bool
        if True, fill in missing edge geometry fields using nodes u and v
    Returns
    -------
    geopandas.GeoDataFrame or tuple
        gdf_nodes or gdf_edges or tuple of (gdf_nodes, gdf_edges). gdf_nodes
        is indexed by osmid and gdf_edges is multi-indexed by u, v, key
        following normal MultiDiGraph structure.
    """
    crs = G.graph["crs"]

    if nodes:

        if not G.nodes:  # pragma: no cover
            raise ValueError("graph contains no nodes")

        nodes, data = zip(*G.nodes(data=True))

        if node_geometry:
            # convert node x/y attributes to Points for geometry column
            geom = (Point(d["x"], d["y"]) for d in data)
            gdf_nodes = gpd.GeoDataFrame(data, index=nodes, crs=crs, geometry=list(geom))
        else:
            gdf_nodes = gpd.GeoDataFrame(data, index=nodes)

        gdf_nodes.index.rename("osmid", inplace=True)
        osmnx.utils.log("Created nodes GeoDataFrame from graph")

    if edges:

        if not G.edges:  # pragma: no cover
            raise ValueError("graph contains no edges")

        if isinstance(G, nx.classes.multidigraph.MultiDiGraph):
            u, v, k, data = zip(*G.edges(keys=True, data=True))
        else:
            u, v, data = zip(*G.edges(data=True))
            k = [0]*len(u)


        if fill_edge_geometry:

            # subroutine to get geometry for every edge: if edge already has
            # geometry return it, otherwise create it using the incident nodes
            x_lookup = nx.get_node_attributes(G, "x")
            y_lookup = nx.get_node_attributes(G, "y")

            def make_geom(u, v, data, x=x_lookup, y=y_lookup):
                if "geometry" in data:
                    return data["geometry"]
                else:
                    return LineString((Point((x[u], y[u])), Point((x[v], y[v]))))

            geom = map(make_geom, u, v, data)
            gdf_edges = gpd.GeoDataFrame(data, crs=crs, geometry=list(geom))

        else:
            gdf_edges = gpd.GeoDataFrame(data)
            if "geometry" not in gdf_edges.columns:
                # if no edges have a geometry attribute, create null column
                gdf_edges["geometry"] = np.nan
            gdf_edges.set_geometry("geometry")
            gdf_edges.crs = crs

        # add u, v, key attributes as index
        gdf_edges["u"] = u
        gdf_edges["v"] = v
        gdf_edges["key"] = k
        gdf_edges.set_index(["u", "v", "key"], inplace=True)

        osmnx.utils.log("Created edges GeoDataFrame from graph")

    if nodes and edges:
        return gdf_nodes, gdf_edges
    elif nodes:
        return gdf_nodes
    elif edges:
        return gdf_edges
    else:  # pragma: no cover
        raise ValueError("you must request nodes or edges or both")


def nodes_gdf_from_edges_gdf(gdf_edges, u, v):
    '''Given a geo data frame of network edges with LineString geometries, etract the start and end points of the 
    LineString geometries and create a nodes geo data frame from these.

    Set the u and v columns of the edges geo data frame to the corresponding node ids

    ----------
    edges_gdf : geopandas.GeoDataFrame
        input edges
    u : str
        column to use for u node id
    v : str
        column to use for v node id
    Returns
    -------
    tuple
        (gdf_nodes, gdf_edges)
    '''
    nodes_data = {'id':[], 'geometry':[]}
    u_data = []
    v_data = []

    node_id = 0
    for row_index, row in gdf_edges.iterrows():
        g = row['geometry']

        nodes_data['id'].append(node_id)
        nodes_data['geometry'].append(Point(g.coords[0]))
        u_data.append(node_id)

        node_id+=1

        nodes_data['id'].append(node_id)
        nodes_data['geometry'].append(Point(g.coords[-1]))
        v_data.append(node_id)

        node_id+=1

    gdf_nodes = gpd.GeoDataFrame(nodes_data)
    gdf_nodes.crs = gdf_edges.crs

    gdf_edges[u] = u_data
    gdf_edges[v] = v_data

    return gdf_nodes, gdf_edges

######################
#
#
# Initialise variables and paths to data inputs and outputs
#
#
#####################
projectCRS = "epsg:27700"

with open("config.json") as f:
    config = json.load(f)

gis_data_dir = config['gis_data_dir']

itn_directory = os.path.join(gis_data_dir, config['mastermap_itn_name'])
itn_link_file = os.path.join(itn_directory, "mastermap-itn RoadLink", "mastermap-itn RoadLink.shp")
itn_node_file = os.path.join(itn_directory, "mastermap-itn RoadNode", "mastermap-itn RoadNode.shp")

open_roads_directory = os.path.join(gis_data_dir, config['open_roads_dir'])
or_link_file = os.path.join(open_roads_directory, config['open_roads_link_file'])
or_node_file = os.path.join(open_roads_directory, config['open_roads_node_file'])

poi_file = os.path.join(gis_data_dir, config['poi_file'])

output_directory = os.path.join(gis_data_dir, "processed_gis_data")

if os.path.isdir(output_directory) == False:
    os.mkdir(output_directory)

selection_layer_file = os.path.join(gis_data_dir, config['clip_file'])

output_itn_link_file = os.path.join(output_directory, config["mastermap_link_processed_file"])
output_itn_node_file = os.path.join(output_directory, config["mastermap_node_processed_file"])

output_or_link_file = os.path.join(output_directory, config["openroads_link_processed_file"])
output_or_node_file = os.path.join(output_directory, config["openroads_node_processed_file"])


###########################
#
#
# Load Data
# 
#
##########################

# Mastermap ITN data - for road network
gdfITNLink = gpd.read_file(itn_link_file)
gdfITNLink.crs = projectCRS

gdfITNNode = gpd.read_file(itn_node_file)
gdfITNNode.crs = projectCRS

# OS Open Road - for ped road network
gdfORLink = gpd.read_file(or_link_file)
gdfORLink.crs = projectCRS

gdfORNode = gpd.read_file(or_node_file)
gdfORNode.crs = projectCRS

gdfPOIs = gpd.read_file(poi_file)

# Study area polygon - to select data within the study area
'''
gdfSelect = gpd.read_file(selection_layer_file)
gdfSelect.crs = projectCRS
SelectPolygon = gdfSelect.loc[0,'geometry']
'''
centre_poi = gdfPOIs.loc[gdfPOIs['ref_no'] == config['centre_poi_ref']] 
centre_poi_geom = centre_poi['geometry'].values[0]

seriesStudyArea = centre_poi.buffer(3000)
seriesStudyArea.to_file(os.path.join(gis_data_dir, "study_area.shp"))
gsStudyAreaWSG84 = seriesStudyArea.to_crs(epsg=4326)

studyPolygon = seriesStudyArea.values[0]
studyPolygonWSG84 = gsStudyAreaWSG84.values[0]

################################
#
# Select the ITN Road Network that lies in the study area
#
# Need to use study polygon rather than network distance method because at this stage the ITN network has not been created.
# This is done in the script makeITNdirectional.py
#
################################

# Select only the polygons that intersect or lie within the junc clip area
gdfITNLink = gdfITNLink.loc[ (gdfITNLink.geometry.intersects(studyPolygon)) | (gdfITNLink.geometry.within(studyPolygon))]
gdfITNNode = gpd.sjoin(gdfITNNode, gdfITNLink.loc[:,['fid','geometry']], op = 'intersects', lsuffix = 'node', rsuffix = 'line')


# Clean up
gdfITNNode.drop(['fid_line', 'index_line'], axis = 1, inplace=True)
gdfITNNode.drop_duplicates(inplace = True)
gdfITNNode.rename(columns = {'fid_node':'fid'}, inplace = True)


##############################
#
# Select the OS Open Road data that lies in the study area
#
# OS Open Road Data is used for pedestrian routing. Less detail on lanes and roundabouts so more suitable for peds
#
##############################

# Get into format required for osmnx compliant graph
gdfORLink = gdfORLink.rename(columns = {"identifier": "osmid", 'startNode':'u', 'endNode':'v'})
gdfORLink['key'] = 0
gdfORLink.set_index(['u','v','key'], inplace=True)
gdfORLink['geometry'] = gdfORLink['geometry'].map(make_linestring_coords_2d)

gdfORNode = gdfORNode.rename(columns = {"identifier": "osmid"})
gdfORNode['x'] = gdfORNode.loc[:, 'geometry'].map(lambda g: g.x)
gdfORNode['y'] = gdfORNode.loc[:, 'geometry'].map(lambda g: g.y)
gdfORNode.set_index('osmid', inplace=True)

# Get largest connected component
'''
edges = gdfORLink.loc[:,['startNode','endNode','length']].values
edges_ids = gdfORLink['fid'].values
G = nx.Graph()
G.add_weighted_edges_from(edges, weight='length', fid = edges_ids)
'''

# Makes sense to set up graph as osmnx compliant object. But need to make sure I can keep track of edge ids
G = osmnx.graph_from_gdfs(gdfORNode, gdfORLink, graph_attrs=None)

# Convert to undirected
U = G.to_undirected()

# Find the or node nearest the centre poi
gdfORNode['dist_to_centre'] = gdfORNode.distance(centre_poi_geom)
nearest_node_id = gdfORNode.sort_values(by = 'dist_to_centre', ascending=True).index[0]

reachable_nodes = largest_connected_component_nodes_within_dist(U, nearest_node_id, 2500, 'length')

U = U.subgraph(reachable_nodes)
G = U.to_directed() # osmnx expected MultiDiGraph. Setting to directed from undirected should maintain undirected nnature but make this explicit

###################################
#
#
# Now that study area network has been selected, clean network by:
# - simplify intersections
# - split lines based on angular deviation
#
####################################

# simplify intersections - not working for some reason
G = osmnx.simplification.consolidate_intersections(G, tolerance=10, rebuild_graph=True, dead_ends=True, reconnect_edges=True)


# Convert to undirected for next bit of cleaning. Keep multi edge representation though. Need to think about this - but think it makes sense to retain most general structure
U = G.to_undirected()
gdfORNode, gdfORLink = osmnx.graph_to_gdfs(U)

gdfORLink.reset_index(inplace=True)
gdfORNode.reset_index(inplace=True)


# Clean data to ensure minimum angular deviation along road link
assert gdfORLink['geometry'].type.unique().size == 1
assert gdfORLink['osmid'].unique().size == gdfORLink.shape[0]


gdfORLink_simplified = simplify_line_gdf_by_angle(gdfORLink, 10, "osmid", "fid")

assert gdfORLink_simplified['fid'].duplicated().any() == False

gdfORNode_simplified, gdfORLink_simplified = nodes_gdf_from_edges_gdf(gdfORLink_simplified, 'startNode','endNode')

# Rename fid columns and node columns to match other road network data columns
gdfORNode_simplified = gdfORNode_simplified.rename(columns = {'id':'node_fid'})
gdfORLink_simplified = gdfORLink_simplified.rename(columns = {"osmid":"old_fid", "startNode":"MNodeFID", "endNode":"PNodeFID"})

# Checking that all node ids in link data match with a node id in nodes data
assert gdfORLink_simplified.loc[ ~(gdfORLink_simplified['MNodeFID'].isin(gdfORNode_simplified['node_fid']))].shape[0] == 0
assert gdfORLink_simplified.loc[ ~(gdfORLink_simplified['PNodeFID'].isin(gdfORNode_simplified['node_fid']))].shape[0] == 0

assert gdfORLink_simplified['fid'].duplicated().any() == False
assert gdfORNode_simplified['node_fid'].duplicated().any() == False


gdfORLink_simplified.crs = projectCRS
gdfORNode_simplified.crs = projectCRS


################################
#
# alternatively use osmnx to select open road network
#
################################

# This point based mathod doesn't worth for some reason
'''
gdfPOIsWSG84 = gdfPOIs.to_crs(epsg=4326)
centre_point = gdfPOIsWSG84.loc[gdfPOIsWSG84['ref_no'] == config['centre_poi_ref'], 'geometry'].values[0].coords[0]
open_network = osmnx.graph.graph_from_point(centre_point, dist=2500, dist_type='bbox', network_type='all', simplify=True, retain_all=False, truncate_by_edge=True, clean_periphery=True, custom_filter=None)
'''

# Try using polygon instead
open_network = osmnx.graph.graph_from_polygon(studyPolygonWSG84, network_type='all', simplify=True, retain_all=False, truncate_by_edge=True, clean_periphery=True, custom_filter=None)

# Get undirected non multi graph version
D = osmnx.get_digraph(open_network) # Converts from multi di graph to di graph
U = D.to_undirected()

gdf_nodes, gdf_edges = graph_to_gdfs(U)
gdf_edges.reset_index(inplace = True)
gdf_nodes = gdf_nodes.to_crs(projectCRS)
gdf_edges = gdf_edges.to_crs(projectCRS)

# Find node closest to centre OR node
gdf_nodes['dist_to_centre'] = gdf_nodes.distance(centre_poi_geom)
nearest_node_id = gdf_nodes.sort_values(by = 'dist_to_centre', ascending=True).index[0]

reachable_nodes = largest_connected_component_nodes_within_dist(U, nearest_node_id, 2500, 'length')

gdf_nodes = gdf_nodes.loc[reachable_nodes]
gdf_edges = gdf_edges.loc[ ( gdf_edges['u'].isin(reachable_nodes)) & (gdf_edges['v'].isin(reachable_nodes))]

gdf_nodes.reset_index(inplace=True)

# osmid col contains list, need to convert to single string
for col in gdf_edges.columns:
    gdf_edges.loc[gdf_edges[col].map(lambda v: isinstance(v, list)), col] = gdf_edges.loc[gdf_edges[col].map(lambda v: isinstance(v, list)), col].map(lambda v: "-".join(str(i) for i in v))

#############################
#
#
# Save the processed data
#
#
#############################

gdfITNLink.to_file(output_itn_link_file)
gdfITNNode.to_file(output_itn_node_file)

gdfORLink_simplified.to_file(output_or_link_file)
gdfORNode_simplified.to_file(output_or_node_file)

gdf_edges.to_file(os.path.join(output_directory, "osmnx_edges.shp"))
gdf_nodes.to_file(os.path.join(output_directory, "osmnx_nodes.shp"))