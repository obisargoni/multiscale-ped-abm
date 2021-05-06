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
from shapely.geometry import Point, MultiPoint, Polygon, MultiPolygon, LineString, MultiLineString
from shapely import ops
import itertools
import copy

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

def _is_point_in_nodes(p, nodes):
    # Search new nodes for this coordinate
    point_in_node = False
    for v in nodes.values():
        if (p.x == v['x']) & (p.y == v['y']):
            # This coordiante has already been recorded as a new node
            point_in_node = False
            break
    return point_in_node

def _node_data_from_point(p):
    d = {}
    d['x'] = p.x
    d['y'] = p.y
    d['geometry'] = p
    return d

def _match_nodes_to_geometry(g, u, v, graph_nodes, other_nodes):
    """Find start and end nodes for the new edge geometry. 

    First check if input u and v nodes match start and end point of geometry. Then search through input list of nodes.
    Allows distingusing between nodes that already are in a graph (u and v) and node that need to be added.
    """

    new_u = None
    new_u_data = None
    if (g.coords[0][0] == graph_nodes[u]['x']) & (g.coords[0][1] == graph_nodes[u]['y']):
        new_u = u
        new_u_data = graph_nodes[new_u]
    elif (g.coords[0][0] == graph_nodes[v]['x']) & (g.coords[0][1] == graph_nodes[v]['y']):
        new_u = v
        new_u_data = graph_nodes[new_u]
    else:
        # Find which newly created nodes matched the end of this line string
        for node_id, node_data in other_nodes.items():
            if (g.coords[0][0] == node_data['x']) & (g.coords[0][1] == node_data['y']):
                new_u = node_id
                new_u_data = node_data
                break
    
    new_v = None
    new_v_data = None
    if (g.coords[-1][0] == graph_nodes[u]['x']) & (g.coords[-1][1] == graph_nodes[u]['y']):
        new_v = u
        new_v_data = graph_nodes[new_v]
    elif (g.coords[-1][0] == graph_nodes[v]['x']) & (g.coords[-1][1] == graph_nodes[v]['y']):
        new_v = v
        new_v_data = graph_nodes[new_v]
    else:
        # Find which newly created nodes matched the end of this line string
        for node_id, node_data in other_nodes.items():
            if (g.coords[-1][0] == node_data['x']) & (g.coords[-1][1] == node_data['y']):
                new_v = node_id
                new_v_data = node_data
                break

    return (new_u, new_u_data), (new_v, new_v_data)
                
def largest_connected_component_nodes_within_dist(G, source_node, dist, weight):
    if G.is_directed():
        lccNodes = max(nx.weakly_connected_components(G), key=len)
    else:
        lccNodes = max(nx.connected_components(G), key=len)

    lccG = G.subgraph(lccNodes).copy()
    
    shortest_paths = nx.single_source_dijkstra(lccG, source_node, target=None, cutoff=dist, weight=weight)   
    reachable_nodes = shortest_paths[0].keys()

    return reachable_nodes


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

    gdf_edges['c1'] = gdf_edges['geometry'].map(lambda g: Point(g.coords[0]))
    gdf_edges['c2'] = gdf_edges['geometry'].map(lambda g: Point(g.coords[1]))

    node_coords = pd.concat([gdf_edges['c1'], gdf_edges['c2']]).drop_duplicates()
    node_ids = ['or_node_{}'.format(i) for i in np.arange(len(node_coords))]
    gdf_nodes = gpd.GeoDataFrame({'node_id': node_ids, 'geometry':node_coords})
    gdf_nodes.crs = gdf_edges.crs

    # Join nodes to edges on coordinate
    gdf_edges = gdf_edges.set_geometry("c1")
    gdf_edges = gpd.geopandas.sjoin(gdf_edges, gdf_nodes, how='inner', op='intersects', lsuffix='left', rsuffix='right')
    assert gdf_edges['node_id'].isnull().any() == False
    gdf_edges.rename(columns={'node_id':u}, inplace=True)
    gdf_edges = gdf_edges.drop(['index_right'], axis = 1)

    gdf_edges = gdf_edges.set_geometry("c2")
    gdf_edges = gpd.geopandas.sjoin(gdf_edges, gdf_nodes    , how='inner', op='intersects', lsuffix='left', rsuffix='right')
    assert gdf_edges['node_id'].isnull().any() == False 
    gdf_edges.rename(columns={'node_id':v}, inplace=True)
    gdf_edges = gdf_edges.drop(['index_right'], axis = 1)

    # Tidy up
    gdf_edges = gdf_edges.set_geometry("geometry")
    gdf_edges = gdf_edges.drop(["c1", "c2"], axis = 1)

    return gdf_nodes, gdf_edges

def simplify_graph(G, strict=True, remove_rings=True, rebuild_geoms = False):
    """
    Simplify a graph's topology by removing interstitial nodes.
    Simplifies graph topology by removing all nodes that are not intersections
    or dead-ends. Create an edge directly between the end points that
    encapsulate them, but retain the geometry of the original edges, saved as
    a new `geometry` attribute on the new edge. Note that only simplified
    edges receive a `geometry` attribute. Some of the resulting consolidated
    edges may comprise multiple OSM ways, and if so, their multiple attribute
    values are stored as a list.
    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    strict : bool
        if False, allow nodes to be end points even if they fail all other
        rules but have incident edges with different OSM IDs. Lets you keep
        nodes at elbow two-way intersections, but sometimes individual blocks
        have multiple OSM IDs within them too.
    remove_rings : bool
        if True, remove isolated self-contained rings that have no endpoints
    Returns
    -------
    G : networkx.MultiDiGraph
        topologically simplified graph, with a new `geometry` attribute on
        each simplified edge
    """
    if "simplified" in G.graph and G.graph["simplified"]:  # pragma: no cover
        raise Exception("This graph has already been simplified, cannot simplify it again.")

    osmnx.utils.log("Begin topologically simplifying the graph...")

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    initial_node_count = len(G)
    initial_edge_count = len(G.edges)
    all_nodes_to_remove = []
    all_edges_to_add = []

    # generate each path that needs to be simplified
    for path in osmnx.simplification._get_paths_to_simplify(G, strict=strict):

        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        edge_attributes = dict()
        for u, v in zip(path[:-1], path[1:]):

            # there should rarely be multiple edges between interstitial nodes
            # usually happens if OSM has duplicate ways digitized for just one
            # street... we will keep only one of the edges (see below)
            if G.number_of_edges(u, v) != 1:
                osmnx.utils.log(f"Found multiple edges between {u} and {v} when simplifying")

            # get edge between these nodes: if multiple edges exist between
            # them (see above), we retain only one in the simplified graph
            edge = G.edges[u, v, 0]
            for key in edge:
                if key in edge_attributes:
                    # if this key already exists in the dict, append it to the
                    # value list
                    if isinstance(edge[key], (list, tuple)):
                        for item in edge[key]:
                            edge_attributes[key].append(item)
                    else:
                        edge_attributes[key].append(edge[key])
                else:
                    # if this key doesn't already exist, but the value is a list set the value to the edge value
                    # otherwise set the value to a list containing the one value
                    if isinstance(edge[key], (list, tuple)):
                        edge_attributes[key] = list(edge[key])
                    else:
                        edge_attributes[key] = [edge[key]]

        for key in edge_attributes:

            # don't touch the length or geometry attribute, we'll sum it at the end
            if key in ["length", "geometry"]:
                continue
            elif len(set(edge_attributes[key])) == 1:
                # if there's only 1 unique value in this attribute list,
                # consolidate it to the single value (the zero-th)
                edge_attributes[key] = edge_attributes[key][0]
            else:
                # otherwise, if there are multiple values, keep one of each value
                edge_attributes[key] = list(set(edge_attributes[key]))

        # construct the geometry and sum the lengths of the segments
        if rebuild_geoms:
            edge_attributes["geometry"] = LineString(
                [Point((G.nodes[node]["x"], G.nodes[node]["y"])) for node in path]
            )
            edge_attributes["length"] = sum(edge_attributes["length"])
        else:
            # Create single geometry from the coordinates of the component geometries
            merged_line = ops.linemerge(MultiLineString(edge_attributes["geometry"]))
            edge_attributes["geometry"] = merged_line
            edge_attributes["length"] = merged_line.length

        # add the nodes and edges to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append(
            {"origin": path[0], "destination": path[-1], "attr_dict": edge_attributes}
        )

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge["origin"], edge["destination"], **edge["attr_dict"])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    if remove_rings:
        # remove any connected components that form a self-contained ring
        # without any endpoints
        wccs = nx.weakly_connected_components(G)
        nodes_in_rings = set()
        for wcc in wccs:
            if not any(osmnx.simplification._is_endpoint(G, n) for n in wcc):
                nodes_in_rings.update(wcc)
        G.remove_nodes_from(nodes_in_rings)

    # mark graph as having been simplified
    G.graph["simplified"] = True

    msg = (
        f"Simplified graph: {initial_node_count} to {len(G)} nodes, "
        f"{initial_edge_count} to {len(G.edges)} edges"
    )
    osmnx.utils.log(msg)
    return G

def break_overlapping_edges(G, id_attr = 'strg_id'):
    """Used to deal with edges geometries that overlap. Create nodes at intersection between edges
    and break edge geometries at these points. Rebuild graph by assigning nodes to start and end of graph.

    --------------
    Returns:
        geopandas.GeoDataFrame, geopandas.GeoDataFrame

        gdfNodes, gdfLinks
    """
    print("Breaking overlapping edges")
    H = nx.MultiGraph()
    H.graph = G.graph

    new_nodes = {} # Record the points of intersections between edge geometries and use these as addition nodes
    nodes = G.nodes(data=True)
    new_node_index = 0
    for n in nodes:
        u = n[0]
        # Now loop through pairs of edges attached to this node
        for e1, e2 in itertools.combinations(G.edges(u,data=True, keys=True), 2):

            # Unpack for ease
            u1,v1,k1,d1 = e1
            u2,v2,k2,d2 = e2

            g1 = copy.deepcopy(d1['geometry'])
            g2 = copy.deepcopy(d2['geometry'])

            # Remove overlaps between geometries
            g1 = g1.difference(g2)

            if (g1.is_empty) & ~(g2.is_empty):
                print("Empty line. e1:{}, e2:{}".format((u1,v1), (u2,v2)))
                continue
            elif ~(g1.is_empty) & (g2.is_empty):
                print("Empty line. e1:{}, e2:{}".format((u1,v1), (u2,v2)))
                continue

            # Find intersecting points of geometries
            intersection = g1.intersection(g2)

            if isinstance(intersection, Point):
                #node_data['geometry'].append(intersection)
                intersection = MultiPoint([intersection])
            elif isinstance(intersection, MultiPoint):
                pass
            else:
                print("Unexpected intersection. e1:{}, e2:{}".format((u1,v1,k1), (u2,v2,k2)))

            # Add intersection points to dict of new nodes
            for p in intersection:

                point_in_nodes = _is_point_in_nodes(p, new_nodes)
                
                if point_in_nodes:
                    continue
                else:
                    node_id = "intersection_node_{}".format(new_node_index)
                    new_node_index += 1
                    new_nodes[node_id] = _node_data_from_point(p)

    # Create multipoint geometry containing all the nodes, use this to split edges
    points = [Point(n[-1]['x'],n[-1]['y']) for n in nodes]
    points += [v['geometry'] for k,v in new_nodes.items()]
    mp = MultiPoint(points)

    # Split edge geometries by the intersecting points
    for e in G.edges(data=True, keys=True):
        data = copy.deepcopy(e[-1])
        g = data['geometry']
        new_edge_index = 0
        for new_geom in ops.split(g, mp):
            # Find start and end nodes for the new edge geometry
            # Although start and end coords might match points of intersection, retain original node IDs where possible
            # Means that only intersection node sthat are used should be added to the graph, otherwise will get duplicated node geometries
            (new_u, new_u_data), (new_v, new_v_data) = _match_nodes_to_geometry(new_geom, e[0], e[1], nodes, new_nodes)

            # Finally add the edge to the new graph
            if (new_u is None) | (new_v is None):
                print("New nodes not found for edge {}".format(e))
            else:
                H.add_node(new_v, **new_v_data) # If node has already been added the graph will be unchanged
                H.add_node(new_u, **new_u_data) # If node has already been added the graph will be unchanged
                
                new_data = data
                new_data['geometry'] = new_geom
                new_data['sub_'+id_attr] = "sub_{}".format(new_edge_index)
                new_edge_index += 1
                H.add_edge(new_u, new_v, **new_data)

    return H

def remove_duplicated_edges(G):
    """Breaking up edges where they overlap can create multiple links between the same two nodes with the same geometry.
    These are considered duplicate link. This functions removes these links from the graph.
    """
    D = G.copy()

    # Loop through node pairs, if multiple keys check if geometry is duplicated and delete if so
    for u in D.nodes():

        if isinstance(G, nx.MultiGraph):
            neighbours = set(D.neighbors(u))
        elif isinstance(G, nx.MultiDiGraph):
            neighbours = set(list(D.predecessors(u)) + list(D.successors(u)))

        for v in neighbours:

            keys = list(D[u][v].keys())

            if len(keys)==1:
                # no multi edge between these nodes
                continue
            else:
                # Check for duplicate geometries
                geoms = [D[u][v][k]['geometry'] for k in keys]

                for k1, k2 in itertools.combinations(keys, 2):
                    try:
                        g1 = D[u][v][k1]['geometry']
                        g2 = D[u][v][k2]['geometry']
                    except KeyError as ke:
                        # Get key error if an edge has been removed
                        continue

                    if g1.equals(g2):
                        D.remove_edge(u,v,key = k2)

                # Check there is still an edge between u and v
                try:
                    e = D[u][v]
                except KeyError:
                    assert False
    return D



def duplicate_geometry_row_ids(gdf, geometry = 'geometry'):
    dup_ids = []
    for ix, row in gdf.iterrows():
        g1 = row[geometry]
        temp = [ix]
        for ix_, row_ in gdf.loc[ix:,].iterrows():
            if ix_ == ix:
                continue
            g2 = row_[geometry]
            if g1.equals(g2):
                temp.append(ix_)
        if len(temp) > 1:
            dup_ids.append(temp)
    return dup_ids

def drop_duplicate_geometries(gdf, geometry = 'geometry', **kwargs):
    duplicated_ids = duplicate_geometry_row_ids(gdf, geometry = geometry)
    for id_group in duplicated_ids:
        gdf = gdf.drop(id_group[1:], axis=0, errors = 'ignore') # Keep the first entry, ignore errors since possible that this will try to drop the same row multiple times
    return gdf


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

seriesStudyArea = centre_poi.buffer(config['study_area_dist']+500)
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
gdfORLink = gdfORLink[ gdfORLink['geometry'].type == "LineString"]

# Handle multi links
gdfORLink['key'] = None
gdfORLink['key'] = gdfORLink.groupby(['startNode','endNode'])['key'].transform(lambda df: np.arange(df.shape[0]))
assert gdfORLink.loc[:, ['startNode','endNode', 'key']].duplicated().any() == False

# Represent undirected network as directed graph
gdfORLinkReversed = gdfORLink.copy()
gdfORLinkReversed = gdfORLink.rename(columns = {'startNode':'endNode', 'endNode':'startNode'})
gdfORLinkReversed['geometry'] = gdfORLinkReversed['geometry'].map(lambda g: LineString(g.coords[::-1]))
gdfORLink = pd.concat([gdfORLink, gdfORLinkReversed])

# Format for osmnx
gdfORLink = gdfORLink.rename(columns = {"identifier": "osmid", 'startNode':'u', 'endNode':'v'})
gdfORLink.set_index(['u','v','key'], inplace=True)
gdfORLink['geometry'] = gdfORLink['geometry'].map(make_linestring_coords_2d)

gdfORNode["geometry"] = gdfORNode["geometry"].map(lambda g: g[0])
gdfORNode = gdfORNode.rename(columns = {"identifier": "osmid"})
gdfORNode['x'] = gdfORNode.loc[:, 'geometry'].map(lambda g: g.x)
gdfORNode['y'] = gdfORNode.loc[:, 'geometry'].map(lambda g: g.y)
gdfORNode.set_index('osmid', inplace=True)


# Makes sense to set up graph as osmnx compliant object. But need to make sure I can keep track of edge ids
G = osmnx.graph_from_gdfs(gdfORNode, gdfORLink, graph_attrs=None)

# Find the or node nearest the centre poi
gdfORNode['dist_to_centre'] = gdfORNode.distance(centre_poi_geom)
nearest_node_id = gdfORNode.sort_values(by = 'dist_to_centre', ascending=True).index[0]

# Get largest connected component
reachable_nodes = largest_connected_component_nodes_within_dist(G, nearest_node_id, config['study_area_dist'], 'length')

G = G.subgraph(reachable_nodes).copy()

# Remove dead ends by removing nodes with degree 1 continually  until no degree 1 nodes left
U = G.to_undirected()
dfDegree = pd.DataFrame(U.degree(), columns = ['nodeID','degree'])
dead_end_nodes = dfDegree.loc[dfDegree['degree']==1, 'nodeID'].values
removed_nodes = []
while(len(dead_end_nodes)>0):
    U.remove_nodes_from(dead_end_nodes)
    removed_nodes = np.concatenate([removed_nodes, dead_end_nodes])
    dfDegree = pd.DataFrame(U.degree(), columns = ['nodeID','degree'])
    dead_end_nodes = dfDegree.loc[dfDegree['degree']==1, 'nodeID'].values

G.remove_nodes_from(removed_nodes)

###################################
#
#
# Now that study area network has been selected, clean network by:
# - simplify intersections
# - split lines based on angular deviation
#
####################################

# simplify topology before breaking up edges based on angular deviation
G_simp = simplify_graph(G, strict=True, remove_rings=True, rebuild_geoms = False)

# simplify intersections
G_simp = osmnx.consolidate_intersections(G_simp, tolerance=15, rebuild_graph=True, dead_ends=True, reconnect_edges=True)

# At this stage reset the road link IDs. These road links correspond to strtegic routing.
# This code keeps the old attributes but not sure that's needed.
new_attributes = {}
for i, edge in enumerate(G_simp.edges(data = True, keys = True)):
    edge_id = "strategic_{}".format(i)
    edge_attributes = {}
    edge_attributes['strg_id'] = edge_id
    new_attributes[edge[:-1]] = edge_attributes
nx.set_edge_attributes(G_simp, new_attributes)

# Convert to undirected for next bit of cleaning. Keep multi edge representation though. Need to think about this - but think it makes sense to retain most general structure
U = G_simp.to_undirected()
U_clip = U.copy()
U_clip = break_overlapping_edges(U_clip)


U_clip = remove_duplicated_edges(U_clip)

U_clip = U_clip.to_directed()
U_clip.graph['simplified'] = False
U_clip = simplify_graph(U_clip, strict=True, remove_rings=True, rebuild_geoms = False)
U_clip = U_clip.to_undirected()


gdfNodes, gdfLinks = osmnx.graph_to_gdfs(U_clip)

# Reset indexes and convert ids to string data
gdfNodes.reset_index(inplace=True)
gdfNodes['osmid'] = gdfNodes['osmid'].map(lambda x: str(x))

gdfLinks.reset_index(inplace=True)
gdfLinks['u'] = gdfLinks['u'].map(lambda x: str(x))
gdfLinks['v'] = gdfLinks['v'].map(lambda x: str(x))
gdfLinks['key'] = gdfLinks['key'].map(lambda x: str(x))

gdfLinks = gdfLinks.reindex(columns = ['u','v','key','osmid','u_original','v_original', 'strg_id', 'sub_strg_id', 'geometry'])
gdfLinks['link_id'] = gdfLinks['strg_id'].map(str) + "_" + gdfLinks['sub_strg_id'].map(str)


for col in gdfLinks.columns:
    gdfLinks.loc[gdfLinks[col].map(lambda v: isinstance(v, list)), col] = gdfLinks.loc[gdfLinks[col].map(lambda v: isinstance(v, list)), col].map(lambda v: "_".join(str(i) for i in v))

for col in gdfNodes.columns:
    gdfNodes.loc[gdfNodes[col].map(lambda v: isinstance(v, list)), col] = gdfNodes.loc[gdfNodes[col].map(lambda v: isinstance(v, list)), col].map(lambda v: "_".join(str(i) for i in v))

assert gdfLinks.loc[:, ['u','v','key']].duplicated().any() == False
assert gdfLinks.loc[:, ['link_id']].duplicated().any() == False

gdfNodes.crs = projectCRS
gdfLinks.crs = projectCRS
gdfNodes.to_file(os.path.join("or_roads", "or_node_cleaned_split.shp"))
gdfLinks.to_file(os.path.join("or_roads", "or_link_cleaned_split.shp"))


gdfORLink_simplified = simplify_line_gdf_by_angle(gdfLinks, 10, "link_id", "fid")

assert gdfORLink_simplified['fid'].duplicated().any() == False

gdfORNode_simplified, gdfORLink_simplified = nodes_gdf_from_edges_gdf(gdfORLink_simplified, 'startNode','endNode')

# Rename fid columns and node columns to match other road network data columns
gdfORNode_simplified = gdfORNode_simplified.rename(columns = {'node_id':'node_fid'})
gdfORLink_simplified = gdfORLink_simplified.rename(columns = {"osmid":"old_fid", "startNode":"MNodeFID", "endNode":"PNodeFID"})


# At this stage can have some duplicated geometries. Check for this and delete duplications
gdfORLink_simplified.index = np.arange(gdfORLink_simplified.shape[0])
gdfORLink_simplified = drop_duplicate_geometries(gdfORLink_simplified)

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
open_network = osmnx.graph.graph_from_point(centre_point, dist=config['study_area_dist'], dist_type='bbox', network_type='all', simplify=True, retain_all=False, truncate_by_edge=True, clean_periphery=True, custom_filter=None)
'''

# Try using polygon instead
open_network = osmnx.graph.graph_from_polygon(studyPolygonWSG84, network_type='all', simplify=True, retain_all=False, truncate_by_edge=True, clean_periphery=True, custom_filter=None)

# Get undirected non multi graph version
#D = osmnx.get_digraph(open_network) # Converts from multi di graph to di graph
U = open_network.to_undirected()

gdf_nodes, gdf_edges = osmnx.graph_to_gdfs(U)
gdf_edges.reset_index(inplace = True)
gdf_nodes = gdf_nodes.to_crs(projectCRS)
gdf_edges = gdf_edges.to_crs(projectCRS)

# Find node closest to centre OR node
gdf_nodes['dist_to_centre'] = gdf_nodes.distance(centre_poi_geom)
nearest_node_id = gdf_nodes.sort_values(by = 'dist_to_centre', ascending=True).index[0]

reachable_nodes = largest_connected_component_nodes_within_dist(U, nearest_node_id, config['study_area_dist'], 'length')

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