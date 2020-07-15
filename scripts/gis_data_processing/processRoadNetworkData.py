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
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString

######################
#
#
# Functions
#
#
######################
import math
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
        angle+=ang(la,lb)

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

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"

itn_directory = os.path.join(gis_data_dir, "mastermap-itn_2903030")
itn_link_file = os.path.join(itn_directory, "mastermap-itn RoadLink", "mastermap-itn 2903030_0 RoadLink.shp")
itn_node_file = os.path.join(itn_directory, "mastermap-itn RoadNode", "mastermap-itn RoadNode.shp")

open_roads_directory = os.path.join(gis_data_dir, "OS Open Roads\\open-roads_subsection")
or_link_file = os.path.join(open_roads_directory, "RoadLink.shp")
or_node_file = os.path.join(open_roads_directory, "RoadNode.shp")

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

# Study area polygon - to select data within the study area
gdfSelect = gpd.read_file(selection_layer_file)
gdfSelect.crs = projectCRS

# Get the area to filter geometries by
SelectPolygon = gdfSelect.loc[0,'geometry']



################################
#
# Select the ITN Road Network that lies in the study area
#
################################

# Select only the polygons that intersect or lie within the junc clip area
gdfITNLink = gdfITNLink.loc[ (gdfITNLink.geometry.intersects(SelectPolygon)) | (gdfITNLink.geometry.within(SelectPolygon))]
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

# Rename identifier to fid to match the id col used in the ITN
gdfORLink = gdfORLink.rename(columns = {"identifier": "fid"})
gdfORNode = gdfORNode.rename(columns = {"identifier": "fid"})


gdfORLink = gdfORLink.loc[ (gdfORLink.geometry.intersects(SelectPolygon)) | (gdfORLink.geometry.within(SelectPolygon))]
gdfORNode = gpd.sjoin(gdfORNode, gdfORLink.loc[:,['fid','geometry']], op = 'intersects', lsuffix = 'node', rsuffix = 'line')

# Clean up
gdfORNode.drop(['fid_line', 'index_line'], axis = 1, inplace=True)
gdfORNode.drop_duplicates(inplace = True)
gdfORNode.rename(columns = {'fid_node':'fid'}, inplace = True)

# Clean data to ensure minimum angular deviation along road link

assert gdfORLink['geometry'].type.unique().size == 1
assert gdfORLink['fid'].unique().size == gdfORLink.shape[0]

'''
# Load the node gis data
gdfORNode = gpd.read_file(os.path.join(gis_data_dir, node_shapefile))
gdfORNode.crs = projectCRS
gdfORNode = gdfORNode.reindex(columns=['identifier','geometry'])

# Merge the nodes with the links
gdfORLink = gdfORLink.merge(gdfORNode, left_on = 'startNode', right_on = 'identifier', how = 'left', suffixes = ('','_start_node'), indicator= True)
print(gdfORLink['_merge'].value_counts())
assert gdfORLink.loc[ gdfORLink['_merge'] != 'both'].shape[0] == 0 #- not all nodes are in the road link clipped data due to the way it was clipped. should change this
gdfORLink.rename(columns={'_merge':'_merge_plus_node'}, inplace=True)


gdfORLink = gdfORLink.merge(gdfORNode, left_on = 'endNode', right_on = 'identifier', how = 'left', suffixes = ('', '_end_node'), indicator=True)
print(gdfORLink['_merge'].value_counts())
assert gdfORLink.loc[ gdfORLink['_merge'] != 'both'].shape[0] == 0 #- - not all nodes are in the road link clipped data due to the way it was clipped. should change this
gdfORLink.rename(columns={'_merge':'_merge_minus_node'}, inplace=True)
'''


###################################
#
#
# # Simplify linestrings
# identify new nodes that will result from splitting lines
# collect these in a df
# split linestrings re-ID components
# find which linestrings' start and end coords dont match the node coordinate
# add these coords as new nodes
#
####################################
gdfORLink_simplified = simplify_line_gdf_by_angle(gdfORLink, 10, "fid", "new_fid")

assert gdfORLink_simplified['new_fid'].duplicated().any() == False

gdfORLink_simplified['startCoord'] = gdfORLink_simplified['geometry'].map(lambda x: Point(x.coords[0]))
gdfORLink_simplified['endCoord'] = gdfORLink_simplified['geometry'].map(lambda x: Point(x.coords[-1]))

# Collect all start and end nodes
coords = pd.concat([gdfORLink_simplified['startCoord'], gdfORLink_simplified['endCoord']]).drop_duplicates()
gdfORNode_simplified = gpd.GeoDataFrame({'geometry':coords})
gdfORNode_simplified['node_fid'] = ["node_id_{}".format(i) for i in gdfORNode_simplified.index]


# Drop the old node identifierd and merge to nodes df to get new identifiers/ids
gdfORLink_simplified = gdfORLink_simplified.drop(['startNode','endNode'], axis = 1)

# Merge the nodes with the links
gdfORLink_simplified = gdfORLink_simplified.set_geometry("startCoord")
gdfORLink_simplified = gpd.geopandas.sjoin(gdfORLink_simplified, gdfORNode_simplified, how='inner', op='intersects', lsuffix='left', rsuffix='right')
assert gdfORLink_simplified['node_fid'].isnull().any() == False
gdfORLink_simplified.rename(columns={'node_fid':'startNode'}, inplace=True)
gdfORLink_simplified = gdfORLink_simplified.drop(['index_right'], axis = 1)


gdfORLink_simplified = gdfORLink_simplified.set_geometry("endCoord")
gdfORLink_simplified = gpd.geopandas.sjoin(gdfORLink_simplified, gdfORNode_simplified, how='inner', op='intersects', lsuffix='left', rsuffix='right')
assert gdfORLink_simplified['node_fid'].isnull().any() == False 
gdfORLink_simplified.rename(columns={'node_fid':'endNode'}, inplace=True)
gdfORLink_simplified = gdfORLink_simplified.drop(['index_right'], axis = 1)

gdfORLink_simplified = gdfORLink_simplified.set_geometry("geometry")
# Clean up and save data
gdfORLink_simplified = gdfORLink_simplified.drop(["startCoord", "endCoord"], axis = 1)

# Rename fid columns and node columns to match other road network data columns
gdfORLink_simplified = gdfORLink_simplified.rename(columns = {"fid":"old_fid", "new_fid":"fid", "startNode":"MNodeFID", "endNode":"PNodeFID"})
assert gdfORLink_simplified['fid'].duplicated().any() == False


gdfORLink_simplified.crs = projectCRS
gdfORNode_simplified.crs = projectCRS


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