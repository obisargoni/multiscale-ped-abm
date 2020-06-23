'''
This is a script for editing OS's ITN Road Network shapefiles so that:
- the attributeds include an id for the 'to' and 'from' nodes
- line strings are duplicated along links that are bidirectional

Purpose of this script is to use the extracted orientation information from the gml data to edit the roads linesting shape file such that 
two way road are duplicated with one linestring orientated in one direction and the other orientated in another direction.
'''


import numpy as np
import pandas as pd
import geopandas as gpd
import os
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString


######################
#
#
# Initialise variables and paths to data inputs and outputs
#
#
#####################
project_crs = "epsg:27700"

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\processed_gis_data"

link_shapefile = "open-roads RoadLink Intersect Within.shp"
node_shapefile = "open-roads RoadNode Intersect Within.shp"

output_dir = gis_data_dir
output_link_shapefile = os.path.join(output_dir, "open_roads RoadLink Intersect Within simplify angles.shp")
output_node_shapefile = os.path.join(output_dir, "open_roads RoadNode Intersect Within simplify angles.shp")


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


###########################
#
#
# Load Data
# 
# merge with nodes on ids
#
##########################

# Read in data. Check geometry types. Check fields
gdfORLink = gpd.read_file(os.path.join(gis_data_dir, link_shapefile))
gdfORLink.crs = project_crs

assert gdfORLink['geometry'].type.unique().size == 1
assert gdfORLink['identifier'].unique().size == gdfORLink.shape[0]

# Check that each link has id of plus and minus nodes
assert gdfORLink['startNode'].isnull().any() == False
assert gdfORLink['endNode'].isnull().any() == False

'''
# Load the node gis data
gdfORNode = gpd.read_file(os.path.join(gis_data_dir, node_shapefile))
gdfORNode.crs = project_crs
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
gdfORLink_simplified = simplify_line_gdf_by_angle(gdfORLink, 10, "identifier", "new_identifier")

assert gdfORLink_simplified['new_identifier'].duplicated().any() == False 

gdfORLink_simplified['startCoord'] = gdfORLink_simplified['geometry'].map(lambda x: Point(x.coords[0]))
gdfORLink_simplified['endCoord'] = gdfORLink_simplified['geometry'].map(lambda x: Point(x.coords[-1]))

# Collect all start and end nodes
coords = pd.concat([gdfORLink_simplified['startCoord'], gdfORLink_simplified['endCoord']])
gdfORNode_simplified = gpd.GeoDataFrame({'geometry':coords})
gdfORNode_simplified['node_identifier'] = ["node_id_{}".format(i) for i in gdfORNode_simplified.index]


# Drop the old node identifierd and merge to nodes df to get new identifiers/ids
gdfORLink_simplified = gdfORLink_simplified.drop(['startNode','endNode'], axis = 1)

# Merge the nodes with the links
gdfORLink_simplified = gdfORLink_simplified.set_geometry("startCoord")
gdfORLink_simplified = gpd.geopandas.sjoin(gdfORLink_simplified, gdfORNode_simplified, how='inner', op='intersects', lsuffix='left', rsuffix='right')
assert gdfORLink_simplified['node_identifier'].isnull().any() == False
gdfORLink_simplified.rename(columns={'node_identifier':'startNode'}, inplace=True)
gdfORLink_simplified = gdfORLink_simplified.drop(['index_right'], axis = 1)


gdfORLink_simplified = gdfORLink_simplified.set_geometry("endCoord")
gdfORLink_simplified = gpd.geopandas.sjoin(gdfORLink_simplified, gdfORNode_simplified, how='inner', op='intersects', lsuffix='left', rsuffix='right')
assert gdfORLink_simplified['node_identifier'].isnull().any() == False 
gdfORLink_simplified.rename(columns={'node_identifier':'endNode'}, inplace=True)
gdfORLink_simplified = gdfORLink_simplified.drop(['index_right'], axis = 1)

gdfORLink_simplified = gdfORLink_simplified.set_geometry("geometry")
# Clean up and save data
gdfORLink_simplified = gdfORLink_simplified.drop(["startCoord", "endCoord"], axis = 1)

gdfORLink_simplified.crs = project_crs
gdfORNode_simplified.crs = project_crs

gdfORLink_simplified.to_file(output_link_shapefile)
gdfORNode_simplified.to_file(output_node_shapefile)