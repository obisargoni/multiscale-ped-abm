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


######################
#
#
# Initialise variables and paths to data inputs and outputs
#
#
#####################
crs = "epsg:27700"

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\processed_gis_data"

link_shapefile = "mastermap-itn RoadLink Intersect Within.shp"
node_shapefile = "mastermap-itn RoadNode Intersect Within.shp"

route_info_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\itn_route_info"
road_route_info_path = os.path.join(route_info_dir, "extracted_RRI.csv")
road_node_info_path = os.path.join(route_info_dir, "extracted_RLNodes.csv")

output_dir = gis_data_dir
output_shapefile = "mastermap-itn RoadLink Intersect Within with orientation.shp"


###########################
#
#
# Load Data
#
#
##########################

# Read in data. Check geometry types. Check fields
gdf_itn = gpd.read_file(os.path.join(gis_data_dir, link_shapefile))

assert gdf_itn['geometry'].type.unique().size == 1
assert gdf_itn['fid'].unique().size == gdf_itn.shape[0]

# Read in node information
dfRLNode = pd.read_csv(road_node_info_path)

# Join the itn road links with the node info to check that the -,+ nodes correspond to the first/last nodes of the linestring
gdf_itn = pd.merge(gdf_itn, dfRLNode, left_on="fid", right_on="RoadLinkFID", how = "left", indicator=True)

assert gdf_itn.loc[ gdf_itn['_merge'] != 'both'].shape[0] == 0
gdf_itn.drop('_merge', axis=1, inplace=True)

# Load the node gis data
gdf_itn_node = gpd.read_file(os.path.join(gis_data_dir, node_shapefile))
gdf_itn_node = gdf_itn_node.reindex(columns=['fid','geometry'])

# Merge the nodes with the links
gdf_itn = gdf_itn.merge(gdf_itn_node, left_on = 'PlusNodeFID', right_on = 'fid', how = 'left', suffixes = ('','_plus_node'), indicator= True)
print(gdf_itn['_merge'].value_counts())
assert gdf_itn.loc[ gdf_itn['_merge'] != 'both'].shape[0] == 0 #- not all nodes are in the road link clipped data due to the way it was clipped. should change this
gdf_itn.rename(columns={'_merge':'_merge_plus_node'}, inplace=True)


gdf_itn = gdf_itn.merge(gdf_itn_node, left_on = 'MinusNodeFID', right_on = 'fid', how = 'left', suffixes = ('', '_minus_node'), indicator=True)
print(gdf_itn['_merge'].value_counts())
assert gdf_itn.loc[ gdf_itn['_merge'] != 'both'].shape[0] == 0 #- - not all nodes are in the road link clipped data due to the way it was clipped. should change this
gdf_itn.rename(columns={'_merge':'_merge_minus_node'}, inplace=True)

# Drop rows where plus or minus coords have not been merged in - not required since all nodes merged
#gdf_itn = gdf_itn.loc[(gdf_itn['_merge_minus_node'] == 'both') & (gdf_itn['_merge_plus_node']=='both')]

gdf_itn['line_first_coord'] = gdf_itn['geometry'].map(lambda x: x.coords[0])
gdf_itn['line_last_coord'] = gdf_itn['geometry'].map(lambda x: x.coords[-1])

# Check that the -,+ nodes match the first / last line string coords
gdf_itn['first_coords_match_minus_node'] = gdf_itn['line_first_coord'] == gdf_itn['geometry_minus_node'].map(lambda x: x.coords[0])
gdf_itn['last_coords_match_plus_node'] = gdf_itn['line_last_coord'] == gdf_itn['geometry_plus_node'].map(lambda x: x.coords[0])

assert gdf_itn['first_coords_match_minus_node'].all() == True
gdf_itn['last_coords_match_plus_node'].all() == True

gdf_itn = None
gdf_itn_node = None


#######################################################################################


# With this check done, can now add the orientation to the itn road link data
gdf_itn_link = gpd.read_file(os.path.join(gis_data_dir, link_shapefile))

dfRRI = pd.read_csv(road_route_info_path)
gdf_itn_link = pd.merge(gdf_itn_link, dfRRI, left_on='fid', right_on='DirectedLinkFID',how = 'left', indicator=True)
gdf_itn_link = pd.merge(gdf_itn_link, dfRLNode, left_on='fid', right_on='RoadLinkFID',how = 'left', indicator=False)
#assert gdf_itn_link.loc[ gdf_itn_link['_merge'] != 'both'].shape[0] ==0 # This fails, ok cos not all road links have orientation routing info


gdf_itn_link['direction'] = gdf_itn_link['DirectedLinkOrientation']
gdf_itn_link['OneWay'] = np.nan
gdf_itn_link.loc[ gdf_itn_link['_merge'] == 'both', 'OneWay'] = True
gdf_itn_link.loc[ gdf_itn_link['_merge'] == 'left_only', 'OneWay'] = False
gdf_itn_link.drop("_merge", axis=1, inplace=True)

# Rename the nodes fields to be less than 10 chars
gdf_itn_link.rename(columns = {"PlusNodeFID":"PNodeFID", "MinusNodeFID":"MNodeFID"}, inplace = True)

gdf_itn_link.loc[ gdf_itn_link['OneWay'] == False, 'direction'] = '-'
gdf_itn_link_two_way = gdf_itn_link.loc[ gdf_itn_link['OneWay'] == False]
gdf_itn_link_two_way['direction'] = '+'
gdf_itn_link = pd.concat([gdf_itn_link, gdf_itn_link_two_way], join = 'inner')

print(gdf_itn_link['OneWay'].value_counts())
print(gdf_itn_link['direction'].value_counts())
print(gdf_itn_link['DirectedLinkOrientation'].value_counts())


gdf_itn_link.to_file(os.path.join(output_dir, output_shapefile))
