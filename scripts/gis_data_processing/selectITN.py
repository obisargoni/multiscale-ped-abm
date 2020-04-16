import os
import geopandas as gpd
import pandas as pd
#from operator import itemgetter

clip_layer_file = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\JunctClip.shp"

data_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\mastermap-itn_2903030\\mastermap-itn RoadLink"
fileITN = 'mastermap-itn 2903030_0 RoadLink.shp'

node_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\mastermap-itn_2903030\\mastermap-itn RoadNode"
node_file = "mastermap-itn RoadNode.shp"

output_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\Intersect Within layers"
itn_file = "mastermap-itn RoadLink Intersect Within.shp"
out_node_file = "mastermap-itn RoadNode Intersect Within.shp"

# Read in the data
gdfITN = gpd.read_file(os.path.join(data_directory,fileITN))
gdfITNNode = gpd.read_file(os.path.join(node_directory, node_file))
gdfClip = gpd.read_file(clip_layer_file)

# Set crs
gdfITN.crs = {'init' :'epsg:27700'}
gdfClip.crs = {'init' :'epsg:27700'}
gdfITNNode.crs = {'init' :'epsg:27700'}

ClipArea = gdfClip.loc[0,'geometry']

# Select only the polygons that intersect or lie within the junc clip area
gdfITN = gdfITN.loc[ (gdfITN.geometry.intersects(ClipArea)) | (gdfITN.geometry.within(ClipArea))]

gdfITNNode = gpd.sjoin(gdfITNNode, gdfITN.loc[:,['fid','geometry']], op = 'intersects', lsuffix = 'node', rsuffix = 'line')

# Clean up
gdfITNNode.drop('fid_line', axis = 1, inplace=True)
gdfITNNode.rename(columns = {'fid_node':'fid'}, inplace = True)

# Could do some cleaning here - multipart to singlepart

# Save the selected ITN geometries
gdfITN.to_file(os.path.join(output_directory, itn_file))
gdfITNNode.to_file(os.path.join(output_directory, out_node_file))