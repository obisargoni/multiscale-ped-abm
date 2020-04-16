# IPython log file

import os
import geopandas as gpd
import pandas as pd
#from operator import itemgetter

clip_layer_file = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\JunctClip.shp"

data_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\mastermap-topo_2903032\\mastermap-topo_2903032_0 TopographicArea"
fileTopo = 'mastermap TopographicArea.shp'

output_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\Intersect Within layers"
vehicle_file = "topographicAreaVehicle.shp"
pedestrian_file = "topographicAreaPedestrian.shp"

# Read in the data
gdfTopoArea = gpd.read_file(os.path.join(data_directory,fileTopo))
gdfClip = gpd.read_file(clip_layer_file)

# Set the crs
gdfTopoArea.crs = {'init' :'epsg:27700'}
gdfClip.crs = {'init':'epsg:27700'}

# Get the area to filter geometries by
ClipArea = gdfClip.loc[0,'geometry']

# Select only the polygons that intersect or lie within the junc clip area
gdfTopoArea = gdfTopoArea.loc[ (gdfTopoArea.geometry.intersects(ClipArea)) | (gdfTopoArea.geometry.within(ClipArea))]


# Could do some cleaning here - multipart to singlepart

# Now select just the polygons for pedestran or vehicle movement
gdfPedestrian = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin(['(1:Path)',
																'(2:Path,Structure)',
																'(2:Path,Tidal Water)',
																'(2:Roadside,Structure)',
																'(1:Roadside)'])]

gdfVehicle = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin([	'(1:Road Or Track)',
																'(2:Road Or Track,Structure)'])]


# Save these polygon layers
gdfPedestrian.to_file(os.path.join(output_directory, pedestrian_file))
gdfVehicle.to_file(os.path.join(output_directory, vehicle_file))