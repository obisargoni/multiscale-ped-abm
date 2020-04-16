# IPython log file

import os
import geopandas as gpd
import pandas as pd
#from operator import itemgetter

data_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\mastermap-topo_2903032\\mastermap-topo_2903032_0 TopographicLine"
fileTopoLine = 'mastermap-topo_2903032_0 TopographicLine.shp'

output_directory = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\Intersect Within layers"
output_file = "topographicLineObstructing_VehiclePedestrianIntersect.shp"

pedestrian_area_file = "topographicAreaPedestrian.shp"
vehicle_area_file = "topographicAreaVehicle.shp"

# Read in the data
gdfTopoLine = gpd.read_file(os.path.join(data_directory,fileTopoLine))
#gdfClip = gpd.read_file(clip_layer_file)
gdfPedestrian = gpd.read_file(os.path.join(output_directory, pedestrian_area_file))
gdfVehicle = gpd.read_file(os.path.join(output_directory, vehicle_area_file))

# Set the crs
gdfTopoLine.crs = {'init' :'epsg:27700'}
gdfPedestrian.crs = {'init' :'epsg:27700'}
gdfVehicle.crs = {'init' :'epsg:27700'}

# Combine pedestrian and vehicle areas into a single geo series
gsPedVeh = pd.concat([gdfVehicle['geometry'], gdfPedestrian['geometry']])
gdfPedVeh = gpd.GeoDataFrame({'geometry':gsPedVeh})
gdfPedVeh.crs = {'init':'epsg:27700'}

# Now select just the polygons for pedestran or vehicle movement
gdfTopoLine = gdfTopoLine.loc[gdfTopoLine['physicalPr'] == "Obstructing"]

# Select only the Obstructing lines that boarder the pedestrian or vehicle areas
gdfTopoLineSJ = gpd.sjoin(gdfTopoLine, gdfPedVeh, how = 'inner', op='intersects')

# Could do some cleaning here - multipart to singlepart

# Save these polygon layers
gdfTopoLineSJ.to_file(os.path.join(output_directory, output_file))