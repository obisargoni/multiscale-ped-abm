import os
import geopandas as gpd
import pandas as pd


################################
#
# Initalise data directories and file names
#
###############################

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
topographic_data_dir = os.path.join(gis_data_dir, "mastermap-topo_2903032\\mastermap-topo_2903032_0 TopographicArea")
output_directory = os.path.join(gis_data_dir, "Intersect Within layers")

if os.path.isdir(output_directory) == False:
	os.mkdirs(output_directory)


selection_layer_file = os.path.join(gis_data_dir, "JunctClip.shp")
topographic_data_file = os.path.join(topographic_data_dir, 'mastermap TopographicArea.shp')
output_vehicle_file = os.path.join(output_directory, "topographicAreaVehicle.shp")
output_pedestrian_file = os.path.join(output_directory, "topographicAreaPedestrian.shp")

projectCRS = {'init' :'epsg:27700'}

priority_column = "priority"

#################################
#
# Read in the data
#
#################################
gdfTopoArea = gpd.read_file(topographic_data_file)
gdfSelect = gpd.read_file(selection_layer_file)

# Set the crs
gdfTopoArea.crs = projectCRS
gdfSelect.crs = projectCRS

# Get the area to filter geometries by
SelectPolygon = gdfSelect.loc[0,'geometry']


################################
#
# Select just the pologyons of interest
#
################################

# Select only the polygons that intersect or lie within the junc clip area
gdfTopoArea = gdfTopoArea.loc[ (gdfTopoArea.geometry.intersects(SelectPolygon)) | (gdfTopoArea.geometry.within(SelectPolygon))]


# Could do some cleaning here - multipart to singlepart

# Now select just the polygons for pedestran or vehicle movement
gdfPedestrian = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin(['(1:Path)',
																'(2:Path,Structure)',
																'(2:Path,Tidal Water)',
																'(2:Roadside,Structure)',
																'(1:Roadside)'])]

gdfVehicle = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin([	'(1:Road Or Track)',
																'(2:Road Or Track,Structure)'])]



###########################
#
# Save processed data
#
###########################
# Save these polygon layers
gdfPedestrian.to_file(os.path.join(output_directory, output_pedestrian_file))
gdfVehicle.to_file(os.path.join(output_directory, output_vehicle_file))