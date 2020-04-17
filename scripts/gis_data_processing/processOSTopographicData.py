import os
import geopandas as gpd
import pandas as pd


################################
#
# Initalise data directories and file names
#
###############################

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
topographic_area_dir = os.path.join(gis_data_dir, "mastermap-topo_2903032\\mastermap-topo_2903032_0 TopographicArea")
topographic_line_dir = os.path.join(gis_data_dir, "mastermap-topo_2903032\\mastermap-topo_2903032_0 TopographicLine")


itn_directory = os.path.join(gis_data_dir, "mastermap-itn_2903030")
itn_link_file = os.path.join(itn_directory, "mastermap-itn RoadLink", "mastermap-itn 2903030_0 RoadLink.shp")
itn_node_file = os.path.join(itn_directory, "mastermap-itn RoadNode", "mastermap-itn RoadNode.shp")


output_directory = os.path.join(gis_data_dir, "processed_gis_data")

if os.path.isdir(output_directory) == False:
    os.mkdir(output_directory)


selection_layer_file = os.path.join(gis_data_dir, "JunctClip.shp")

topographic_area_file = os.path.join(topographic_area_dir, 'mastermap TopographicArea.shp')
topographic_line_file = os.path.join(topographic_line_dir, 'mastermap-topo_2903032_0 TopographicLine.shp')

output_vehicle_file = os.path.join(output_directory, "topographicAreaVehicle.shp")
output_pedestrian_file = os.path.join(output_directory, "topographicAreaPedestrian.shp")
output_line_file = os.path.join(output_directory, "topographicLineObstructing_VehiclePedestrianIntersect.shp")

output_itn_link_file = os.path.join(output_directory, "mastermap-itn RoadLink Intersect Within.shp")
output_itn_node_file = os.path.join(output_directory, "mastermap-itn RoadNode Intersect Within.shp")


projectCRS = {'init' :'epsg:27700'}

priority_column = "priority"


#################################
#
# Read in the data
#
#################################

# Mastermap topographic area data  - for pavement and carriageway polygons
gdfTopoArea = gpd.read_file(topographic_area_file)
gdfTopoArea.crs = projectCRS

# Mastermap ITN data - for road network
gdfITNLink = gpd.read_file(itn_link_file)
gdfITNLink.crs = projectCRS

gdfITNNode = gpd.read_file(itn_node_file)
gdfITNNode.crs = projectCRS

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
gdfITNNode.drop(['fid_line'], axis = 1, inplace=True)
gdfITNNode.rename(columns = {'fid_node':'fid'}, inplace = True)

# Could do some cleaning here - multipart to singlepart

# Now select just the polygons for pedestran or vehicle movement
gdfPedestrian = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin(['(1:Path)',
																'(2:Path,Structure)',
																'(2:Path,Tidal Water)',
																'(2:Roadside,Structure)',
																'(1:Roadside)'])]

gdfVehicle = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin([	'(1:Road Or Track)',
																'(2:Road Or Track,Structure)'])]


# Add in priority field
gdfPedestrian[priority] = "pedestrian"
gdfVehicle[priority] = "vehicle"

##################################
#
# Process topographic line data
#
##################################


# Read in the data
gdfTopoLine = gpd.read_file(topographic_line_file)

# Could do some cleaning here - multipart to singlepart


# Set the crs
gdfTopoLine.crs = projectCRS

# Combine pedestrian and vehicle areas into a single geo series
gsPedVeh = pd.concat([gdfVehicle['geometry'], gdfPedestrian['geometry']])
gdfPedVeh = gpd.GeoDataFrame({'geometry':gsPedVeh})

# Now select just the polygons for pedestran or vehicle movement
gdfTopoLine = gdfTopoLine.loc[gdfTopoLine['physicalPr'] == "Obstructing"]

# Select only the Obstructing lines that boarder the pedestrian or vehicle areas
gdfTopoLineSJ = gpd.sjoin(gdfTopoLine, gdfPedVeh, how = 'inner', op='intersects')


###########################
#
# Save processed data
#
###########################
# Save these polygon layers
gdfPedestrian.to_file(output_pedestrian_file)
gdfVehicle.to_file( output_vehicle_file)
gdfTopoLineSJ.to_file(output_line_file)# Save the selected ITN geometries
gdfITNLink.to_file(output_itn_link_file)
gdfITNNode.to_file(output_itn_node_file)
