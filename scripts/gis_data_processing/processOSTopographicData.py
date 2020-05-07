import os
import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString


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

qgis_workings_dir = os.path.join(gis_data_dir, "qgis_workings")

if os.path.isdir(qgis_workings_dir) == False:
    os.mkdir(qgis_workings_dir)

output_directory = os.path.join(gis_data_dir, "processed_gis_data")

if os.path.isdir(output_directory) == False:
    os.mkdir(output_directory)


selection_layer_file = os.path.join(gis_data_dir, "JunctClip.shp")

topographic_area_file = os.path.join(topographic_area_dir, 'mastermap TopographicArea.shp')
topographic_line_file = os.path.join(topographic_line_dir, 'mastermap-topo_2903032_0 TopographicLine.shp')

output_vehicle_file = os.path.join(output_directory, "topographicAreaVehicle.shp")
output_pedestrian_file = os.path.join(output_directory, "topographicAreaPedestrian.shp")
output_line_file = os.path.join(output_directory, "boundaryPedestrianVehicleArea.shp")

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



################################
#
# Process topographic area data
#
################################

# Select only the polygons that intersect or lie within the junc clip area
gdfTopoAreaFiltered = gdfTopoArea.loc[ (gdfTopoArea.geometry.intersects(SelectPolygon)) | (gdfTopoArea.geometry.within(SelectPolygon))]


# Select vehicle areas as those that intersect the road network
# Results in 252 polygons, many more that previous method
gdfVehicle = gpd.sjoin(gdfTopoArea, gdfITNLink.loc[:,['geometry']], op = 'intersects', rsuffix = 'itn')
gdfVehicle.drop('index_itn', axis = 1, inplace=True)
gdfVehicle.drop_duplicates(inplace = True)


#gdfVehicle = gdfTopoArea.loc[ (gdfTopoArea.geometry.intersects(gdfITNLink.geometry)) | (gdfTopoArea.geometry.intersects(gdfITNNode.geometry)) ]

# Select polygons with descriptive types that correspons to areas where pedestrians typically have priority
gdfPedestrian = gdfTopoArea.loc[gdfTopoArea['descriptiv'].isin(['(1:Path)',
                                                                '(2:Path,Structure)',
                                                                '(2:Path,Tidal Water)',
                                                                '(2:Roadside,Structure)',
                                                                '(1:Roadside)'])]

# Filter pedestrian polygons to just be those within the study area or that touch a vehicle polygon
gdfPedestrianA = gdfPedestrian.loc[ (gdfPedestrian.geometry.intersects(SelectPolygon)) | (gdfPedestrian.geometry.within(SelectPolygon))]

gdfPedestrianB = gdfPedestrian.copy()

# Loop through geometries and include them if they touch a vehicle polygon
keep_indices = []
for i, row in gdfPedestrianB.iterrows():
    poly = row['geometry']
    for veh_poly in gdfVehicle['geometry']:
        if poly.touches(veh_poly):
            keep_indices.append(i)
            break

gdfPedestrianB = gdfPedestrianB.loc[keep_indices]

# Overwrite initial ped polygon gdf with the polygons that are within the study area
gdfPedestrian = pd.concat([gdfPedestrianA, gdfPedestrianB]).drop_duplicates()

gdfPedestrianA = gdfPedestrianB = None 


# Add in priority field
gdfPedestrian[priority_column] = "pedestrian"
gdfVehicle[priority_column] = "vehicle"


########################################
#
# Select only polygons that are within the largest connected component, avoids inaccessible islands
#
########################################

# Combine pedestrian and vehicle areas into a single geo series
gdfPedVeh = pd.concat([gdfVehicle, gdfPedestrian])



# initialise df to recod the neighbours for each ped-vehicle space polygon. use to identify clusters
dfPedVehNeighbours = pd.DataFrame()
for index, row in gdfPedVeh.iterrows():
    # Touches identifies polygons with at least one point in common but interiors don't intersect. So will work as long as none of my topographic polygons intersect
    neighborFIDs = gdfPedVeh[gdfPedVeh.geometry.touches(row['geometry'])]['fid'].tolist() 
    #neighborFIDs = neighborFIDs.remove(row['fid']) polygons don't touch them selves because they intersect
    df = pd.DataFrame({'neighbourfid':neighborFIDs})
    df['fid'] = row['fid']
    dfPedVehNeighbours = pd.concat([dfPedVehNeighbours, df])


g = nx.from_pandas_edgelist(dfPedVehNeighbours, source='fid', target = 'neighbourfid')
ccs = list(nx.connected_components(g))
ccs.sort(key=len)

cc_largest = ccs.pop(-1)


# Ideally want to restrict analysis to just the largest connected component since other ccs are likely to be court yard or other areas separate from roads
# For the smaller ccs check that they dont intersect with road network, in which case they can be excluded
# (might also want to check the road network is a single connected component)

check_ok = True
for cc in ccs:
    gdfCC = gdfPedVeh.loc[ gdfPedVeh['fid'].isin(cc)]

    # Check that none of these intersect with road links
    if gpd.sjoin(gdfCC, gdfITNLink, op='intersects').shape[0] != 0:
        print("Connected component at index {} intersects the road network".format(ccs.index(cc)))
        check_ok = False

assert check_ok


# Filter pedestrian and vehicle polygons to just include those in the largest connected component
gdfPedVeh = gdfPedVeh.loc[ gdfPedVeh['fid'].isin(cc_largest)]

# Overwrite previous pedestrian and vehicle polygon geo dataframes as will want to replace witht he filtered data
gdfPedestrian = gdfVehicle = None

# Could do some cleaning here - multipart to singlepart


gdfPedestrian = gdfPedVeh.loc[ gdfPedVeh[priority_column] == 'pedestrian']
gdfVehicle = gdfPedVeh.loc[ gdfPedVeh[priority_column] == 'vehicle']


##################################
#
# Process topographic line data
#
#
# Have decided to use the perimiter of the pedestrian and vehicle polygons combined. This is less restrictive that the topographic lines
# that are categorised as obstructing, which in some places wouls unrealistically constrict pedestrian movement.
#
##################################

'''
# Read in the data
gdfTopoLine = gpd.read_file(topographic_line_file)

# Set the crs
gdfTopoLine.crs = projectCRS

# Now select just the polygons for pedestran or vehicle movement
gdfTopoLine = gdfTopoLine.loc[gdfTopoLine['physicalPr'] == "Obstructing"]

# Select only the Obstructing lines that boarder the pedestrian or vehicle areas
gdfTopoLineSJ = gpd.sjoin(gdfTopoLine, gdfPedVeh, how = 'inner', op='intersects')
gdfTopoLineSJ.drop_duplicates(inplace=True)
'''

# Get perimiter(s) of pedestrian and vehicle areas and include this with the obstructing lines (will need to make sure that this perimiter includes all ITN network)

# Disolve pedestrian and vehicle polygons into single polygon and extract perimiter
gdfPedVeh['dissolve_key'] = 1
gdfDissolved = gdfPedVeh.dissolve(by = "dissolve_key")

# Disolved geometries are multi polygons, explode to single polygons
def explode(indf, single_type = Polygon, multi_type = MultiPolygon):
    outdf = gpd.GeoDataFrame(columns=indf.columns)
    for idx, row in indf.iterrows():
        if type(row.geometry) == single_type:
            outdf = outdf.append(row,ignore_index=True)
        if type(row.geometry) == multi_type:
            multdf = gpd.GeoDataFrame(columns=indf.columns)
            recs = len(row.geometry)
            multdf = multdf.append([row]*recs,ignore_index=True)
            for geom in range(recs):
                multdf.loc[geom,'geometry'] = row.geometry[geom]
            outdf = outdf.append(multdf,ignore_index=True)
    return outdf

gdfDissolved = explode(gdfDissolved)

# Get linstrings of exterior and interior of the dissolved pedestrian + vehicle polygons. These will mark the perimiters of the space
gdfPerimiter = gpd.GeoDataFrame(columns = ['type','geometry'])
gdfPerimiter.crs = gdfDissolved.crs
for geom in gdfDissolved['geometry']:
    exterior = LineString(geom.exterior.coords)
    gdfPerimiter = gdfPerimiter.append({"type":"exterior", "geometry":exterior}, ignore_index = True)

    for i in geom.interiors:
        interior = LineString(i)
        gdfPerimiter = gdfPerimiter.append({"type":"interior", "geometry":interior}, ignore_index = True)

gdfPerimiter.crs = projectCRS
gdfPerimiter[priority_column] = "pedestrian_obstruction"

# Check how many exterior boundaries - ideally want one by can have more than this is largest connected component of ped&veh polys includes a polygon
# that only touches at one coordinate
print(gdfPerimiter['type'].value_counts())


###########################
#
#
# Linking Pedestrian and Vehicle polygons with Road Link IDs
#
#
###########################

# Buffer road node
node_buffer_dist = 1
gdfITNNodeBuffer = gdfITNNode.copy()
gdfITNNodeBuffer['geometry'] = gdfITNNode.buffer(node_buffer_dist)

# Clip road links to exclude the junctions
gdfITNLinkClipped = gpd.overlay(gdfITNLink, gdfITNNodeBuffer, how = 'difference')

# Save to file so can be read into QGIS layer
output_itn_link_clipped_file = os.path.join(qgis_workings_dir, "itn_link_clipped.shp")
gdfITNLinkClipped.to_file(output_itn_link_clipped_file)


# Initialise QGIS API
import sys
from qgis.core import *

from qgis.analysis import QgsNativeAlgorithms
from qgis.PyQt.QtCore import QVariant

# Initialise QGIS
qgis_dir = "C:\\Program Files\\QGIS 3.4\\"
QgsApplication.setPrefixPath(os.path.join(qgis_dir, "bin\\qgis-ltr-bin-g7.exe"), True)
qgs = QgsApplication([], False)
qgs.initQgis()

# Append the path where processing plugin can be found
sys.path.append(os.path.join(qgis_dir, '\\apps\\qgis-ltr\\python\\plugins'))


from qgis import processing
from processing.core.Processing import Processing
Processing.initialize()
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())



itn_links_clipped = QgsVectorLayer(output_itn_link_clipped_file, "itn_links_clipped", "ogr")

'''
# Loop through features and densify
densify_distance = 5 # Points separated by 5m
for f in itn_links_clipped.getFeatures():
    g = f.geometry()
    newG = g.densifyByDistance(1)
    f.setGeometry(newG)
'''

# Extract points
points_path = os.path.join(qgis_workings_dir, "itn_links_dense_clipped_points.shp")
fields = itn_links_clipped.fields()
fields.append(QgsField("point_id",  QVariant.Int))
writer = QgsVectorFileWriter(points_path, "UTF-8", fields, QgsWkbTypes.Point, driverName="ESRI Shapefile")

if writer.hasError() != QgsVectorFileWriter.NoError:
	print("Error when creating shapefile: ",  writer.errorMessage())

# Loop through the densified cliped linestrings and extract the points
point_id = 0
for f in itn_links_clipped.getFeatures():
    g = f.geometry()
    gD = g.densifyByDistance(1)
    for p in gD.vertices():
        output_feature = QgsFeature()
        attrs = f.attributes()
        attrs.append(point_id)
        output_feature.setAttributes(attrs)
        output_feature.setGeometry(p)
        writer.addFeature(output_feature)
del writer

itn_links_dense_clipped_points = QgsVectorLayer(points_path, "itn_links_dense_clipped_points", "ogr")

# Now get voronoi polygons from these points
outpath = os.path.join(qgis_workings_dir, "itn_links_voronoi.shp")
result = processing.run("qgis:voronoipolygons", {"INPUT":itn_links_dense_clipped_points, "OUTPUT":outpath})

# Read in voronoi polygons as geodataframe
gdfLinksVoronoi = gpd.read_file(outpath)
gdfLinksVoronoi.crs = projectCRS

# Disolve by road link fid to get road link voroni regions
gdfLinksVoronoi = gdfLinksVoronoi.dissolve(by = 'fid')
gdfLinksVoronoi.reset_index(inplace=True)
gdfLinksVoronoi.to_file(os.path.join(qgis_workings_dir, "itn_links_voronoi_dissolve.shp"))


# Intersect pedestrian and vehicle priority polygons with links voronoi and disolve by road link fid
gdfPedestrianLinks = gpd.overlay(gdfPedestrian, gdfLinksVoronoi, how = 'intersection')

# Duplicate columns result in suffixes, edit these to make more sense
rename_dict = {i:i.replace('_1','_topo').replace('_2','_itn') for i in gdfPedestrianLinks.columns}
gdfPedestrianLinks.rename(columns = rename_dict, inplace = True)
gdfPedestrianLinks = gdfPedestrianLinks.dissolve(by = 'fid_itn').reset_index()
gdfPedestrianLinks = explode(gdfPedestrianLinks)
gdfPedestrianLinks.crs = gdfPedestrian.crs

gdfVehicleLinks = gpd.overlay(gdfVehicle, gdfLinksVoronoi, how = 'intersection')
rename_dict = {i:i.replace('_1','_topo').replace('_2','_itn') for i in gdfVehicleLinks.columns}
gdfVehicleLinks.rename(columns = rename_dict, inplace = True)
gdfVehicleLinks = gdfVehicleLinks.dissolve(by = 'fid_itn').reset_index()
gdfVehicleLinks = explode(gdfVehicleLinks)
gdfVehicleLinks.crs = gdfVehicle.crs

# Rename the ITN Road Link FID column to match the name expected by the Repast Model
gdfVehicleLinks.rename(columns = {"fid_itn":"roadLinkID"}, inplace = True)
gdfPedestrianLinks.rename(columns = {"fid_itn":"roadLinkID"}, inplace = True)

###########################
#
# Save processed data
#
###########################
# Save these polygon layers
gdfPedestrianLinks.to_file(output_pedestrian_file)
gdfVehicleLinks.to_file( output_vehicle_file)
gdfPerimiter.to_file(output_line_file)

# Save the selected ITN geometries
gdfITNLink.to_file(output_itn_link_file)
gdfITNNode.to_file(output_itn_node_file)

'''
# Now densify the clipped linestrings
def densify(l, d):
    length = l.length

    if l.length < d:
        return l

    # Find minimum number of points needed between start and end to achieve max separation of d
    np = int(length / d)
    dp = length / (np + 1) # total length divided by number of gaps

    points = []
    points.append(l.coords[0])
    dtot = dp

    # number of points to add in is np
    for i in range(np):
        points.append(l.interpolate(dp*i))

    # Add final point
    points.append(l.coords[-1])
    
    # Return densified linestring
    return LineString(points)


'''