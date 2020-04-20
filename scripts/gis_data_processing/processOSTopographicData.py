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
# Select only polygons that are withing the largest connected component, avoids unaccessible islands
#
########################################

# Combine pedestrian and vehicle areas into a single geo series
gdfPedVeh = pd.concat([gdfVehicle, gdfPedestrian])



# initialise df to recod the neighbours for each ped-vehicle space polygon. use to identify cluders
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
##################################


# Read in the data
gdfTopoLine = gpd.read_file(topographic_line_file)

# Set the crs
gdfTopoLine.crs = projectCRS

# Now select just the polygons for pedestran or vehicle movement
gdfTopoLine = gdfTopoLine.loc[gdfTopoLine['physicalPr'] == "Obstructing"]

# Select only the Obstructing lines that boarder the pedestrian or vehicle areas
gdfTopoLineSJ = gpd.sjoin(gdfTopoLine, gdfPedVeh, how = 'inner', op='intersects')
gdfTopoLineSJ.drop_duplicates(inplace=True)


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
gdfPerimiter.crs = projectCRS
for geom in gdfDissolved['geometry']:
    exterior = LineString(geom.exterior.coords)
    gdfPerimiter = gdfPerimiter.append({"type":"exterior", "geometry":exterior}, ignore_index = True)

    for i in geom.interiors:
        interior = LineString(i)
        gdfPerimiter = gdfPerimiter.append({"type":"interior", "geometry":interior}, ignore_index = True)

# Combine perimiter and obstruction linestrings
gdfLines = pd.concat([gdfTopoLineSJ, gdfPerimiter])

# Disolve together to avoid duplicated geometries
gdfLines['dissolve_key'] = 1
gdfLinesDissolved = gdfLines.dissolve(by = 'dissolve_key')
gdfLinesDissolved = explode(gdfLinesDissolved, single_type = LineString, multi_type = MultiLineString)

# Could do some cleaning here - multipart to singlepart - I think this has been taken care of by explode()
gdfLinesDissolved[priority_column] = "pedestrian_obstruction"




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
