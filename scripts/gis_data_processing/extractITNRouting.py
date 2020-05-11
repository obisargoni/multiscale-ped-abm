# Get lookup from RoadLinkID to the Plus and Minus Node IDs for each RoadLink in the ITN .gml file
# Get lookup from DirectedLinkID (same as RoadLinkID but for only directed links) to the orientaition (direction) of that road link

from lxml import etree
import geopandas as gpd
import pandas as pd
import numpy as np
import os

########################
#
# Initialise paths
#
########################
itnDir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo\\mastermap-itn_2903030\\"
itnF = "mastermap-itn_2903030_0.gml"

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
output_directory = os.path.join(gis_data_dir, "itn_route_info")

if os.path.isdir(output_directory) == False:
    os.mkdir(output_directory)

output_route_info_path = os.path.join(output_directory, "extracted_RRI.csv") 
output_node_info_path = os.path.join(output_directory, "extracted_RLNodes.csv") 

#######################
#
# Load Data
#
#######################
elemTree = etree.parse(os.path.join(itnDir, itnF))
itnTree = elemTree.getroot()
rl_elements = itnTree.findall(".//osgb:RoadLink", itnTree.nsmap)


dfRLNodes = pd.DataFrame(columns = ['RoadLinkFID', 'PlusNodeFID','MinusNodeFID'])

def get_road_link_nodes(road_link_element, namespaces, node_orientation = None):
    '''
    Returns a list of directed node elements. Optionally returns only those with a certain orientation.
    '''
    if node_orientation is None:
        xpath = ".//osgb:directedNode"
    elif node_orientation == '-' or node_orientation == '+':
        xpath = ".//osgb:directedNode[@orientation='{}']".format(node_orientation)
    else:
        return []
    return road_link_element.findall(xpath, namespaces)

for rl_element in rl_elements:
    rl_fid = rl_element.get('fid')
    plus_node_list = get_road_link_nodes(rl_element, itnTree.nsmap, node_orientation='+')
    assert len(plus_node_list) == 1
    minus_node_list = get_road_link_nodes(rl_element, itnTree.nsmap, node_orientation='-')
    assert len(minus_node_list) == 1
    plus_node_fid = plus_node_list[0].get('{}href'.format("{"+itnTree.nsmap['xlink']+"}")).replace("#","")
    minus_node_fid = minus_node_list[0].get('{}href'.format("{"+itnTree.nsmap['xlink']+"}")).replace("#","")
    dfRLNodes = dfRLNodes.append({'RoadLinkFID':rl_fid, 'PlusNodeFID':plus_node_fid, 'MinusNodeFID':minus_node_fid}, ignore_index = True)


# Check data frame for nulls

ri_elements = itnTree.findall(".//osgb:RoadRouteInformation", itnTree.nsmap)

dfRRI = pd.DataFrame(columns=["RoadRouteInformationFID", "DirectedLinkFID", "DirectedLinkOrientation"])

for ri_element in ri_elements:
	ri_fid = ri_element.get("fid")

	# Get all directed link elements
	dl_elements = ri_element.findall(".//osgb:directedLink", itnTree.nsmap)
	for dl_element in dl_elements:
		dl_fid = dl_element.get('{}href'.format("{"+itnTree.nsmap['xlink']+"}")).replace("#","") # Remove the preceeding hash
		dl_orientation = dl_element.get("orientation")

		# Add to the dataframe
		dfRRI = dfRRI.append({"RoadRouteInformationFID":ri_fid, "DirectedLinkFID":dl_fid, "DirectedLinkOrientation":dl_orientation}, ignore_index=True)

assert dfRLNodes.isnull().any().any() == False
assert dfRRI.isnull().any().any() == False


#########################
#
# Save data
#
#########################
dfRRI.to_csv(output_route_info_path, index= False)
dfRLNodes.to_csv(output_node_info_path, index = False)