import pandas as pd
import numpy as np
import os
import geopandas as gpd 

##################
#
# Config
#
##################
gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
output_directory = os.path.join(gis_data_dir, "processed_gis_data")

ped_ods_file = os.path.join(output_directory, "OD_pedestrian_nodes.shp")
output_file = os.path.join(output_directory, "pedestrianODFlows.csv")


#################
#
# Generate Flows
#
#################

# Load ped ODs to get number of origins/destinations
gdf_od = gpd.read_file(ped_ods_file)

n_origins = gdf_od['id'].unique().shape[0]
n_destinations = n_origins

random_flows = np.random.rand(n_origins,n_destinations)
df1 = pd.DataFrame(random_flows)

# Set diagonal elements to 0, don't want an agent's destination to be the same as its origin
for i in range(max(n_origins,n_destinations)):
	df1.loc[i,i] = 0
df1.to_csv(output_file, index=False, header=False)
