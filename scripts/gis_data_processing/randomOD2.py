import pandas as pd
import numpy as np
import os

##################
#
# Config
#
##################
n_origins = 4
n_destinations = 4

gis_data_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\CoventGardenWaterloo"
output_directory = os.path.join(gis_data_dir, "processed_gis_data")
output_file = os.path.join(output_directory, "pedestrianODFlows.csv")

#################
#
# Generate Flows
#
#################
random_flows = np.random.rand(n_origins,n_destinations)
df1 = pd.DataFrame(random_flows)

# Set diagonal elements to 0, don't want an agent's destination to be the same as its origin
for i in range(8):
	df1.loc[i,i] = 0
df1.to_csv(output_file, index=False, header=False)
