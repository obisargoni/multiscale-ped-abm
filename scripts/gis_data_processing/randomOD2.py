import pandas as pd
import numpy as np

##################
#
# Config
#
##################
n_origins = 4
n_destinations = 4
output_file = "pedestrianODFlows.csv"

#################
#
# Generate Flows
#
#################
random_flows = np.random.rand(n_origins,n_destinations)
df1 = pd.DataFrame(random_flows)
df1

# Set diagonal elements to 0, don't want an agent's destination to be the same as its origin
for i in range(8):
	df1.loc[i,i] = 0
df1.to_csv(output_file, index=False, header=False)
