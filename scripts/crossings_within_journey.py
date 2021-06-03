# Script to analyse where along their journey pedestrians cross

import json
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
import networkx as nx
from datetime import datetime as dt

import batch_data_utils as bd_utils

#####################################
#
# Globals
#
#####################################
project_crs = {'init': 'epsg:27700'}

with open(".//gis_data_processing//config.json") as f:
    config = json.load(f)

gis_data_dir = os.path.abspath("..\\data\\model_gis_data")
data_dir = config['batch_data_dir']
img_dir = ".\\output\\img\\"
l_re = re.compile(r"(\d+\.\d+),\s(\d+\.\d+)")

pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])

# Model output data
file_datetime_string = "2021.May.28.17_57_01"
file_datetime  =dt.strptime(file_datetime_string, "%Y.%b.%d.%H_%M_%S")
file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = None)
ped_crossings_file = os.path.join(data_dir, bd_utils.most_recent_directory_file(data_dir, file_re))

file_re = bd_utils.get_file_regex("pedestrian_pave_link_crossings", file_datetime = None, suffix = 'batch_param_map')
batch_file = bd_utils.most_recent_directory_file(data_dir, file_re)

# Output paths
img_dir = "..\\output\\img\\"
crossings_in_journey_fig = os.path.join(img_dir, "cross_in_journey.png")


#####################################
#
#
# Functions
#
#
#####################################

def get_path_index_of_crossing(row, path_col = 'FullStrategicPathString', cross_link_col = 'pedRLID'):
    try:
        index = row[path_col].index(row[cross_link_col])
    except ValueError as e:
        index = None
    return index


#####################################
#
#
# Load Data
#
#
#####################################

# Data from model run
dfPedCrossings = pd.read_csv(ped_crossings_file)
dfRun = pd.read_csv(os.path.join(data_dir, batch_file))

gdfPaveNetwork = gpd.read_file(pavement_links_file)

#######################################
#
#
# Process data
#
#
########################################

# Filter to just include pavement links used for crossing
dfPedCrossings = dfPedCrossings.loc[ ~dfPedCrossings['CurrentPavementLinkID'].isnull()]
dfPedCrossings = dfPedCrossings.merge(gdfPaveNetwork, left_on = 'CurrentPavementLinkID', right_on = 'fid', how = 'left', indicator=True)
assert dfPedCrossings.loc[ dfPedCrossings['_merge']!='both'].shape[0]==0
dfPedCrossings.drop('_merge', axis=1, inplace=True)

# Filter columns and drop duplicates to retain only each crossing pavement link in the journey Extract the strategic route from the data
dfCrossingJourney = dfPedCrossings.reindex(columns = ['run','ID','CurrentPavementLinkID', 'ChosenCrossingTypeString', 'pedRLID', 'FullStrategicPathString']).drop_duplicates()
dfCrossingJourney = dfCrossingJourney.loc[ dfCrossingJourney['pedRLID'].notnull()]
dfCrossingJourney['FullStrategicPathString'] = dfCrossingJourney['FullStrategicPathString'].map(lambda s: s.strip('-').split('-'))

# Identify the strategic path link index crosses at and aggregate to get number of crossings made on 1st, 2nd, etc link
dfCrossingJourney['strategicCrossIndex'] = dfCrossingJourney.apply(get_path_index_of_crossing, axis=1)

# Check crossing index values
# 1. Where the road link crossed on is not in strategic path the chosen crossing type should be none. These are secondary crossings.
assert dfCrossingJourney.loc[(dfCrossingJourney['strategicCrossIndex'].isnull()), 'ChosenCrossingTypeString'].unique()[0] == 'none'

# 2. Otherwise there should be no instances where the index is null
# Fails. These are the cases where ped chooses not to cross, and therefore has to make secondary crossing the start of next tactical route.
#assert dfCrossingJourney.loc[ (dfCrossingJourney['ChosenCrossingTypeString']=='none') & (dfCrossingJourney['strategicCrossIndex'].notnull())].shape[0] == 0

# Remove secondary crossings
dfCrossingJourney = dfCrossingJourney.loc[ dfCrossingJourney['ChosenCrossingTypeString']!='none']

# Aggregate data
dfIndexCounts = dfCrossingJourney.groupby(['run', 'strategicCrossIndex']).apply(lambda g: g.shape[0]).reset_index()
dfIndexCounts.rename(columns = {0:'index_count'}, inplace=True)

# Join with run data
#gdfIndexCounts = pd.merge(gdfORLinks, dfIndexCounts, left_on = 'fid', right_on = 'pedRLID', how = 'left')
#gdfIndexCounts['cross_count'] = gdfIndexCounts['cross_count'].fillna(0)
#dfIndexCounts = pd.merge(dfIndexCounts, dfRun, on = 'run')
dfIndexCounts = dfIndexCounts.set_index(["run","strategicCrossIndex"]).unstack()
dfIndexCounts.columns = [c[1] for c in dfIndexCounts.columns]
dfIndexCounts['9+'] = dfIndexCounts.iloc[:, 9:].sum(axis=1)
dfIndexCounts = dfIndexCounts.iloc[:, [0,1,2,3,4,5,6,7,8,-1]]
dfIndexCounts.fillna(0, inplace=True)

##########################################
#
#
# Visualise
#
#
##########################################

from matplotlib import pyplot as plt

groups = ["Run {}".format(i) for i in dfIndexCounts.index]
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
labels = [ordinal(n+1) for n in range(len(dfIndexCounts.columns))]
x = np.arange(len(labels))  # the label locations
width = 0.8 / len(labels)  # the width of the bars, leave 0.2 spacing between groups

fig, ax = plt.subplots(figsize = (15,15))
rects = []
for ix, i in enumerate(dfIndexCounts.index):
    data = dfIndexCounts.loc[i].values
    data_rects = ax.bar(x + (ix * width), data, width, label=groups[ix])
    rects.append(data_rects)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of crossings')
ax.set_xlabel('Strategic Link Crossed On')
ax.set_title('Where along the journey peds cross')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

fig.tight_layout()

plt.show()
plt.savefig(crossings_in_journey_fig)