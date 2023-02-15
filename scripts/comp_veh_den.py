import json

import os

import pandas as pd

import numpy as np

import geopandas as gpd

import batch_data_utils as bd_utils


#####################################

#

#

# Load config file

#

#

#####################################

with open(".//config.json") as f:

    config = json.load(f)

file_datetime_string = config['file_datetime_string']

setting = config['setting']

ug_results = '2023.Feb.15.08_28_59'
qg_results = '2023.Feb.15.08_49_00'
cc_results = '2023.Feb.15.09_10_11'

v_param = 'avNVehicles'

'''
# old results, using t_veh
ug_results = '2022.Nov.10.18_17_11'
qg_results = '2022.Nov.10.16_45_50'
cc_results = '2022.Nov.10.15_21_09'

v_param = 'addVehicleTicks'
'''


data_dir = config['batch_data_dir']
data_dir = "C:\\Users\\obisargoni\\eclipse-workspace\\repastInterSim\\output\\batch\\model_run_data"
img_dir = "..\\output\\img\\"

env_col = 'environment'
palette = ['#1b9e77', '#d95f02', '#7570b3']


cc_gis_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\Clapham Common\\processed_gis_data\\"
qg_gis_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\QuadGrid145NodesBuffer\\processed_gis_data\\"
ug_gis_dir = "S:\\CASA_obits_ucfnoth\\1. PhD Work\\GIS Data\\ToyGrid141\\processed_gis_data\\"

'''
pavement_links_file = os.path.join(gis_data_dir, config['pavement_links_file'])
pavement_nodes_file = os.path.join(gis_data_dir, config['pavement_nodes_file'])
or_links_file = os.path.join(gis_data_dir, config['openroads_link_processed_file'])
or_nodes_file = os.path.join(gis_data_dir, config['openroads_node_processed_file'])
itn_links_file = os.path.join(gis_data_dir, config['mastermap_itn_processed_direction_file'])
itn_nodes_file = os.path.join(gis_data_dir, config['mastermap_node_processed_file'])
crossing_alternatives_file = os.path.join(gis_data_dir, config['crossing_alternatives_file'])
ped_ods_file = os.path.join(gis_data_dir, config['pedestrian_od_file'])
'''

bin_dist=2

cc_data_paths = bd_utils.get_data_paths(cc_results, data_dir)
batch_file = cc_data_paths["batch_file"]
cc_av_vehicle_counts_file = cc_data_paths['av_vehicle_counts']

qg_data_paths = bd_utils.get_data_paths(qg_results, data_dir)
qg_av_vehicle_counts_file = qg_data_paths['av_vehicle_counts']
ug_data_paths = bd_utils.get_data_paths(ug_results, data_dir)
ug_av_vehicle_counts_file = ug_data_paths['av_vehicle_counts']


ug_output_paths = bd_utils.get_ouput_paths(ug_results, data_dir, nbins = bin_dist)
ug_output_path = ug_output_paths["output_route_data"]


qg_output_paths = bd_utils.get_ouput_paths(qg_results, data_dir, nbins = bin_dist)
qg_output_path = qg_output_paths["output_route_data"]


cc_output_paths = bd_utils.get_ouput_paths(cc_results, data_dir, nbins = bin_dist)
cc_output_path = cc_output_paths["output_route_data"]




# load road networks and calculate vehicle capacity

gdfccitn = gpd.read_file(os.path.join(cc_gis_dir, config['mastermap_itn_processed_direction_file']))
gdfqgitn = gpd.read_file(os.path.join(qg_gis_dir, config['mastermap_itn_processed_direction_file']))
gdfugitn = gpd.read_file(os.path.join(ug_gis_dir, config['mastermap_itn_processed_direction_file']))

veh_length = 3.5
for gdf in [gdfccitn, gdfqgitn, gdfugitn]:
    gdf['cap'] = gdf['geometry'].map(lambda x: max(1, (x.length/veh_length)//1))


# calculate vehicle density
dfVehCountsCC = bd_utils.get_road_link_vehicle_density_from_vehicle_counts(gdfccitn, cc_av_vehicle_counts_file, cc_output_paths['output_vehicle_density_file'])
dfVehCountsCC.rename(columns = {'fid':'itn_fid'}, inplace=True)

dfVehCountsQG = bd_utils.get_road_link_vehicle_density_from_vehicle_counts(gdfqgitn, qg_av_vehicle_counts_file, qg_output_paths['output_vehicle_density_file'])
dfVehCountsQG.rename(columns = {'fid':'itn_fid'}, inplace=True)

dfVehCountsUG = bd_utils.get_road_link_vehicle_density_from_vehicle_counts(gdfugitn, ug_av_vehicle_counts_file, ug_output_paths['output_vehicle_density_file'])
dfVehCountsUG.rename(columns = {'fid':'itn_fid'}, inplace=True)


dfccv = pd.read_csv(cc_output_paths['output_vehicle_density_file'])
dfqgv = pd.read_csv(qg_output_paths['output_vehicle_density_file'])
dfugv = pd.read_csv(ug_output_paths['output_vehicle_density_file'])
dfRun = pd.read_csv(os.path.join(data_dir, batch_file))
dfccv = pd.merge(dfRun, dfccv, on='run')
dfqgv = pd.merge(dfRun, dfqgv, on='run')
dfugv = pd.merge(dfRun, dfugv, on='run')


dfccv = pd.merge(dfccv, gdfccitn.reindex(columns = ['fid', 'cap']), on='fid')
dfqgv = pd.merge(dfqgv, gdfqgitn.reindex(columns = ['fid', 'cap']), on='fid')
dfugv = pd.merge(dfugv, gdfugitn.reindex(columns = ['fid', 'cap']), on='fid')

dfccv['cap_pcnt'] = (dfccv['AvVehCount'] / dfccv['cap']) * 100.00
dfqgv['cap_pcnt'] = (dfqgv['AvVehCount'] / dfqgv['cap']) * 100.00
dfugv['cap_pcnt'] = (dfugv['AvVehCount'] / dfugv['cap']) * 100.00


dfccagg = dfccv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].quantile([0.9,0.95,0.99]).unstack()
dfqgagg = dfqgv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].quantile([0.9,0.95,0.99]).unstack()
dfugagg = dfugv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].quantile([0.9,0.95,0.99]).unstack()

dfccaggd = dfccv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].describe()
dfqgaggd = dfqgv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].describe()
dfugaggd = dfugv.groupby([v_param,'addPedTicks','randomSeed'])['cap_pcnt'].describe()


dfccagg = pd.merge(dfccagg, dfccaggd, left_index=True, right_index=True).reset_index()
dfccagg = dfccagg.groupby([v_param,'addPedTicks']).mean().reset_index()
dfccagg['environment'] = 'CC'

dfqgagg = pd.merge(dfqgagg, dfqgaggd, left_index=True, right_index=True).reset_index()
dfqgagg = dfqgagg.groupby([v_param,'addPedTicks']).mean().reset_index()
dfqgagg['environment'] = 'QG'

dfugagg = pd.merge(dfugagg, dfugaggd, left_index=True, right_index=True).reset_index()
dfugagg = dfugagg.groupby([v_param,'addPedTicks']).mean().reset_index()
dfugagg['environment'] = 'UG'


dfagg = pd.concat([dfccagg, dfqgagg, dfugagg]).rename(columns = {0.9:'90%', 0.95:'95%', 0.99:'99%'})

dfagg.reindex(columns = ['environment',v_param,'addPedTicks','50%','75%','90%','95%','99%','max']).to_csv(os.path.join(data_dir, 'av_veh_den.csv'), index=False)

#dfccagg.to_csv(os.path.join(data_dir, 'cc_vehden_agg_{}.csv'.format(config['cc_results'])))
#dfqgagg.to_csv(os.path.join(data_dir, 'qg_vehden_agg_{}.csv'.format(config['qg_results'])))
#dfugagg.to_csv(os.path.join(data_dir, 'ug_vehden_agg_{}.csv'.format(config['ug_results'])))