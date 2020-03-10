# Python script for visualising arrays of data produced by the repast simulation
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import csv
import os

def grid_mask_array(grid, mask_coordiantes):
	array = np.empty(grid.shape)
	for coord in mask_coordiantes:
		c = tuple(reversed(coord.astype(int)))
		array[c] = 1
	maskArray = np.ma.masked_not_equal(array, 1)
	return maskArray

def plot_grid_image(array, colour_map, alpha = 1, mask_coordiantes = None, mask_cm = None, output_path = None):
	fig, ax = plt.subplots(figsize = (15,15))
	im = ax.imshow(array, cmap=colour_map, alpha = alpha, origin='lower', extent=[0, array.shape[1], 0, array.shape[0]])
	plt.colorbar(im)

	if(mask_coordiantes is not None):
		mask = grid_mask_array(array, mask_coordiantes)
		im = ax.imshow(mask, mask_cm)

	if output_path is not None:
		plt.savefig(output_path)
	else:
		plt.show()

	return fig, ax

output_dir = "..\\output\\"
export_dir = os.path.join(output_dir,"export\\")

input_files = {	"grid_path":"export_grid_coverage_path.csv",
				"pruned_path":"export_pruned_grid_coverage_path.csv",
				"path_crossings":"export_grid_coverage_path_crossings.csv"}

output_files = {"grid_path":"floodValues_path.png",
				"pruned_path":"gridValues_pruned.png",
				"path_crossings":"gridValues_crossing.png"}

# Input paths
gridValuesFile = os.path.join(export_dir, "export_grid_coverage_values.csv")
floodFillValuesFile = os.path.join(export_dir,"export_flood_fill_values.csv")
# Output paths
gridOutput = os.path.join(output_dir,"gridValues.png")
floodOutput = os.path.join(output_dir,"floodValues.png")


# Data input and cleaning
gridValues = np.genfromtxt(gridValuesFile, delimiter=',')
#floodValues = np.genfromtxt(floodFillValuesFile, delimiter = ',')

if(np.isnan(gridValues[:,-1]).all()):
	gridValues = gridValues[:,:-1]
'''
if(np.isnan(floodValues[:,-1]).all()):
	floodValues = floodValues[:,:-1]

max_int = 2147483647
floodValues = np.where(floodValues == max_int, np.nan, floodValues)
'''

# Plot grid values and flood fill
plot_grid_image(gridValues, cm.Blues, output_path = gridOutput)
#plot_grid_image(floodValues, cm.Blues, output_path = floodOutput)


for prefix in ["initial_", "final_"]:
	in_path_file = prefix + input_files['grid_path']
	in_prune_file = prefix + input_files['pruned_path']
	in_crossing_file = prefix + input_files['path_crossings']

	out_path_file = prefix + output_files['grid_path']
	out_prune_file = prefix + output_files['pruned_path']
	out_crossing_file = prefix + output_files['path_crossings']

	gridPathFile = os.path.join(export_dir,in_path_file)
	prunedPathFile = os.path.join(export_dir,in_prune_file)
	pathCrossingsFile = os.path.join(export_dir,in_crossing_file)

	floodPathOutput = os.path.join(output_dir,out_path_file)
	gridPrunedOutput = os.path.join(output_dir,out_prune_file)
	gridCrossingOutput = os.path.join(output_dir,out_crossing_file)

	gridPath = np.genfromtxt(gridPathFile, delimiter = ',')
	prunedPath = np.genfromtxt(prunedPathFile, delimiter = ',')
	pathCrossings = np.genfromtxt(pathCrossingsFile, delimiter = ',')

	# Make sure pruned path is nested coords, required if only one coord in pruned path
	if len(prunedPath.shape) == 1:
		prunedPath = np.array([prunedPath])

	# Plot flood with with path
	plot_grid_image(gridValues, cm.Blues, 0.5, gridPath, 'cool_r', output_path = floodPathOutput)
	# Plot grid values with pruned route and crossing points
	plot_grid_image(gridValues, cm.Blues, 0.5, prunedPath, 'cool_r', output_path = gridPrunedOutput)
	plot_grid_image(gridValues, cm.Blues, 0.5, pathCrossings, 'cool_r', output_path = gridCrossingOutput)