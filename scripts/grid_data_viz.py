# Python script for visualising arrays of data produced by the repast simulation
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import csv

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



# Input paths
gridValuesFile = "..\\output\\export\\export_grid_coverage_values.csv"
floodFillValuesFile = "..\\output\\export\\export_flood_fill_values.csv"

gridPathFile = "..\\output\\export\\export_grid_coverage_path.csv"
gridPathSectionFile = "..\\output\\export\\export_grid_coverage_path_section.csv"
prunedPathFile = "..\\output\\export\\export_pruned_grid_coverage_path.csv"
pathCrossingsFile = "..\\output\\export\\export_grid_coverage_path_crossings.csv"


# Output paths
gridOutput = "..\\output\\gridValues.png"
floodOutput = "..\\output\\floodValues.png"
floodPathOutput = "..\\output\\floodValues_path.png"
floodPathSectionOutput = "..\\output\\floodValues_path_section.png"
gridPrunedOutput = "..\\output\\gridValues_pruned.png"
gridCrossingOutput = "..\\output\\gridValues_crossing.png"

# Data input and cleaning
gridValues = np.genfromtxt(gridValuesFile, delimiter=',')
#floodValues = np.genfromtxt(floodFillValuesFile, delimiter = ',')

gridPath = np.genfromtxt(gridPathFile, delimiter = ',')
gridPathSection = np.genfromtxt(gridPathSectionFile, delimiter = ',')
prunedPath = np.genfromtxt(prunedPathFile, delimiter = ',')
pathCrossings = np.genfromtxt(pathCrossingsFile, delimiter = ',')

if(np.isnan(gridValues[:,-1]).all()):
	gridValues = gridValues[:,:-1]
'''
if(np.isnan(floodValues[:,-1]).all()):
	floodValues = floodValues[:,:-1]

max_int = 2147483647
floodValues = np.where(floodValues == max_int, np.nan, floodValues)
'''

# Make sure pruned path is nested coords, required if only one coord in pruned path
if len(prunedPath.shape) == 1:
	prunedPath = np.array([prunedPath])

# Plot grid values and flood fill
plot_grid_image(gridValues, cm.Blues, output_path = gridOutput)
#plot_grid_image(floodValues, cm.Blues, output_path = floodOutput)

# Plot flood with with path
plot_grid_image(gridValues, cm.Blues, 0.5, gridPath, 'cool_r', output_path = floodPathOutput)
plot_grid_image(gridValues, cm.Blues, 0.5, gridPathSection, 'cool_r', output_path = floodPathSectionOutput)

# Plot grid values with pruned route and crossing points
plot_grid_image(gridValues, cm.Blues, 0.5, prunedPath, 'cool_r', output_path = gridPrunedOutput)
plot_grid_image(gridValues, cm.Blues, 0.5, pathCrossings, 'cool_r', output_path = gridCrossingOutput)