# Python script for visualising arrays of data produced by the repast simulation
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import csv

# Input paths
gridValuesFile = "..\\data\\output_grid_coverage_values.csv"
floodFillValuesFile = "..\\data\\output_flood_fill_values.csv"

# Output paths
gridOutput = "..\\output\\gridValues.png"
floodOutput = "..\\output\\floodValues.png"

# Data input and cleaning
gridValues = np.genfromtxt(gridValuesFile, delimiter=',')
floodValues = np.genfromtxt(floodFillValuesFile, delimiter = ',')

assert np.isnan(gridValues[:,-1]).all()
assert np.isnan(floodValues[:,-1]).all()

gridValues = gridValues[:,:-1]
floodValues = floodValues[:,:-1]

max_int = 2147483647
floodValues = np.where(floodValues == max_int, np.nan, floodValues)

def plot_grid_image(array, colour_map, output_path = None):
	fig, ax = plt.subplots()
	im = ax.imshow(array, cmap=colour_map,origin='lower', extent=[0, array.shape[1], 0, array.shape[0]])
	plt.colorbar(im)
	if output_path is None:
		plt.show()
	else:
		plt.savefig(output_path, bbox_inches='tight')

	return fig, ax

plot_grid_image(gridValues, cm.YlOrRd, gridOutput)
plot_grid_image(floodValues, cm.YlOrRd, floodOutput)