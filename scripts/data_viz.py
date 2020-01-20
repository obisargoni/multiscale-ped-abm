# Python script for visualising arrays of data produced by the repast simulation
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np
import csv

gridValuesFile = "..\\data\\output_grid_coverage_values.csv"
floodFillValuesFile = "..\\data\\output_flood_fill_values.csv"

gridValues = np.genfromtxt(gridValuesFile, delimiter=',')

# Create 2d arrays to represent all possible dealer and player hands
x = np.arange(gridValues[0].size) * np.ones(gridValues.size)[:, np.newaxis]
y = np.ones([gridValues.size,gridValues[0].size]) * np.arange(gridValues.size)[:, np.newaxis]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(x, y, z, rstride=1, cstride=1)
surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,linewidth=0, antialiased=False)
plt.show()
