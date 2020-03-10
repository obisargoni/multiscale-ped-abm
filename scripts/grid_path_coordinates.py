import pandas as pd
import numpy as np
import geopandas as gpd

# Read in the grid cell to coordinates lookup
df_gridcoords = pd.read_csv("..\\output\\output_grid_coordinates.csv", header = None)
df_gridcoords.columns = ['cellX','cellY','coordX','coordY']

# Filter to just select the coordinates of grid cells in the route
gridPathFile = "..\\output\\export\\export_grid_coverage_path.csv"
gridPath = np.genfromtxt(gridPathFile, delimiter = ',')
df_gridpath = pd.DataFrame(gridPath, columns = ['cellX','cellY'])

df_gridpath = pd.merge(df_gridpath, df_gridcoords, how = 'inner')
assert df_gridpath.isnull().any().any() == False

# Now convert to a geo data frame and save as a shapefile
gdf_path = gpd.GeoDataFrame(df_gridpath, geometry=gpd.points_from_xy(df_gridpath.coordX, df_gridpath.coordY))
gdf_path.crs = {'init': 'epsg:27700'}

gdf_path.to_file("..\\output\\path_coordinates.shp")