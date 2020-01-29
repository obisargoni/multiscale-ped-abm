import numpy as np

oldFloodFillValuesFile = "..\\data\\export\\post-summand-refactor\\export_flood_fill_values.csv"
newFloodFillValuesFile = "..\\data\\export\\export_flood_fill_values.csv"

oldFloodFillValues = np.genfromtxt(oldFloodFillValuesFile, delimiter=",")
newFloodFillValues = np.genfromtxt(newFloodFillValuesFile, delimiter=",")

# 99.3% the same
# Check difference in values
diff = oldFloodFillValues - newFloodFillValues

print(np.nanmax(diff))
print(np.nanmin(diff))