package repastInterSim.main;

import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Enumeration;
import java.util.List;
import java.util.Properties;

import javax.imageio.ImageIO;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;

import com.opencsv.CSVReader;
import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.WritableGridCoverage2D;
import repastInterSim.environment.GISFunctions;

public class IO {
	
	public IO() {
		
	}
	
	public static void twodDoubleArrayToCSV(double [][] array, String path) {
		
		int width = array[0].length;
		int height = array.length;
		
		double arrayElem;
	    try {
			FileWriter writer = new FileWriter(path);			
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					arrayElem = array[j][i];
					writer.append(String.valueOf(arrayElem));
					
					// If at end of row don't add delimiter
					if (i != width-1) {
						writer.append(",");
					}
				}
				writer.append("\n");
			}
			writer.close();
	    }
			
	    catch (IOException e) {
		// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void gridCoordiantesIterableToCSV(Iterable<GridCoordinates2D> coordinates, String path) {
		
	    try {
			FileWriter writer = new FileWriter(path);			
			for (GridCoordinates2D c: coordinates) {
				boolean delim = false;
				for (int cValue: c.getCoordinateValues()) {
					if(delim) {
						writer.append(",");
					}
					writer.append(String.valueOf(cValue));
					delim = true;
				}
				writer.append("\n");
			}
			writer.close();
	    }
			
	    catch (IOException e) {
		// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void gridCoverageToImage(GridCoverage2D grid, String path) {
	    RenderedImage gridImage = grid.getRenderedImage();
	    File output = new File(path);
	    try {
			ImageIO.write(gridImage, "png", output);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void gridCoverageValuesToCSV(GridCoverage2D grid, String path) {
		int width = grid.getRenderedImage().getTileWidth();
		int height = grid.getRenderedImage().getTileHeight();
		double[] gridValue = null;
		double[][] gridValues = new double[height][width];			
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				GridCoordinates2D cell = new GridCoordinates2D(i,j);					
				gridValue = grid.evaluate(cell, gridValue);
				gridValues[j][i] = gridValue[0];
			}
		}
		
		twodDoubleArrayToCSV(gridValues, path);
	}
	
	public static void gridCoverageCoordinatesToCSV(GridCoverage2D grid, String path) {
		int width = grid.getRenderedImage().getTileWidth();
		int height = grid.getRenderedImage().getTileHeight();
		
		// 4 entries for each grid cell: grid cell coords, gis coords
		double[][] gridCoordiantes = new double[height*width][4];
		int rowi = 0;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				GridCoordinates2D cell = new GridCoordinates2D(i,j);
				Coordinate cellCoord = GISFunctions.gridCellToCoordinate(grid, cell);
				double[] row = {cell.x, cell.y, cellCoord.x, cellCoord.y};
				gridCoordiantes[rowi] = row;
				rowi++;
			}
		}
		twodDoubleArrayToCSV(gridCoordiantes, path);
	}
	
	public static void coordiantesIterableToCSV(Iterable<Coordinate> coordinates, String path) {
		String header = "x,y";
		
	    try {
			FileWriter writer = new FileWriter(path);
			writer.append(header);
			writer.append("\n");
			for (Coordinate c: coordinates) {
				writer.append(String.valueOf(c.x) + "," + String.valueOf(c.y));
				writer.append("\n");
			}
			writer.close();
	    }
			
	    catch (IOException e) {
		// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
    public static String getCoordinateListString(List<Coordinate> pedPrimaryRoute) {
    	String cString = "";
    	for (Coordinate c: pedPrimaryRoute) {
    		cString = cString + c.toString()+",";
    	}
    	return cString;
    }
    
    public static void exportGridCoverageData(WritableGridCoverage2D grid) {
    	String gridValueFile =  GlobalVars.exportDir + "export_grid_coverage_values.csv";
    	String gridImageFile = GlobalVars.outputDir + "output_grid_vales.png";
    	String griCoordinatesFile = GlobalVars.outputDir + "output_grid_coordinates.csv";
    			
		IO.gridCoverageToImage(grid, gridImageFile);
		IO.gridCoverageValuesToCSV(grid, gridValueFile);
		IO.gridCoverageCoordinatesToCSV(grid, griCoordinatesFile);
	}
    
    /*
    public static void exportPedGridRouteData(Ped p, String filePrefix, Boolean idSuffix) {
    	if (filePrefix == null) {
    		filePrefix = "";
    	}
    	String suffix = "";
    	if (idSuffix) {
    		suffix = "_pedID_"+p.getID();
    	}
    	
    	String prunedGridPathFile = GlobalVars.exportDir + filePrefix +  "export_pruned_grid_coverage_path"+suffix+".csv";
    	String gridPathCrossingsFile = GlobalVars.exportDir + filePrefix +  "export_grid_coverage_path_crossings"+suffix+".csv";
    	String gridPathFile = GlobalVars.exportDir + filePrefix + "export_grid_coverage_path"+suffix+".csv";
    	
		// TODO Auto-generated method stub
    	// Because this method should be called once ped has completed their route,
    	// the grid path is the one actually used by the pedestrian, accounts for dynamic updates
		List<GridCoordinates2D> gridPath = new ArrayList<GridCoordinates2D>();
		Map<GridCoordinates2D, List<GridCoordinates2D>> groupedPath = p.getRoute().getGroupedGridPath();
		for (GridCoordinates2D key: groupedPath.keySet()) {
			gridPath.addAll(groupedPath.get(key));
		}
		List<GridCoordinates2D> prunedGridPath = p.getRoute().getPrunedGridPath();
		List<GridCoordinates2D> gridPathCrossings = p.getRoute().getGridPathCrossings();
		
		IO.gridCoordiantesIterableToCSV(gridPath, gridPathFile);
		IO.gridCoordiantesIterableToCSV(prunedGridPath, prunedGridPathFile);
		IO.gridCoordiantesIterableToCSV(gridPathCrossings, gridPathCrossingsFile);
	}
	*/
    
	/*
	 * Read CSV data into List of String arrays. Catch exceptions.
	 * 
	 * @param String filePath. Path of the CSV file to read
	 * 
	 * @returns List<String[]> The csv data as a list of string arrays
	 */
	public static List<String[]> readCSV(String filePath) {
		// Read in OD matrix data for pedestrians from CSV
		List<String[]> csvData = null;
		CSVReader reader = null;
		try {
			reader = new CSVReader(new FileReader(filePath));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			csvData = reader.readAll();
			reader.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return csvData;
	}
	
	/**
	 * Get the value of a property in the properties file. If the input is empty or null or if there is no property with
	 * a matching name, throw a RuntimeException.
	 * 
	 * @param property
	 *            The property to look for.
	 * @return A value for the property with the given name.
	 */
	public static String getProperty(String property) {
		if (property == null || property.equals("")) {
			throw new RuntimeException("getProperty() error, input parameter (" + property + ") is "
					+ (property == null ? "null" : "empty"));
		} else {
			String val = GlobalVars.properties.getProperty(property);
			if (val == null || val.equals("")) { // No value exists in the
													// properties file
				throw new RuntimeException("checkProperty() error, the required property (" + property + ") is "
						+ (property == null ? "null" : "empty"));
			}
			return val;
		}
	}

	/**
	 * Read the properties file and add properties. Will check if any properties have been included on the command line
	 * as well as in the properties file, in these cases the entries in the properties file are ignored in preference
	 * for those specified on the command line.
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void readProperties() throws FileNotFoundException, IOException {

		File propFile = new File("./data/respastInterSim.properties");
		if (!propFile.exists()) {
			throw new FileNotFoundException("Could not find properties file in the default location: "
					+ propFile.getAbsolutePath());
		}

		//LOGGER.log(Level.FINE, "Initialising properties from file " + propFile.toString());

		GlobalVars.properties = new Properties();

		FileInputStream in = new FileInputStream(propFile.getAbsolutePath());
		GlobalVars.properties.load(in);
		in.close();

		// See if any properties are being overridden by command-line arguments
		for (Enumeration<?> e = GlobalVars.properties.propertyNames(); e.hasMoreElements();) {
			String k = (String) e.nextElement();
			String newVal = System.getProperty(k);
			if (newVal != null) {
				// The system property has the same name as the one from the
				// properties file, replace the one in the properties file.
				/*
				LOGGER.log(Level.INFO, "Found a system property '" + k + "->" + newVal
						+ "' which matches a NeissModel property '" + k + "->" + properties.getProperty(k)
						+ "', replacing the non-system one.");
						*/
				GlobalVars.properties.setProperty(k, newVal);
			}
		} // for
		return;
	} // readProperties
}

	



