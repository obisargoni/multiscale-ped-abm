package repastInterSim.main;

import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.WritableGridCoverage2D;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.agent.Ped;
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
    
    public static void exportFinalGridRouteData(Ped p, String filePrefix, Boolean idSuffix) {
    	if (filePrefix == null) {
    		filePrefix = "";
    	}
    	String suffix = "";
    	if (idSuffix) {
    		suffix = "_pedID_"+p.getID();
    	}
    	
    	String prunedGridPathFile = GlobalVars.exportDir + filePrefix +  "export_pruned_grid_coverage_path"+suffix+".csv";
    	String gridPathCrossingsFile = GlobalVars.exportDir + filePrefix +  "export_grid_coverage_path_crossings"+suffix+".csv";
    	String pedGridPathFile = GlobalVars.exportDir + filePrefix + "export_grid_coverage_path_final"+suffix+".csv";
    	
		// TODO Auto-generated method stub
    	// Because this method should be called once ped has completed their route,
    	// the grid path is the one actually used by the pedestrian, accounts for dynamic updates
		List<GridCoordinates2D> pedGridPath = p.getRoute().getGridPath();
		List<GridCoordinates2D> prunedGridPath = p.getRoute().getPrunedGridPath();
		List<GridCoordinates2D> gridPathCrossings = p.getRoute().getGridPathCrossings();
		
		IO.gridCoordiantesIterableToCSV(pedGridPath, pedGridPathFile);
		IO.gridCoordiantesIterableToCSV(prunedGridPath, prunedGridPathFile);
		IO.gridCoordiantesIterableToCSV(gridPathCrossings, gridPathCrossingsFile);
	}
}

	



