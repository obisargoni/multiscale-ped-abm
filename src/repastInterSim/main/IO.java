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

import repast.simphony.util.collections.IndexedIterable;

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

}

	



