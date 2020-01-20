package repastInterSim.main;

import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.geotools.coverage.grid.GridCoverage2D;

public class IO {
	
	public IO() {
		
	}
	
	public static void twodDoubleArrayToCSV(double [][] array, String path) {
		
		int width = array[0].length;
		int height = array.length;
		
		double arrayElem;
	    try {
			FileWriter writer = new FileWriter(path);			
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					arrayElem = array[i][j];
					writer.append(String.valueOf(arrayElem));
					
					// If not at end or row add delimiter
					if (i < width-1) {
						writer.append(",");
					}
					writer.append("\n");
				}
			writer.close();
			}
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

}

	



